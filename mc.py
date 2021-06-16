from time import time as sec
from time import ctime
import numpy as np
import dxfgrabber as dxf
from pandas import DataFrame as df
from pandas import read_csv as rcsv
import core
from os import makedirs, chdir
from shutil import copyfile
import sys
from fdsafir import Thermal, user_config, Logger, progressBar
import fires

'''Read geometry and map it to the fire (choose the most exposed sections)'''


class Single:
    def __init__(self, title):
        self.title = title  # simulation title
        self.prof_type = 'invalid profile type'  # type of steel profile
        self.geometry = self.read_dxf()  # import geometry from DXF file

    # read dxf geometry
    def read_dxf(self):
        t1 = sec()

        print('Reading DXF geometry...', end='\r')
        dxffile = dxf.readfile('{}.dxf'.format(self.title))
        print('[OK] DXF geometry imported ({} s)'.format(round(sec()-t1, 2)))

        beams = []
        columns = []
        x = 0
        t = len(dxffile.entities)
        # assign LINES elements to columns or beams tables
        for ent in dxffile.entities:
            progressBar('Converting lines', x, t)
            if ent.dxftype == 'LINE':
                if ent.start[2] == ent.end[2]:
                    beams.append(ent)
                else:
                    columns.append(ent)
            x += 1
        print('[OK] Lines converted                         ')

        # assign 3DFACE elements to shells table
        shells = [ent for ent in dxffile.entities if ent.dxftype == '3DFACE']

        return {'b': beams, 'c': columns, 's': shells, 'f': []}

    # map fire and geometry - choose fire scenario
    def generate(self, f_coords, element):

        # select the most exposed section among the lines and return its config
        def map_lines():
            # no element exception
            lines = self.geometry['f']
            if lines.__len__() == 0:
                return (*f_coords, None, None, None, None, shell_lvl, None, None, None, None)

            d = 1e10  # infinitely large number
            closest = None

            # return vectors for further calculations (line_start[0], line_end[1], fire[2], es[3], fs[4], fe[5], se[6])
            # find references
            def vectors(line):
                l_start = np.array(line.start)
                l_end = np.array(line.end)
                fire = np.array(f_coords)

                return l_start, l_end, fire, l_end - l_start, fire - l_start, fire - l_end, l_start - l_end

            # iterate over lines to select the closest to the fire
            for l in lines:
                v = vectors(l)

                # orthogonal projection module
                def cos_vec(v1, v2):
                    return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

                if cos_vec(v[3], v[4]) >= 0 and cos_vec(v[6], v[5]) >= 0:  # choose projected point
                    d_iter = np.linalg.norm(np.cross(v[3], v[4])) / np.linalg.norm(v[3])
                else:  # choose the nearest edge if not
                    d_iter = min([np.linalg.norm(v[2] - v[0]), (np.linalg.norm(v[2] - v[0]))])

                # overwrite with analysed line if it is closer to the fire then the already chosen
                if d_iter < d:
                    d = d_iter
                    closest = l
            v = vectors(closest)  # generate vectors for selected line
            section = v[0] + (np.dot(v[4], v[3]) / np.dot(v[3], v[3])) * v[3]  # find the most exposed section coords

            unit_v = v[3] / np.linalg.norm(v[3])  # unit vector of selected line

            # set column's section to the biggest heat flux height (1.2m from the fire base)
            if element == 'c':
                # check if addition 1.2 m to the section Z is possible
                if section[-1] + 1.2 < max([v[1][-1], v[0][-1]]):
                    section += [0, 0, 1.2]
                else:
                    section[-1] = max([v[1][-1], v[0][-1]])
            generated = (*f_coords, *section, d, shell_lvl, closest.layer, *unit_v)

            self.geometry['f'].clear()  # clear temporary layout 'foo'

            # fire coords(list), section coords(list), length of the fire-section vector(float),
            # level of shell above the fire(float), profile(string), unit vector
            return generated

        # remove elements beneath the fire base or above shell level from the lines
        # cut those between the values
        def cut_lines(lines: list):
            for l in lines:
                # start point cannot be higher than end point
                if l.start[2] > l.end[2]:
                    start_rev = l.end
                    end_rev = l.start
                    l.start = start_rev
                    l.end = end_rev

                z1 = l.start[2]
                z2 = l.end[2]
                # do not consider lines beneath the fire base or above the ceiling
                if z2 < f_coords[2] or z1 >= shell_lvl:
                    continue
                # accept lines in (fire base, ceiling) ranges
                elif z1 >= f_coords[2] and z2 < shell_lvl:
                    self.geometry['f'].append(l)
                # cut lines to (fire base, ceiling) ranges
                else:
                    if z1 < f_coords[2]:
                        self.geometry['f'].append(l)
                        self.geometry['f'][-1].start = (l.start[0], l.start[1], f_coords[2])
                    if z2 >= shell_lvl:
                        self.geometry['f'].append(l)
                        self.geometry['f'][-1].end = (l.end[0], l.end[1], shell_lvl - 0.001)

        # checking if point consists in polygon (XY plane only)
        def ray_tracing_method(point: iter, poly: iter) -> bool:
            n = len(poly)
            inside = False
            x = point[0]
            y = point[1]

            p1x = poly[0][0]
            p1y = poly[0][1]
            for i in range(n + 1):
                p2x = poly[i % n][0]
                p2y = poly[i % n][1]
                if y > min(p1y, p2y):
                    if y <= max(p1y, p2y):
                        if x <= max(p1x, p2x):
                            if p1y != p2y:
                                xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                            if p1x == p2x or x <= xints:
                                inside = not inside
                p1x, p1y = p2x, p2y

            return inside

        shell_lvl = 1e5  # atmosphere bounds (here the space begins!)

        # check for shell (plate, ceiling) existing above the fire assign the level if true
        for s in self.geometry['s']:
            lvl = s.points[0][2]  # read level from first point of shell
            if float(f_coords[2]) <= lvl < shell_lvl and ray_tracing_method(f_coords, s.points):
                shell_lvl = lvl

        # initialize mapping
        if element == 'b':  # cut beams accordingly to Z in (fire_z - shell_lvl) range and map to relative
            cut_lines(self.geometry['b'])
        elif element == 'c':  # cut columns accordingly to Z in (fire_z - shell_lvl) range and map to relative
            cut_lines(self.geometry['c'])
        else:
            raise ValueError('[ERROR] {} is not a valid element type (\'b\' or \'c\' required)'.format(element))

        mapped = list(map_lines())

        self.prof_type = mapped[3]

        # (*fire coords, *section coords, length of the fire-section vector, level of shell above the fire, profile,
        # unit vector))
        return mapped


class Generator:
    def __init__(self, t_end, title, fire_type, fuelconfig):
        self.t_end = t_end  # simulation duration time
        self.title = title  # simulation title
        self.f_type = fire_type  # type of fire
        self.fire_coords = []  # to export to Single class

        print('Reading fuel configuration files...', end='\r')
        t = sec()
        if fuelconfig == 'stp':
            self.fuel = fires.Fuel(title).read_fuel()  # import fuel from STEP and FUL config files
        elif fuelconfig == 'obj':
            self.fuel = fires.FuelOBJ(title).read_fuel()  # import fuel from OBJ and FUL config files
        else:
            self.fuel = rcsv('{}.ful'.format(title))

        print('[OK] Fuel configuration imported ({} s)'.format(round(sec()-t, 2)))

    # import fire config
    def fire(self):
        # fire area is limited only by model limitation implemented to the fires.Properties
        mc_params, self.fire_coords = fires.f_localization(self.fuel)
        try:
            type, occupation = self.f_type.split('_')
        except ValueError:
            type = self.f_type
            occupation = None

        if 'alfat2' in type:
            f = fires.Fire(self.t_end, mc_params, self.fire_coords[2], occupation)
        elif 'sprink-eff' in type:
            f = fires.SprinkEff(self.t_end, mc_params, self.fire_coords[2], occupation)
        elif 'sprink-noeff' in type:
            f = fires.SprinkNoEff(self.t_end, mc_params, self.fire_coords[2], occupation)
        else:
            raise KeyError('[ERROR] {} is not a proper fire type'.format(self.f_type))

        return f.burn()

    def locafitxt(self, position, hrr_tab, diam_tab, z_ceil):
        # add data to LOCAFI.txt core
        lcf = core.locafi.copy()

        def ins(arg, index): return lcf[index].split('*')[0] + str(arg) + lcf[index].split('*')[1]

        str_position = '{} {} {}'.format(*position)

        lcf[2] = ins(str_position, 2)
        lcf[3] = ins(z_ceil, 3)

        # add tables of hrr and diameter
        def v2str(vector): return [str(i) for i in vector]

        def add(start, tab): [lcf.insert(i + start, '    '.join(v2str(tab[i])) + '\n') for i in range(len(tab))]

        add(lcf.index('DIAMETER\n') + 2, diam_tab)
        add(lcf.index('RHR\n') + 2, hrr_tab)

        # save locafi.txt to the dir
        with open('locafi.txt', 'w+') as file:
            file.writelines(lcf)


class MultiT2D:
    def __init__(self, time_end):
        self.t_end = time_end

    def dummy(self, chid, section, unit_v):

        # calculate nodes position
        np_section = np.array(section).astype(float)
        node1 = np_section - (unit_v / 1000)
        node2 = np_section + (unit_v / 1000)
        center = np_section
        # add perpendicular vector to section point
        if unit_v[0] != 0:
            lax = np_section + np.array([-unit_v[2] / unit_v[0], 0, 1])
        elif unit_v[1] != 0:
            lax = np_section + np.array([0, -unit_v[2] / unit_v[1], 1])
        elif unit_v[2] != 0:
            lax = np_section + np.array([1, 0, -unit_v[0] / unit_v[2]])
        else:
            raise ValueError('[ERROR] Zero length unit vector')

        # save nodes to a dummy.IN file
        lines = core.dummy.copy()

        def v2str(vector):
            return [str(i) for i in vector]

        def ins(index, arg):
            return lines[index].split('*')[0] + ' '.join(v2str(arg)) + lines[index].split('*')[1]

        for n in [(18, node1), (19, node2), (20, center), (21, lax)]:
            lines[n[0]] = ins(*n)

        # change T_END
        for n in (36, 41):
            lines[n] = str(self.t_end).join(lines[n].split('&T_END&'))

        with open('{}.in'.format(chid), 'w+') as file:
            file.writelines(lines)

    def profile(self, chid, profile_type):
        try:
            copyfile('{1}/{0}.gid/{0}.in'.format(profile_type, config['config_path']), '{}.in'.format(profile_type))
        except FileNotFoundError:
            raise FileNotFoundError('There is no {}.gid directory among configuration'.format(profile_type))

        # change profile to LCF and set chid.in as S3D
        Thermal(chid, 'LCF', frame_chid=chid, profile_pth='{}.in'.format(profile_type), time_end=self.t_end).change_in()

    # generate initial files (elem.in, prof.in, locafi.txt) based on DataFrame row (title_set.csv)
    def prepare(self, data_row):
        # create dummy.in (for one element S3D - just to map fire to element in SAFIR)
        unit_vector = np.array([data_row['u_x'], data_row['u_y'], data_row['u_z']]).astype(float)
        self.dummy(str(data_row['ID']), (data_row['x_s'], data_row['y_s'], data_row['z_s']), unit_vector)

        # copy profile GiD directory and prepare T2D files
        self.profile(str(data_row['ID']), str(data_row['profile']))


# generates a set of n scenarios
def generate_set(n, title, t_end, fire_type, config_path, results_path, fuelconfig):
    def create_df():
        return df(columns=('ID', 'element_type', 'time', 'x_f', 'y_f', 'z_f', 'x_s', 'y_s', 'z_s', 'distance',
                           'ceiling_lvl', 'profile', 'u_x', 'u_y', 'u_z', 'HRRPUA', 'alpha'))

    # append DataFrame to CSV file
    def df2csv(df, path='{}\{}_set.csv'.format(results_path, title)):
        try:
            with open(path):
                header = False
        except FileNotFoundError:
            header = True

        df.to_csv(path, mode='a', header=header)

    # create locafi.txt file
    def locafi(row, fire):
        chdir(config_path)

        # create simulation directory
        try:
            makedirs('{}\{}'.format(results_path, str(row['ID'])))
        except FileExistsError:
            pass
        chdir('{}\{}'.format(results_path, str(row['ID'])))

        # create locafi.txt fire file
        gen.locafitxt((row['x_f'], row['y_f'], row['z_f']), *fire, row['ceiling_lvl'])

        chdir(config_path)

    csvset = create_df()
    df2csv(csvset)
    simid_core = int(sec())

    if simid_core % 2 != 0:  # check if odd
        simid_core += 1

    print('[OK] User configuration imported')

    sing = Single(title)
    gen = Generator(t_end, title, fire_type, fuelconfig)

    t = sec()
    # draw MC input samples
    for i in range(0, int(n) * 2, 2):
        progressBar('Preparing fire scenarios', i, n * 2)
        fire = list(gen.fire())  # draw fire

        # draw localization of the most exposed beam
        csvset.loc[i] = [simid_core + i, 'b', ctime(sec())] + sing.generate(gen.fire_coords, 'b') + fire[2:]
        locafi(csvset.loc[i], fire[:2])  # generate locafi.txt

        # draw localization of the most exposed column
        csvset.loc[i + 1] = [simid_core + i + 1, 'c', ctime(sec())] + sing.generate(gen.fire_coords, 'c') + fire[2:]
        locafi(csvset.loc[i + 1], fire[:2])  # generate locafi.txt

        # write rows every 8 records (4 fire scenarios)
        if (i + 2) % 8 == 0:
            df2csv(csvset)
            del csvset
            csvset = create_df()

    # write unwritten rows
    try:
        df2csv(csvset)
        del csvset
    except ValueError:
        pass

    return '[OK] {} scenarios (2 simulations each) generated ({} s)'.format(int(n), round(sec() - t, 2))


# generate files for multisimulation
def generate_sim(data_path):
    t = sec()
    chdir(config['results_path'])
    data_set = rcsv(data_path)
    for i, r in data_set.iterrows():
        if r['profile'] != r['profile']:
            with open('{0}\{0}.err'.format(r['ID']), 'w') as err:
                mess = '[WARNING] There are no elements above the fire base in scenario {}'.format(r['ID'])
                err.write('{}\nMax element temperature in te scenario is equal to the ambient temperature'.format(mess))
            print(mess)
            continue
        chdir(str(r['ID']))
        MultiT2D(config['time_end']).prepare(r)
        chdir('..')

    return '[OK] {} simulation files created ({} s)'.format(len(data_set.index), round(sec() - t, 2))


if __name__ == '__main__':
    sys.stdout = Logger('mc.log')

    print('Reading user configuration...', end='\r')
    config = user_config(sys.argv[1])  # import multisimulation config
    try:
        makedirs(config['results_path'])  # create results directory
    except FileExistsError:
        pass

    chdir(config['config_path'])  # change to config

    print(generate_set(config['max_iterations'], config['case_title'], config['time_end'], config['fire_type'],
                       config['config_path'], config['results_path'], config['fuel']))
    print(generate_sim('{}\{}_set.csv'.format(config['results_path'], config['case_title'])))
