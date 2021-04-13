from time import time as current_seconds
from time import ctime
import numpy as np
import ezdxf
from pandas import DataFrame as df
from pandas import read_csv as rcsv
import core
from os import mkdir, chdir
from shutil import copyfile
from sys import argv
from fdsafir import Thermal, user_config
from fires import f_localization, Properties, Fuel

'''Read geometry and map it to the fire (choose the most exposed sections)'''


class Single:
    def __init__(self, title, fire_type, fire_coords):
        self.title = title     # simulation title
        self.f_type = fire_type    # type of fire
        self.prof_type = 'invalid profile type'     # type of steel profile
        self.fire_coords = fire_coords        # dummy fire coordinates
        self.geometry = self.read_dxf()     # import geometry from DXF file

    # read dxf geometry
    def read_dxf(self):
        def dxf_to_list(dxfshell: ezdxf.entities.solid.Face3d) -> list:
            nodes = []
            i = 0
            try:
                while True:
                    nodes.append(tuple([dxfshell[i][j] for j in range(3)]))
                    i += 1
            except IndexError:
                return nodes

        dxffile = ezdxf.readfile('{}.dxf'.format(self.title))
        print('Reading DXF geometry...')
        msp = dxffile.modelspace()
        columns = []
        beams = []
        for l in msp.query('LINE'):     # assign LINES elements to columns or beams tables
            if l.dxf.start[0] == l.dxf.end[0] and l.dxf.start[1] == l.dxf.end[1]:
                columns.append(l)
            else:
                beams.append(l)

        # assign 3DFACE elements to shells table
        shells = []
        for s in msp.query('3DFACE'):
            pass
        [shells.append(dxf_to_list(s)) for s in msp.query('3DFACE')]

        print(shells)

        print('[OK] DXF geometry imported')

        return beams, columns, shells

    # map fire and geometry - choose fire scenario
    def map(self, f_coords, elements, element):
        # select the most exposed section among the lines and return its config
        def map_lines(lines):
            d = 1e9  # infinitely large number
            index = None

            # return vectors for further calculations (line_start[0], line_end[1], fire[2], es[3], fs[4], fe[5], se[6])
            # find references
            def vectors(line):
                l_start = np.array(line.dxf.start)
                l_end = np.array(line.dxf.end)
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
                    index = lines.index(l)

            v = vectors(lines[index])  # generate vectors for selected line
            section = v[0] + (np.dot(v[4], v[3]) / np.dot(v[3], v[3])) * v[3]  # find the most exposed section coords

            unit_v = v[3] / np.linalg.norm(v[3])    # unit vector of selected line

            # set column's section to the biggest heat flux height (1.2m from the fire base)
            if element == 'c':
                # check if addition 1.2 m to the section Z is possible
                if section[-1] + 1.2 < max([v[1][-1], v[0][-1]]):
                    section += [0, 0, 1.2]
                else:
                    section[-1] = max([v[1][-1], v[0][-1]])

            # fire coords(list), section coords(list), length of the fire-section vector(float),
            # level of shell above the fire(float), profile(string), unit vector
            return (*f_coords, *section, d, shell_lvl, lines[index].dxf.layer, *unit_v)

        # remove elements beneath the fire base or above shell level from the lines
        def cut_lines(lines):
            for l in lines:
                if l.dxf.start[2] > shell_lvl or l.dxf.end[2] < f_coords[2]:
                    lines.remove(l)
            return lines

        # checking if point consists in polygon (XY plane only)
        def ray_tracing_method(point: iter, poly: list) -> bool:
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

        # import fire and structure data
        beams, columns, shells = elements  # shells -> dict{Z_level:Plygon} | lines
        shell_lvl = 1e9     # infinitely large number

        # check for shell (plate, ceiling) existing above the fire assign the level if true
        for s in shells:
            lvl = s[0][2]   # read level from first point of shell
            if float(f_coords[2]) <= lvl < shell_lvl and ray_tracing_method(f_coords, s):
                shell_lvl = lvl
        if shell_lvl == 1e6:  # set shell_lvl as -1 when there is no plate above the fire (besides compartment ceiling)
            shell_lvl = -1

        # initialize mapping
        if element == 'b':  # cut beams accordingly to Z in (fire_z - shell_lvl) range and map to relative
            mapped = map_lines(cut_lines(beams))
        elif element == 'c':  # cut columns accordingly to Z in (fire_z - shell_lvl) range and map to relative
            mapped = map_lines(cut_lines(columns))
        else:
            raise ValueError('{} is not a proper element type (\'b\' or \'c\' required)'.format(element))

        self.prof_type = mapped[3]

        # (*fire coords, *section coords, length of the fire-section vector, level of shell above the fire, profile,
        # unit vector))
        return mapped

    def generate(self, element_type='b'):
        return list(self.map(self.fire_coords, self.geometry, element=element_type))


class Generator:
    def __init__(self, t_end, title, fire_type, fuelconfig='fuel&stp'):
        self.t_end = t_end      # simulation duration time
        self.title = title     # simulation title
        self.f_type = fire_type    # type of fire

        print('Reading fuel configuration files...')

        if fuelconfig == 'fuel&stp':
            self.fuel = Fuel(title).read_fuel()  # import fuel from config files
        else:
            self.fuel = rcsv('{}.ful'.format(title))
        self.fire_props, self.fire_coords = f_localization(self.fuel)

        print('[OK] Fuel configuration files imported')

# import fire config
    def fire(self):
        # fire area is limited only by model limitation implemented to the fires.Properties
        f = Properties(self.t_end, self.fire_props)
        if self.f_type == 'alfat2':
            f_hrr, f_diam, hrrpua, alpha = f.alfa_t2()
        elif self.f_type == 'alfat2_store':
            f_hrr, f_diam, hrrpua, alpha = f.alfa_t2(property='store')
        elif self.f_type == 'sprink-eff':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_eff()
        elif self.f_type == 'sprink-eff_store':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_eff(property='store')
        elif self.f_type == 'sprink-noeff':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_noeff()
        elif self.f_type == 'sprink-noeff_store':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_noeff(property='store')
        else:
            raise KeyError('{} is not a proper fire type'.format(self.f_type))

        print('alpha={}\nHRRPUA={}'.format(alpha, hrrpua))

        return f_hrr, f_diam, hrrpua, alpha

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
    def __init__(self):
        pass

    def dummy(self, chid, section, unit_v):

        # calculate nodes position
        np_section = np.array(section).astype(float)
        node1 = np_section - (unit_v / 1000)
        node2 = np_section + (unit_v / 1000)
        center = np_section
        # add perpendicular vector to section point
        if unit_v[0] != 0:
            lax = np_section + np.array([-unit_v[2]/unit_v[0], 0, 1])
        elif unit_v[1] != 0:
            lax = np_section + np.array([0, -unit_v[2]/unit_v[1], 1])
        elif unit_v[2] != 0:
            lax = np_section + np.array([1, 0, -unit_v[0] / unit_v[2]])
        else:
            raise ValueError('Zero length unit vector')

        # save nodes to a dummy.IN file
        lines = core.dummy.copy()
        # print(core.dummy)
        # print(lines)
        def v2str(vector): return [str(i) for i in vector]
        def ins(index, arg): return lines[index].split('*')[0] + ' '.join(v2str(arg)) + lines[index].split('*')[1]
        for n in [(18, node1), (19, node2), (20, center), (21, lax)]:
            # print(lines[n[0]])
            lines[n[0]] = ins(*n)

        with open('{}.in'.format(chid), 'w+') as file:
            file.writelines(lines)

    def profile(self, chid, profile_type):
        try:
            copyfile('{1}/{0}.gid/{0}.in'.format(profile_type, config['config_path']), '{}.in'.format(profile_type))
        except FileNotFoundError:
            raise FileNotFoundError('There is no {}.gid directory among configuration'.format(profile_type))

        # change profile to LCF and set chid.in as S3D
        Thermal(chid, 'LCF', frame_chid=chid, profile_pth='{}.in'.format(profile_type)).change_in()

    # generate initial files (elem.in, prof.in, locafi.txt) based on DataFrame row (title_set.csv)
    def prepare(self, data_row):
        # create dummy.in (for one element S3D - just to map fire to element in SAFIR)
        unit_vector = np.array([data_row['u_x'], data_row['u_y'], data_row['u_z']]).astype(float)
        self.dummy(str(data_row['ID']), (data_row['x_s'], data_row['y_s'], data_row['z_s']), unit_vector)

        # copy profile GiD directory and prepare T2D files
        self.profile(str(data_row['ID']), str(data_row['profile']))


# generates a set of n scenarios
def generate_set(n, title, t_end, fire_type, config_path, results_path):
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
        print('{} records written to {}_set.csv'.format(len(csvset), title))

    # create locafi.txt file
    def locafi(row, fire):
        chdir(config_path)

        # create simulation directory
        try:
            mkdir('{}\{}'.format(results_path, str(row['ID'])))
        except FileExistsError:
            print('directory already exists')
        chdir('{}\{}'.format(results_path, str(row['ID'])))

        # create locafi.txt fire file
        gen.locafitxt((row['x_f'], row['y_f'], row['z_f']), *fire, row['ceiling_lvl'])

        chdir(config_path)

    csvset = create_df()
    df2csv(csvset)
    simid_core = int(current_seconds())
    gen = Generator(t_end, title, fire_type)
    sing = Single(title, fire_type, gen.fire_coords)

    # draw MC input samples
    for i in range(0, int(n)*2, 2):
        fire = list(gen.fire())   # draw fire

        # draw localization of the most exposed beam
        csvset.loc[i] = [simid_core+i, 'b', ctime(current_seconds())] + sing.generate() + fire[2:]
        locafi(csvset.loc[i], fire[:2])     # generate locafi.txt

        # draw localization of the most exposed column
        csvset.loc[i+1] = [simid_core+i+1, 'c', ctime(current_seconds())] + sing.generate(element_type='c') + fire[2:]
        locafi(csvset.loc[i + 1], fire[:2])    # generate locafi.txt

        # write rows every 8 records (4 fire scenarios)
        if (i+2) % 8 == 0:
            df2csv(csvset)
            del csvset
            csvset = create_df()

    # write unwritten rows
    try:
        df2csv(csvset)
        del csvset
    except ValueError:
        print('No more data to be written')
        pass

    return 0


# generate files for multisimulation
def generate_sim(data_path):
    chdir(config['results_path'])
    data_set = rcsv(data_path)
    for i, r in data_set.iterrows():
        chdir(str(r['ID']))
        MultiT2D().prepare(r)
        chdir('..')

    return 0


if __name__ == '__main__':
    config = user_config(argv[1])   # import multisimulation config
    try:
        mkdir(config['results_path'])   # create results directory
    except FileExistsError:
        pass

    chdir(config['config_path'])    # change to config

    generate_set(config['max_iterations'], config['case_title'], config['time_end'], config['fire_type'],
                 config['config_path'], config['results_path'])
    generate_sim('{}\{}_set.csv'.format(config['results_path'], config['case_title']))
