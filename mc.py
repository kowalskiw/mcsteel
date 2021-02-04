from time import time as current_seconds
import numpy as np
from fires import f_localization, Properties
import ezdxf
import shapely.geometry as sh
from pandas import DataFrame as df
from pandas import read_csv as rcsv
import core
from os import mkdir, chdir, getcwd
from shutil import copyfile
from fdsafir import Thermal


'''Read geometry and map it to the fire (choose the most exposed sections)'''


class Single:
    def __init__(self, title, element_type='b'):
        self.title = title     # simulation title
        self.f_type = 'alfat2_store'    # type of fire
        self.t_end = 1800      # simulation duration time
        self.element_type = element_type  # 'b' - beam, 'c' - column
        self.prof_type = 'invalid profile type'     # type of steel profile
        self.fire_properties = f_localization(title)[1]

    # read dxf geometry
    def read_dxf(self):
        def dxf_to_shapely(dxfshells):    # convert shells to shapely objects
            shshells = []

            for s in dxfshells:
                shshell = sh.Polygon([s[0], s[1], s[2], s[3]])
                shshells.append(shshell)

            return shshells

        dxffile = ezdxf.readfile('{}.dxf'.format(self.title))
        msp = dxffile.modelspace()
        columns = []
        beams = []
        for l in msp.query('LINE'):     # assign LINES elements to columns or beams tables
            if l.dxf.start[:-1] == l.dxf.end[:-1]:
                columns.append(l)
            else:
                beams.append(l)

        # assign 3DFACE elements to shells table
        shells = []
        [shells.append(s) for s in msp.query('3DFACE')]

        return beams, columns, dxf_to_shapely(shells)

    # map fire and geometry - choose fire scenario
    def map(self, f_coords, elements, element):
        # select the most exposed section among the lines and return its properties
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

            # iterate over lines to select the closest to the fire one
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

            # fire_relative = section - fire  # calculate fire-section vector

            # fire coords(list), section coords(list), length of the fire-section vector(float),
            # level of shell above the fire(float), profile(string), unit vector
            return (*tuple(*fire.coords), *section, d, shell_lvl, lines[index].dxf.layer, *unit_v)

        # remove elements beneath the fire base or above shell level from the lines
        def cut_lines(lines):
            for l in lines:
                if l.dxf.start[2] > shell_lvl or l.dxf.end[2] < f_coords[2]:
                    lines.remove(l)
            return lines

        # import fire and structure data
        beams, columns, shells = elements  # shells -> dict{Z_level:Plygon} | lines
        fire = sh.Point(list(f_coords))  # shapely does not support the 3D objects - z coordinate is not used
        shell_lvl = 1e9     # infinitely large number

        # check for shell (plate, ceiling) existing above the fire assign the level if true
        for poly in shells:
            lvl = float(poly.wkt.split()[-1][:-2])
            if float(f_coords[2]) <= lvl < shell_lvl and poly.contains(fire):
                shell_lvl = lvl
        if shell_lvl == 1e6:  # set shell_lvl as -1 when there is no plate above the fire (besides compartment ceiling)
            shell_lvl = -1

        # initialize mapping
        if element == 'b':  # cut beams accordingly to Z in (fire_z - shell_lvl) range and map to relative
            mapped = map_lines(cut_lines(beams))
        elif element == 'c':  # cut columns accordingly to Z in (fire_z - shell_lvl) range and map to relative
            mapped = map_lines(cut_lines(columns))
        else:
            raise ValueError('{} is not a proper element type'.format(element))

        self.prof_type = mapped[3]

        # (*fire coords, *section coords, length of the fire-section vector, level of shell above the fire, profile,
        # unit vector))
        return mapped

    def generate(self, element_type='b'):
        return list(self.map(self.fire_properties, self.read_dxf(), element=element_type))


class Generator:
    def __init__(self):
        self.t_end = 1800      # simulation duration time
        self.title = 'foo'     # simulation title
        self.f_type = 'alfat2_store'    # type of fire
        pass

# import fire properties
    def fire(self):

        # fire area is limited only by model limitation implemented to the fires.Properties
        f = Properties(self.t_end)
        if self.f_type == 'alfat2':
            f_hrr, f_diam, hrrpua, alpha = f.alfa_t2(self.title)
        elif self.f_type == 'alfat2_store':
            f_hrr, f_diam, hrrpua, alpha = f.alfa_t2(self.title, property='store')
        elif self.f_type == 'sprink-eff':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_eff(self.title)
        elif self.f_type == 'sprink-eff_store':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_eff(self.title, property='store')
        elif self.f_type == 'sprink-noeff':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_noeff(self.title)
        elif self.f_type == 'sprink-noeff_store':
            f_hrr, f_diam, hrrpua, alpha = f.sprink_noeff(self.title, property='store')
        else:
            raise KeyError('{} is not a proper fire type'.format(self.f_type))

        print('alpha={}\nHRRPUA={}'.format(alpha, hrrpua))

        return f_hrr, f_diam

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
        # directory == 'simulation/multit2d/chid'
        self.title = 'foo'

    def dummy(self, chid, section, unit_v):

        # calculate nodes position
        np_section = np.array(section)
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
        copyfile('../{0}.gid/{0}.in'.format(profile_type), '{}.in'.format(profile_type))

        # change profile to LCF and set chid.in as S3D
        Thermal(chid, 'LCF', frame_chid=chid, profile_pth='{}.in'.format(profile_type)).change_in()

    # generate initial files (elem.in, prof.in, locafi.txt) based on DataFrame row (title_set.csv)
    def prepare(self, data_row):
        # create dummy.in (for one element S3D - just to map fire to element in SAFIR)
        unit_vector = np.array([data_row['u_x'], data_row['u_y'], data_row['u_z']])
        # try:
        self.dummy(str(data_row['ID']), (data_row['x_s'], data_row['y_s'], data_row['z_s']), unit_vector)
        # except IndexError:
        #     print('list error')
        #     pass

        # copy profile GiD directory and prepare T2D files
        self.profile(str(data_row['ID']), str(data_row['profile']))


# generates a set of n scenarios
def generate_set(n, title):
    csvset = df(columns=('ID', 'element_type', 'x_f', 'y_f', 'z_f', 'x_s', 'y_s', 'z_s', 'distance', 'ceiling_lvl',
                         'profile', 'u_x', 'u_y', 'u_z'))
    simid_core = int(current_seconds())

    for i in range(0, n*2, 2):  # add MC-drawn records to the DataFrame
        sing = Single(title)
        csvset.loc[i] = [simid_core+i, 'b'] + sing.generate()
        csvset.loc[i+1] = [simid_core+i+1, 'c'] + sing.generate(element_type='c')

    csvset.to_csv('{}_set.csv'.format(title))    # save to the csv file

    # generate locafi.txt file
    for i, row in csvset.iterrows():
        gen = Generator()
        fire = gen.fire()
        chid = str(row['ID'])

        # create simulation directory
        try:
            mkdir(chid)
        except FileExistsError:
            print('directory already exists')
        chdir(chid)

        # create locafi.txt fire file
        gen.locafitxt((row['x_f'], row['y_f'], row['z_f']), *fire, row['ceiling_lvl'])

        chdir('..')

    return 0


generate_set(10, 'foo')

data_set = rcsv('foo_set.csv')
for i, r in data_set.iterrows():
    chdir(str(r['ID']))
    print(getcwd())
    MultiT2D().prepare(r)
    chdir('..')

