from time import time as current_seconds
import numpy as np
from fires import Properties
import ezdxf
import shapely.geometry as sh
from pandas import DataFrame as df


'''Read geometry and map it to the fire (choose the most exposed sections)'''


class Single:
    def __init__(self, title, element_type='b'):
        self.title = title     # simulation title
        self.f_type = ''    # type of fire
        self.t_end = 0      # simulation duration time
        self.element_type = element_type  # 'b' - beam, 'c' - column

    # read dxf geometry
    def read_dxf(self):
        def dxf_to_shapely(dxfshells):    # convert shells to shapely objects
            shshells = {}

            for s in dxfshells:
                shshell = sh.Polygon([s[0], s[1], s[2], s[3]])
                shshells[str(s[0][-1])] = shshell

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

    # import fire properties
    def f_properties(self):

        # fire area is limited only by model limitation implemented to the fires.Properties

        f = Properties(self.t_end)
        if self.f_type == 'alfat2':
            f_hrr, f_diam, f_coordinates, hrrpua, alpha = f.alfa_t2(self.title)
        elif self.f_type == 'alfat2_store':
            f_hrr, f_diam, f_coordinates, hrrpua, alpha = f.alfa_t2(self.title, property='store')
        elif self.f_type == 'sprink-eff':
            f_hrr, f_diam, f_coordinates, hrrpua, alpha = f.sprink_eff(self.title)
        elif self.f_type == 'sprink-eff_store':
            f_hrr, f_diam, f_coordinates, hrrpua, alpha = f.sprink_eff(self.title, property='store')
        elif self.f_type == 'sprink-noeff':
            f_hrr, f_diam, f_coordinates, hrrpua, alpha = f.sprink_noeff(self.title)
        elif self.f_type == 'sprink-noeff_store':
            f_hrr, f_diam, f_coordinates, hrrpua, alpha = f.sprink_noeff(self.title, property='store')
        else:
            raise(KeyError, '{} is not a proper fire type'.format(self.f_type))

        return f_coordinates, f_hrr, f_diam

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

            # set column's section to the biggest heat flux height (1.2m from the fire base)
            if element == 'c':
                # check if addition 1.2 m to the section Z is possible
                if section[-1] + 1.2 < max([v[1][-1], v[0][-1]]):
                    section += [0, 0, 1.2]
                else:
                    section[-1] = max([v[1][-1], v[0][-1]])

            fire_relative = section - fire  # calculate fire-section vector

            # fire coords, section coords, length of the fire-section vector, level of shell above the fire, profile
            return fire, section, d, shell_lvl, lines[index].dxf.layer

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
        for lvl, poly in shells.items():
            lvl = float(lvl)
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

        # # returns tuple (x_r, y_r, z_r, [distance3D, LOCAFI_h, 'h', profile, shell height])
        # return (*mapped[0], [mapped[1], mapped[2] - fire_z, element, mapped[3], mapped[2]])

        # (fire coords, section coords, length of the fire-section vector, level of shell above the fire, profile)
        return mapped

    def generate(self): return self.map(self.f_properties()[0], self.read_dxf(), element=self.element_type)


# generates set of n scenarios
def generate_set(n, title):
    csvset = df(columns=('ID', 'x_f', 'y_f', 'z_f', 'x_s', 'y_s', 'z_s', 'distance', 'ceiling_lvl', 'profile'))
    simid_core = int(current_seconds())

    for i in range(n*2, step=2):  # add MC-drawn records to the DataFrame
        csvset.loc[i] = [simid_core+i] + list(Single(title).generate())
        csvset.loc[i+1] = [simid_core+i] + list(Single(title, element_type='c').generate())

    csvset.to_csv('{}_set.csv'.format(title))    # save to the csv file


# generates initial files (elem.in, prof.in, locafi.txt) based on DataFrame row (title_set.csv)
def generate_in(sim_data):
    # create simulation directory
    # copy profile GiD directory
    # create dummy.in (for one element S3D - just to map fire to element in SAFIR)
    # generate locafi.txt file
    pass

