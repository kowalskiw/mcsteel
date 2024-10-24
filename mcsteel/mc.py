import os.path
from os import makedirs, scandir
from time import time as sec
from time import ctime
import numpy as np
import dxfgrabber as dxf
from pandas import DataFrame as df
from pandas import read_csv
from shutil import copy2
import sys
from numpy import random
from statistics import fmean

from mcsteel.utils import ThermalTEM, Config, progress_bar, out, triangular
import mcsteel.core
import mcsteel.fires

global outpth


class FireScenario:
    def __init__(self, config_object: Config, fire_properties, sprinkler_activation):
        self.config = config_object  # config class
        self.fire_location = fire_properties[2]  # [x, y, z]
        self.alpha = fire_properties[0]  # [W/s^2]
        self.hrrpua = fire_properties[1]  # [W/m^2]
        self.sprinklers = sprinkler_activation  # [s]
        self.profiles = []  # [profile1, profile2, profile3] profile1 = fdsafir2.ThermalTEM
        self.fire_curve = [[], []]  # [[0,... time steps ... t_end], [HRR(0), ... HRR ... HRR(t_end)]]
        self.fire_type = config_object.fire_type  # fire curve function type form fires.Fires
        self.mapped = []  # complete set of data for profiles to be calculated in this scenario
        self.ceiling = 1e5  # level of ceiling above the fire source (here the space begins!)
        self.locafi_lines = []  # lines for locafi.txt fire file

    # map fire location with structure to find the most heated profiles to be analysed
    def map(self, structure):
        # (*fire coords, *section coords, length of the fire-section vector, level of shell above the fire, profile,
        # *unit vector))

        # select the most exposed section among the lines and return its config
        def map_lines(element):
            # no element exception
            lins = structure['f']
            if lins.__len__() == 0:
                return [*self.fire_location, None, None, None, None, self.ceiling, None, None, None, None, None, None]

            d = 1e10  # infinitely large number
            closest = None

            # return vectors for further calculations (line_start[0], line_end[1], fire[2], es[3], fs[4], fe[5], se[6])
            # find references
            def vectors(single_line):
                l_start = np.array(single_line.start)
                l_end = np.array(single_line.end)
                fire = np.array(self.fire_location)

                return l_start, l_end, fire, l_end - l_start, fire - l_start, fire - l_end, l_start - l_end

            # iterate over lines to select the closest to the fire
            for line in lins:
                v = vectors(line)

                # here begins orthogonal projection module
                # calculate cosine between two vectors
                def cos_vec(v1, v2):
                    return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

                # check if point of fire source has an orthogonal projection on the line
                if cos_vec(v[3], v[4]) >= 0 and cos_vec(v[6], v[5]) >= 0:
                    # calculate distance from point to its orthogonal projection on line
                    d_iter = np.linalg.norm(np.cross(v[3], v[4])) / np.linalg.norm(v[3])
                else:  # choose the nearest edge if not
                    d_iter = min([np.linalg.norm(v[2] - v[0]), (np.linalg.norm(v[2] - v[0]))])

                # overwrite with analysed line if it is closer to the fire then the already chosen
                if d_iter < d:
                    d = d_iter
                    closest = line
            v = vectors(closest)  # generate vectors for selected line
            section = v[0] + min([1, (np.dot(v[4], v[3]) / np.dot(v[3], v[3]))]) * v[3]

            unit_v = v[3] / np.linalg.norm(v[3])  # unit vector of selected line

            # set column's section to the biggest heat flux height (1.2m from the fire base)
            if element == 'c':
                # check if addition 1.2 m to the section Z is possible
                if section[-1] + 1.2 < max([v[1][-1], v[0][-1]]):
                    section += [0, 0, 1.2]
                else:
                    section[-1] = max([v[1][-1], v[0][-1]])
            generated = [*self.fire_location, *section, d, self.ceiling, closest.layer.split('*')[0], *unit_v,
                         self.hrrpua, self.alpha]

            structure['f'].clear()  # clear temporary layout 'foo'

            # fire coords(list), section coords(list), length of the fire-section vector(float),
            # level of shell above the fire(float), profile(string), unit vector
            return generated

        # remove elements beneath the fire base or above shell level from the lines
        # cut those between the values
        def cut_lines():
            for line in lines:
                # start point cannot be higher than end point
                if line.start[2] > line.end[2]:
                    start_rev = line.end
                    end_rev = line.start
                    line.start = start_rev
                    line.end = end_rev

                z1 = line.start[2]
                z2 = line.end[2]
                # do not consider lines beneath the fire base or above the ceiling
                if z2 <= self.fire_location[2] or z1 >= self.ceiling:
                    continue
                # accept lines in (fire base, ceiling) ranges
                elif z1 > self.fire_location[2] and z2 < self.ceiling:
                    structure['f'].append(line)
                # cut lines to (fire base, ceiling) ranges with 0.01 tolerance
                else:
                    to_save = None
                    if z1 <= self.fire_location[2]:
                        to_save = line
                        to_save.start = (line.start[0], line.start[1], self.fire_location[2] + 0.01)
                    if z2 >= self.ceiling:
                        to_save = line
                        to_save.end = (line.end[0], line.end[1], self.ceiling - 0.01)
                    # check if line has non-zero length
                    if np.linalg.norm(np.array(to_save.start) - np.array(to_save.end)) > 0:
                        structure['f'].append(to_save)

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
                            xints = None
                            if p1y != p2y:
                                xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                            if p1x == p2x or x <= xints:
                                inside = not inside
                p1x, p1y = p2x, p2y

            return inside

        # check for shell (plate, ceiling) existing above the fire assign the level if true
        for s in structure['s']:
            lvl = s.points[0][2]  # read level from first point of shell
            if float(self.fire_location[2]) <= lvl < self.ceiling and ray_tracing_method(self.fire_location, s.points):
                self.ceiling = lvl

        for i, element_type in enumerate(['b', 'c']):
            lines = structure[element_type]  # choose beams or columns as lines
            cut_lines()  # cut beams accordingly to Z in (fire_z - shell_lvl) range and map to relative

            self.mapped.append(map_lines(element_type))

            self.profiles.append(self.mapped[i][3])

    # calculate HRR(t) and D(t) tables
    def create_fire_curve(self):
        # fire area is limited only by model limitation implemented to the fires.Properties
        if 'alfat2' in self.fire_type:
            f = mcsteel.fires.AlfaT2(self.config.time_end, self.alpha, self.hrrpua)
        elif 'sprink-eff' in type:
            f = mcsteel.fires.SprinkEff(self.config.time_end, self.alpha, self.hrrpua, self.sprinklers)
        elif 'sprink-noeff' in type:
            f = mcsteel.fires.SprinkNoEff(self.config.time_end, self.alpha, self.hrrpua, self.sprinklers)
        else:
            raise KeyError(f'[ERROR] {self.fire_type} is not a proper fire type')

        self.fire_curve = f.burn()


class Iteration:
    def __init__(self, scenario: FireScenario, iteration_data: list, chid: str):
        self.chid = chid
        self.fire = scenario
        self.config = scenario.config
        self.dir_path = os.path.join(self.config.results_path, chid)

        # [*self.fire_location, *section, d, self.ceiling, closest.layer.split('*')[0], *unit_v]
        self.data = iteration_data

    # fill locafi.txt template from core.py
    def prepare_locafi(self):
        if not all(self.fire.fire_curve):
            raise IndexError('[ERROR] Empty locafi_lines list')
        str_position = f'{"    ".join([str(c) for c in self.fire.fire_location])}'
        lcf = mcsteel.core.locafi.copy()

        # insert value to the template
        def ins(arg, index): return lcf[index].split('*')[0] + str(arg) + lcf[index].split('*')[1]

        # change vector to the list of string values
        def v2str(vector): return [str(i) for i in vector]

        # add tables of hrr and diameter to the template
        def add(start, tab): [lcf.insert(i + start, '    '.join(v2str(tab[i])) + '\n') for i in range(len(tab))]

        lcf[2] = ins(str_position, 2)
        lcf[3] = ins(str(self.fire.ceiling), 3)
        [add(lcf.index(t) + 2, self.fire.fire_curve[i]) for i, t in enumerate(['RHR\n', 'DIAMETER\n'])]

        with open(os.path.join(self.dir_path, 'locafi.txt'), 'w') as lcffile:
            lcffile.writelines(lcf)

    def copy_section(self, section_chid):
        copy2(self.config.section_path(section_chid), self.dir_path)
        ThermalTEM(1, [section_chid, [], []], self.config.config_path, 'lcf',
                   self.config.time_end, self.dir_path).change_in(os.path.basename(self.dir_path))

    def write_dummy_structural(self, section: list, unit_v: np.array):
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
        lines = mcsteel.core.dummy.copy()

        def v2str(vector):
            return [str(i) for i in vector]

        def ins(index, arg):
            return lines[index].split('*')[0] + ' '.join(v2str(arg)) + lines[index].split('*')[1]

        for n in [(18, node1), (19, node2), (20, center), (21, lax)]:
            lines[n[0]] = ins(*n)

        # change T_END
        for n in (36, 41):
            lines[n] = str(self.config.time_end).join(lines[n].split('&T_END&'))

        with open(os.path.join(self.dir_path, f'{self.chid}.in'), 'w+') as file:
            file.writelines(lines)

    # allows to prepare files in different location (i.e. node)
    def prepare_files(self):
        # check if there are any elements above the fire source
        if self.data[-4] != self.data[-4]:
            with open(os.path.join(self.dir_path, f'{self.chid}.err'), 'w') as err:
                mess = f'[WARNING] There are no structural elements above the fire base in the' \
                       f' {self.chid} fire scenario'
                err.write(f'{mess}\nMax element temperature in this scenario is equal to the ambient temperature')
            out(outpth, mess)

        makedirs(self.dir_path)

        # save fire file to the directory
        self.prepare_locafi()

        # create SAFIR files
        self.copy_section(self.data[-6])
        self.write_dummy_structural(self.data[3:6], np.array(self.data[-3:]).astype(float))


class MCGenerator:
    def __init__(self, config_object: Config):
        self.config = config_object
        self.n = self.config.max_iterations if self.config.max_iterations else 100  # size of the sample
        self.fuel = self._read_fuel() if self.config.fire_type.lower() != 'cfast' else None  # fuel distribution and properties
        self.set = []

    def _read_fuel(self):
        out(outpth, 'Reading fuel configuration files...\r')
        t0 = sec()

        fuel_type = self.config.fuel.lower()
        fuel_full_path = os.path.join(self.config.config_path, f'{self.config.title}.fuel')
        if fuel_type == 'obj':
            fuel = mcsteel.fires.FuelOBJ(fuel_full_path).read_fuel()
        elif fuel_type == 'step':
            fuel = mcsteel.fires.Fuel(fuel_full_path).read_fuel()
        else:
            fuel = mcsteel.fires.OldFuel(fuel_full_path).read_fuel()

        out(outpth, f'[OK] Fuel configuration imported ({round(sec() - t0, 2)} s)                      ')
        return fuel

    @staticmethod
    def _find_hrrpua(fire_z, properties):
        # calculate HRRPUA according to triangular distribution specified by user
        upward = properties.ZB - fire_z
        downward = fire_z - properties.ZA
        # scale HRRPUA acc. to NFPA204 (1/10 downwards)
        reduction = (upward + downward / 10) / properties.hrrpua_height
        return reduction * triangular(properties.hrrpua_min, properties.hrrpua_max, mode=properties.hrrpua_mode)

    # calculate ALPHA according to the experimental log-norm or user's triangular distribution
    def _find_alpha(self, hrrpua, properties):
        if 'store' in {self.config.occupancy, self.config.fire_type}:
            return hrrpua * random.lognormal(-9.72, 0.97)  # [kW/s2]
        else:
            return triangular(properties.alpha_min, properties.alpha_max, mode=properties.alpha_mode)  # [kW/s2]

    # find fire localization and properites of fuel in that place
    def _find_fire_origin(self):
        def random_position(xes, yes, zes):
            coordinates = []
            [coordinates.append(random.randint(int(10 * i[0]), int(10 * i[1])) / 10) for i in [xes, yes, zes]]
            return coordinates

        # find the fuel actual_site within the fuel sites
        def find_site(fuel_sites):
            ases = []  # list with partial factors A of each fuel area
            probs = []  # list with probabilities of ignition in each fuel area

            # calculate partial factor A (area * probability) of each fuel area
            for i, r in fuel_sites.iterrows():
                a = (r['XB'] - r['XA']) * (r['YB'] - r['YA']) * r['MC']
                ases.append(a)

            # calculate probability of ignition in each fuel area
            for a in ases:
                probs.append(a / sum(ases))

            # return drawn fuel area
            return random.choice(len(probs), p=probs)

        site_no = find_site(self.fuel)  # generate fire coordinates from MC function

        site = self.fuel.iloc[site_no]  # information about chosen fuel actual_site

        return site, random_position((site.XA, site.XB), (site.YA, site.YB), zes=(site.ZA, site.ZB))

    # gather all fire properties
    def _find_fire(self):
        fuel_properties, fire_coordinates = self._find_fire_origin()
        hrrpua = self._find_hrrpua(fire_coordinates[2], fuel_properties)
        alpha = self._find_alpha(hrrpua, fuel_properties)
        try:
            sprink_act = fuel_properties.t_sprink
        except KeyError:
            sprink_act = None

        return [alpha, hrrpua, fire_coordinates], sprink_act

    # main MC-sampling function for fire
    def sampling(self):
        t = sec()

        if self.config.fire_type.lower() == 'cfast':
            cfast_subdirs = []
            for i in scandir(os.path.join(self.config.config_path, 'cfast')):
                cfast_subdirs.append(i.path) if i.is_dir() else None

        for i in range(self.n):
            progress_bar('Monte Carlo sampling', i, self.n)
            if self.config.fire_type.lower() == 'cfast':
                try:
                    self.set.append(CFASTScenario(self.config, cfast_subdirs[i]))
                except IndexError:
                    out(outpth, f'[WARNING] Not enough CFAST files. {self.n} fire scenarios were requested in .USER file.'
                                f' Proceeding with {i} scenarios')
            else:
                self.set.append(FireScenario(self.config, *self._find_fire()))
        out(outpth, f'[OK] {self.n} fire scenarios were chosen ({round(sec() - t, 3)}) s                  ')
        return self.set


class CFASTScenario(FireScenario):
    def __init__(self, config_object: Config, path_to_cfast_dir, chid='cfast'):
        self.cfast_dir = path_to_cfast_dir
        self.cfast_chid = chid
        self.cfast_fire_id = str()
        fire_properties, ceiling = self._cfast_extract_in()

        super().__init__(config_object, fire_properties, None)
        self.ceiling = ceiling

    @staticmethod
    def _unify_characters(l: str):
        for char in (',', '=', '/'):
            l = l.replace(char, ' ')
        return l

    @staticmethod
    def _get_fire_xy_location(fireline: str):
        return [float(i) for i in fireline.split('LOCATION')[1].split()[:2]]

    @staticmethod
    def _get_fire_compa_ids(fireline: str):
        fire_id = fireline.split('FIRE_ID')[1].split()[0][1:-1]
        compa_id = fireline.split('COMP_ID')[1].split()[0][1:-1]
        return fire_id, compa_id

    @staticmethod
    def _get_hrr_height_area_record(tablline: str):
        data =  tablline.split('DATA')[1].split()
        return tuple(float(i) for i in data[1:4])

    @staticmethod
    def _get_room_height(compaline: str):
        compa_id = compaline.split(' ID')[1].split()[0][1:-1]    # space to avoid i.e. WALL_MATL_ID
        height = float(compaline.split('HEIGHT')[1].split()[0])
        baseline_z =  float(compaline.split('ORIGIN')[1].split()[2])
        return compa_id, (baseline_z, height)

    @staticmethod
    def _calculate_hrrpua(hrr_area_tabl: list):
        hrrs, areas = zip(*hrr_area_tabl)
        return round(max(hrrs) / max(areas), 1)

    ''' data extracted from chid.in '''
    def _cfast_extract_in(self):
        cfast_heights = dict()
        cfast_fire_location = list()
        fire_base_z_tbl = list()
        cfast_hrrpua = float()
        fire_data = list()
        fire_id = str()
        fire_compa = str()

        with open(os.path.join(self.cfast_dir, f'{self.cfast_chid}.in')) as file:
            lines = file.readlines()

        for line in lines:
            if line.startswith('&COMP'):
                line = self._unify_characters(line)
                k, v = self._get_room_height(line)
                cfast_heights[k] = v

            elif not cfast_fire_location and line.startswith('&FIRE'):
                line = self._unify_characters(line)
                cfast_fire_location = self._get_fire_xy_location(line)
                fire_id, fire_compa = self._get_fire_compa_ids(line)

            elif not cfast_hrrpua and line.startswith('&TABL') and 'LABELS' not in line and fire_id in line:
                line = self._unify_characters(line)
                hrr, height, area = self._get_hrr_height_area_record(line)
                fire_data.append((hrr, area))
                fire_base_z_tbl.append(height)

        cfast_fire_location.append(fmean(fire_base_z_tbl) + cfast_heights[fire_compa][0])

        return (None, self._calculate_hrrpua(fire_data), cfast_fire_location), sum(cfast_heights[fire_compa])

    def create_fire_curve(self):
        # chid_compartments.csv
        # translate cfast fire curve to mcsteel
        # [hrr_tab, diameter_tab]
        df = read_csv(os.path.join(self.cfast_dir, f'{self.cfast_chid}_compartments.csv'))
        df = df.iloc[3:][['Time', 'HRR_1']].astype(float)

        def hrr2diam(hrr): return round(2 * np.sqrt(hrr / self.hrrpua / 1e3 / np.pi), 2)
        df['diam'] = df['HRR_1'].apply(hrr2diam)

        self.fire_curve = [i.values.tolist() for i in (df[['Time', 'HRR_1']], df[['Time', 'diam']])]


class Multisimulation:
    def __init__(self, config_object: Config):
        self.config = config_object
        self.structure = self._read_dxf()
        self.data_frame = df(columns=['cfast_fire_id', 'calc_no', 'time', 'x_f', 'y_f', 'z_f', 'x_s', 'y_s', 'z_s',
                                      'distance', 'ceiling_lvl', 'profile', 'u_x', 'u_y', 'u_z', 'HRRPUA', 'alpha'])

    # read dxf geometry
    def _read_dxf(self):
        t1 = sec()

        out(outpth, 'Reading DXF geometry...\r')
        dxffile = dxf.readfile(os.path.join(self.config.config_path, f'{self.config.title}.dxf'))
        out(outpth, f'[OK] DXF geometry imported ({round(sec() - t1, 2)} s)')

        beams = []
        columns = []
        x = 0
        t = len(dxffile.entities)
        # assign LINES elements to columns or beams tables
        start = sec()
        for ent in dxffile.entities:
            progress_bar('Converting lines', x, t)
            if ent.dxftype == 'LINE':
                if ent.start[2] == ent.end[2]:
                    beams.append(ent)
                else:
                    columns.append(ent)
            x += 1
        out(outpth, f'[OK] Lines converted ({round(sec() - start, 2)} s)                       ')

        # assign 3DFACE elements to shells table
        shells = [ent for ent in dxffile.entities if ent.dxftype == '3DFACE']

        return {'b': beams, 'c': columns, 's': shells, 'f': []}

    # append DataFrame to CSV file
    def _writedf2csv(self, iteration_no: int):
        path = os.path.join(self.config.results_path, f'{self.config.title}_set.csv')
        try:
            with open(path):
                header = False
            to_be_written = self.data_frame[-iteration_no:]
        except FileNotFoundError:
            header = True
            to_be_written = self.data_frame

        to_be_written.to_csv(path_or_buf=path, mode='a', header=header)

    def prepare(self):
        prev_s_no_max = 0
        if os.path.exists(self.config.results_path):
            for d in os.listdir(self.config.results_path):
                try: 
                    temp_s_no = int(d.split('_')[0])
                except:
                    continue
                if temp_s_no > prev_s_no_max:
                    prev_s_no_max = temp_s_no
        else:
            os.makedirs(self.config.results_path)

        gen = MCGenerator(self.config)  # current set
        gen.sampling()

        t = sec()
        dfindex = 0
        save_interval = max([int(gen.n / 20), 2])  # save 20 times
        for s_no, scenario in enumerate(gen.set):
            progress_bar('Preparing files', s_no, gen.n)
            s_no += prev_s_no_max + 1
            scenario.create_fire_curve()
            scenario.map(self.structure)
            for i_no, data in enumerate(scenario.mapped):
                i = Iteration(scenario, data, f'{s_no}_{i_no}')
                i.prepare_files()

                # save this iteration data to the data frame
                self.data_frame.loc[dfindex+prev_s_no_max*2] = [s_no, i_no, ctime(sec())] + data
                dfindex += 1

            if ((s_no-prev_s_no_max) % save_interval == 0) or (dfindex == gen.n*2):
                self._writedf2csv(save_interval*2)   # two iterations for each scenario
        out(outpth, f'[OK] {dfindex} file sets were prepared ({round(sec() - t, 2)}) s           ')

        return self.data_frame


if __name__ == '__main__':
    outpth = './mc.log'
    out(outpth, '========================================================================================' 
                      '\nmc.py  Copyright (C) 2022  Kowalski W.'
                      '\nThis program comes with ABSOLUTELY NO WARRANTY.'
                      '\nThis is free software, and you are welcome to redistribute it under certain conditions.'
                      '\nSee GPLv3.0 for details (https://www.gnu.org/licenses/gpl-3.0.html).'
                      '\n========================================================================================\n')

    cfg = Config(sys.argv[1])

    set_of_simulations = Multisimulation(cfg)
    set_of_simulations.prepare()

    out(outpth, '========================================================================================'
                      '\nThank you for using mcsteel package :)'
                      '\nVisit project GitHub site: https://github.com/kowalskiw/mcsteel and contribute!'
                      '\n========================================================================================\n')
