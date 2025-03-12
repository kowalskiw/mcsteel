import os.path
import numpy as np
import pandas as pd
import sys
from os import makedirs, scandir
from time import time as sec
from time import ctime
from shutil import copy2
from statistics import fmean

from mcsteel.utils import ThermalTEM, Config, progress_bar, out, triangular
import mcsteel.core
import mcsteel.fires
import mcsteel.geom

global outpth


'''Fire scenario class - it basically stores scenario data and generates HRR(t) and D(t) curves'''
class FireScenario:
    def __init__(self, config_object: Config, fire_properties, sprinkler_activation):
        self.config = config_object  # config class
        self.fire_location = fire_properties[2]  # [x, y, z]
        self.alpha = fire_properties[0]  # [W/s^2]
        self.hrrpua = fire_properties[1]  # [W/m^2]
        self.sprinklers = sprinkler_activation  # [s]
        self.fire_curve = [[], []]  # [[0,... time steps ... t_end], [HRR(0), ... HRR ... HRR(t_end)]]
        self.fire_type = config_object.fire_type  # fire curve function type form fires.Fires
        self.mapped = []  # complete set of data for profiles to be calculated in this scenario

        self.locafi_lines = []  # lines for locafi.txt fire file
        self.no_valid_elements = False    # if no valid elements are available for this fire (not expected for well set projects)

        self.ceiling = 1e5  # level of ceiling above the fire source (here the space begins!)
        self.profiles = []  # [profile1, profile2, profile3] profile1 = fdsafir2.ThermalTEM

    # cut those between the values
    # def cut_lines(self, lines, structure):
    #     for line in lines:
    #         # start point cannot be higher than end point
    #         if line.start[2] > line.end[2]:
    #             start_rev = line.end
    #             end_rev = line.start
    #             line.start = start_rev
    #             line.end = end_rev
    #
    #         z1 = line.start[2]
    #         z2 = line.end[2]
    #         # do not consider lines beneath the fire base or above the ceiling
    #         if z2 <= self.fire_location[2] or z1 >= self.ceiling:
    #             continue
    #         # accept lines in (fire base, ceiling) ranges
    #         elif z1 > self.fire_location[2] and z2 < self.ceiling:
    #             structure['temp'].append(line)
    #         # cut lines to (fire base, ceiling) ranges with 0.01 tolerance
    #         else:
    #             to_save = None
    #             if z1 <= self.fire_location[2]:
    #                 to_save = line
    #                 to_save.start = (line.start[0], line.start[1], self.fire_location[2] + 0.01)
    #             if z2 >= self.ceiling:
    #                 to_save = line
    #                 to_save.end = (line.end[0], line.end[1], self.ceiling - 0.01)
    #             # check if line has non-zero length
    #             if np.linalg.norm(np.array(to_save.start) - np.array(to_save.end)) > 0:
    #                 structure['temp'].append(to_save)
    #
    # # checking if point consists in polygon (XY plane only)
    # @staticmethod
    # def ray_tracing_method(point: iter, poly: iter) -> bool:
    #     n = len(poly)
    #     inside = False
    #     x = point[0]
    #     y = point[1]
    #
    #     p1x = poly[0][0]
    #     p1y = poly[0][1]
    #     for i in range(n + 1):
    #         p2x = poly[i % n][0]
    #         p2y = poly[i % n][1]
    #         if y > min(p1y, p2y):
    #             if y <= max(p1y, p2y):
    #                 if x <= max(p1x, p2x):
    #                     xints = None
    #                     if p1y != p2y:
    #                         xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
    #                     if p1x == p2x or x <= xints:
    #                         inside = not inside
    #         p1x, p1y = p2x, p2y
    #
    #     return inside
    #
    # # map fire location with structure to find the most heated profiles to be analysed
    def map(self, structure: mcsteel.geom.Geometry):
        self.mapped  = structure.find_closest_sections(self.fire_location).items()
        if not len(self.mapped):
            self.no_valid_elements = True
        else:
            self.ceiling = structure.get_above_shell_lvl(self.fire_location)
            self.profiles = self.mapped.keys()


    #     # (*fire coords, *section coords, length of the fire-section vector, level of shell above the fire, profile,
    #     # *unit vector))
    #
    #     # check for shell (plate, ceiling) existing above the fire assign the level if true
    #     for s in structure['s']:
    #         lvl = s.points[0][2]  # read level from first point of shell
    #         if float(self.fire_location[2]) <= lvl < self.ceiling and self.ray_tracing_method(self.fire_location, s.points):
    #             self.ceiling = lvl
    #
    #     for i, element_type in enumerate(['b', 'c']):
    #         lines = structure[element_type]  # choose beams or columns as lines
    #         self.cut_lines(lines, structure)  # cut beams accordingly to Z in (fire_z - shell_lvl) range and map to relative
    #
    #         self.mapped.append(self.map_lines(element_type, structure))
    #
    #         self.profiles.append(self.mapped[i][8])

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
        return self.fire_curve


'''Single iteration class - necessary to gather data from samplers and create simulation input files'''
class Iteration:
    def __init__(self, fire_scenario: FireScenario, section_coordinates: list, chid: str):
        self.chid = chid
        self.section_name = '_'.join(chid.split('_')[1:])
        self.fire = fire_scenario
        self.config = fire_scenario.config
        self.dir_path = os.path.join(self.config.results_path, chid)

        self.section_coords = section_coordinates

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

    def copy_section(self):
        copy2(self.config.section_path(self.section_name), self.dir_path)
        ThermalTEM(1, [self.section_name, [], []], self.config.config_path, 'lcf',
                   self.config.time_end, self.dir_path).change_in(os.path.basename(self.dir_path))

    def write_dummy_structural(self, unit_v: np.array = np.array([0, 0, 1])):
        # calculate nodes position
        np_section = np.array(self.section_coords).astype(float)
        node1 = np_section - (unit_v / 1000)
        node2 = np_section + (unit_v / 1000)
        center = np_section
        # THIS PART SHOULD BE REDESIGNED - ASSUMING +Z UNIT VECTOR AND +X LAX
        # add perpendicular vector to section point
        if unit_v[0] != 0:
            lax = np_section + np.array([-unit_v[2] / unit_v[0], 0, 1])
        elif unit_v[1] != 0:
            lax = np_section + np.array([0, -unit_v[2] / unit_v[1], 1])
        elif unit_v[2] != 0:
            lax = np_section + np.array([1, 0, -unit_v[0] / unit_v[2]])
        else:
            # unit vector cannot be [0, 0, 0]
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
        # THIS CHECK SHOULD BE RECONSIDERED
        # check if there are any elements above the fire source
        # if self.data[-4] != self.data[-4]:
        #     with open(os.path.join(self.dir_path, f'{self.chid}.err'), 'w') as err:
        #         mess = f'[WARNING] There are no structural elements above the fire base in the' \
        #                f' {self.chid} fire scenario'
        #         err.write(f'{mess}\nMax element temperature in this scenario is equal to the ambient temperature')
        #     out(outpth, mess)

        # create directories
        makedirs(self.dir_path)

        # save fire file to the directory
        self.prepare_locafi()

        # create SAFIR files
        self.copy_section()
        self.write_dummy_structural()


'''Monte Carlo sampler for fire scenarios'''
# TO BE REDESIGNED ACCORDING TO THE NEW ARCHITECTURE
class MCGenerator:
    def __init__(self, config_object: Config, geometry: mcsteel.geom.Geometry):
        self.config = config_object
        self.n = self.config.max_iterations if self.config.max_iterations else 100  # size of the sample
        self.geom = geometry    # gmsh geometry data - fuel volumes
        self.fuel_db = self._read_fuel() if self.config.fire_type.lower() != 'cfast' else None  # fuel properties
        self.set = []

    def _read_fuel(self):
        t0 = sec()
        out(outpth, 'Importing fuel properites...\r')
        self.fuel_db = mcsteel.fires.FuelDB(self.config, autoopen=True)
        out(outpth, f'[OK] Fuel properties imported ({round(sec() - t0, 2)} s)                      ')
        return self.fuel_db

    # find fire localization and properties of fuel in that place
    def _find_fire_origin(self):
        def random_position(xes, yes, zes):
            coordinates = []
            [coordinates.append(np.random.randint(int(10 * i[0]), int(10 * i[1])) / 10) for i in [xes, yes, zes]]
            return coordinates

        corrected_volumes = []
        for vol in self.geom.volumes:
            corrected_volumes.append(abs(np.prod([vol.bbox[3+i]-vol.bbox[i] for i in range(3)])) * vol.fairshare)

        total_cv = sum(corrected_volumes)
        probs = []
        for cv in corrected_volumes:
            probs.append(cv/total_cv)

        fire_volume = self.geom.volumes[np.random.choice(len(probs), p=probs)]

        fire_position = [np.random.uniform(fire_volume.bbox[i], fire_volume.bbox[i+3]) for i in range(3)]

        return self.fuel_db.get_fuel_properties(fire_volume.name), fire_position

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
            return hrrpua * np.random.lognormal(-9.72, 0.97)  # [kW/s2]
        else:
            return triangular(properties.alpha_min, properties.alpha_max, mode=properties.alpha_mode)  # [kW/s2]


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
                except UnboundLocalError as e:
                    out(outpth, f'[ERROR] Internal error')
                    raise Exception(e)

            else:
                self.set.append(FireScenario(self.config, *self._find_fire()))
        out(outpth, f'[OK] {self.n} fire scenarios have been sampled ({round(sec() - t, 3)}) s                  ')
        return self.set


'''Fire scenario class with data obtained from CFAST output'''
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
        df = pd.read_csv(os.path.join(self.cfast_dir, f'{self.cfast_chid}_compartments.csv'))
        df = df.iloc[3:][['Time', 'HRR_1']].astype(float)

        def hrr2diam(hrr): return round(2 * np.sqrt(hrr / self.hrrpua / 1e3 / np.pi), 2)
        df['diam'] = df['HRR_1'].apply(hrr2diam)

        self.fire_curve = [i.values.tolist() for i in (df[['Time', 'HRR_1']], df[['Time', 'diam']])]

        # update time end in config for this iteration
        self.config.time_end = round(df['Time'].max())


'''Top-tier class to manipulate all multisimulation'''
class Multisimulation:
    def __init__(self, config_object: Config):
        self.config = config_object
        self.structure = self._read_geom()
        self.iterations_per_scenario = len(self.structure.sections)
        self.data_frame = pd.DataFrame(columns=['cfast_fire_id', 'calc_no', 'time', 'x_f', 'y_f', 'z_f', 'x_s', 'y_s',
                                                'z_s', 'distance', 'ceiling_lvl', 'profile', 'u_x', 'u_y', 'u_z',
                                                'HRRPUA', 'alpha'])

    # read GEO file with structure lines, floors surfaces and fuel volumes
    def _read_geom(self):
        t1 = sec()
        out(outpth, 'Reading GEO geometry file...\r')
        geofile = os.path.join(self.config.config_path, f'{self.config.title}.geo')
        geometry_instance = mcsteel.geom.read_geo_geom(geofile)
        out(outpth, f'[OK] DXF geometry imported ({round(sec() - t1, 2)} s)')

        return geometry_instance

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

    # some iterations may be already done
    def _find_the_previous_sim_id(self):
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

        return prev_s_no_max

    def prepare(self):
        previous_sim_id = self._find_the_previous_sim_id()

        gen = MCGenerator(self.config, self.structure)  # MC sampler for current set of iterations
        gen.sampling()    # create the set of fire scenarios

        t = sec()
        unsaved_df_records_no = 0
        save_interval = max([int(gen.n / 20), 2])  # save 20 times
        for s_no, scenario in enumerate(gen.set):
            progress_bar('Preparing files', s_no, gen.n)
            s_no += previous_sim_id + 1    # current iteration id
            scenario.create_fire_curve()   # define fire scenario
            scenario.map(self.structure)
            if scenario.no_valid_elements:
                print(f'[WARNING] Excluding scenario {s_no} from the analysis')
                continue
            for section, section_coords in scenario.mapped.items():
                i = Iteration(scenario, section_coords[0], f'{s_no}_{section}')
                i.prepare_files()

                # save this iteration data to the data frame
                i_data = ([s_no, section, ctime(sec())] + scenario.fire_location + section_coords[0] +
                          [section_coords[1]] + [scenario.ceiling, section, 0, 0, 1, scenario.hrrpua, scenario.alpha])
                self.data_frame.loc[len(self.data_frame)] = i_data
                unsaved_df_records_no += 1

            if (s_no-previous_sim_id) % save_interval == 0:
                self._writedf2csv(unsaved_df_records_no)
                unsaved_df_records_no = 0

        if unsaved_df_records_no:
            self._writedf2csv(unsaved_df_records_no)
        out(outpth, f'[OK] {gen.n * self.iterations_per_scenario} file sets were prepared ({round(sec() - t, 2)}) s           ')

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
