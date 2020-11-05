"""this script will itarate over McOZone results in order to choose simulations for SAFIR simulations

input: number of SAFIR simulations, path to simulation folder
output: fire scenarios (LOCAFI.txt), coordinates of fire

proposed directories tree:
simX/config /frame.gid
            /section1.gid
            /section2.gid ______already calculated with ISO curve
            /simX.cel
            /simX.ful
            /simX.geom
            /simX.mat
            /simX.op
            /simX.par
            /simX.str
            /simX.xel (to be changed to DXF file)
            /node.user
    /results/sim1544xxxx/frame.gid
                        /section1.gid
                        /section2.gid
                        /LOCAFI.txt
            /sim1544yyyy/frame.gid
                        /section1.gid
                        /section2.gid
                        /LOCAFI.txt
            .
            /details
            /stoch_res.csv
            .
            /results.txt
The script also will run chosen scenarios (basic McSAFIR functionality)"""

import pandas as pd
from os import mkdir, scandir, getcwd
from core import core
from shutil import copytree as copy
import fdsafir
import numpy as np

'''open data, import to DataFrame'''
# cwd has to be directory where results CSV is stored
# import to pandas DataFrame

'''iterate over DF to select X simulations or 1% of the worst'''


# sort DF with time of critical temperature
# save 1.percentyl or X simulations (if defined) to another DF
# generate LOCAFI.txt files for each record in chosen DF and save them in different dirs


# converting floated string list to numpy array
def np_flt_tab(tab):
    for i in range(len(tab)):
        tab[i] = float(tab[i][:-1])
    return np.array(tab)


# removing '/n' symbol from the end of each string list (tab) elements
def remove_n(tab):
    for i in range(len(tab)):
        tab[i] = tab[i][:-1]
    return tab


# linear interpolation between two points with given first coordinate x_i, returns y_i
def linear_inter(point1, point2, x_i):
    return point1[1] + (x_i - point1[0]) / (point2[0] - point1[0]) * (point2[1] - point1[1])


class Prepare:
    def __init__(self, sim_number=0):
        self.results = pd.read_csv('results\stoch_rest.csv', sep=',')
        self.sim_number = sim_number

    def percentile(self):
        if self.sim_number == 0:
            len_df = 1000
            perc = int(0.01 * len_df)
        else:
            perc = self.sim_number

        sorted = self.results.sort_values(by=['time_crit'])
        for i, r in sorted.iterrows():
            if r['time_crit'] == 0:
                sorted = sorted.drop([i], axis='index')
                continue
            else:
                chosen = sorted.head(n=perc)
                break

        return chosen

    def data_import(self, df_row):
        # open PRI file for single sim
        with open('results\details\{}{}'.format(str(df_row['ID']), '.ozn')) as file:
            ozn = file.readlines()

        fire_row = ozn.index('Localised\n')

        # import ozone time table [min] -> [s]
        time_table = np_flt_tab(remove_n(ozn[fire_row + 9: fire_row + int(float(ozn[fire_row + 8][:-1])) * 2][::2]))*60

        # import HRR table [MW] -> [W] and convert to floats
        hrr_float = np_flt_tab(ozn[fire_row + 9: fire_row + int(float(ozn[fire_row + 8][:-1])) * 2][1::2])*1e6

        # calculate fire area [m] table
        max_diam = float(ozn[fire_row + 5][:-1])
        hrrpua = max(hrr_float[1::2]) / max_diam
        diam = np.array(hrr_float) / hrrpua

        # import z_ceiling
        z_ceil = float(ozn[fire_row + 1][:-1])

        # import absolute fire position
        # make list from DF
        position = ''
        for i in list(df_row[['abs_x', 'abs_y', 'abs_z']]):
            position = '    '.join([position, str(i)])

        return position, z_ceil, self.convert(time_table, diam), self.convert(time_table, hrr_float)

    # convert data [HRR(t) and D(t)] from OZone format (max 120 records) to SAFIR format (max 108 or 20 recs)
    def convert(self, ozn_time, ozn_data, records=20):
        safir_data = []

        # create default 20-records time list
        step = ozn_time[-1] / (records - 1)
        converted_time = []
        [converted_time.append(step * i) for i in range(records)]

        # linear interpolation of values accordingly to converted time list
        inter_values = []
        for t in converted_time:
            for i in range(len(ozn_time)):
                if ozn_time[i] <= t:
                    inter_values.append(linear_inter((ozn_time[i], ozn_data[i]), (ozn_time[i + 1], ozn_data[i + 1]), t))
                    break
                else:
                    raise OverflowError('There is something wrong with interpolation module')

        # add interpolated values to converted time list
        [safir_data.append('    {}    {}\n'.format(converted_time[i], inter_values[i])) for i in range(records)]

        return safir_data

    def gen_lcf(self, position, z_ceil, diam, hrr):
        # add data to LOCAFI.txt core
        lcf = core.copy()

        def ins(arg, index): return lcf[index].split('*')[0] + str(arg) + core[index].split('*')[1]

        lcf[2] = ins(position, 2)
        lcf[3] = ins(z_ceil, 3)

        # add tables of hrr and diameter

        def add(start, tab): [lcf.insert(i + start, tab[i]) for i in range(len(tab))]

        add(lcf.index('DIAMETER\n') + 2, diam)
        add(lcf.index('RHR\n') + 2, hrr)

        # return LOCAFI.txt as list
        # print(lcf)

        return lcf

    # create simulations' directories and save LOCAFI.txt in proper dir
    def save_in_dir(self):
        chosen = self.percentile()
        for index, row in chosen.iterrows():
            # create directory
            try:
                mkdir('results\sim{}'.format(row['ID']))
            except FileExistsError:
                print('Directory sim{} already exists'.format(row['ID']))

            # generate locafi.txt
            # save locafi.txt to the dir
            with open('results\{}\LOCAFI.txt'.format('sim{}'.format(str(row['ID']))), 'w') as file:
                file.writelines(self.gen_lcf(*self.data_import(row)))

        return 0


'''run chosen scenarios in SAFIR'''


def main():
    Prepare().save_in_dir()
    # copy general model files (already calculated for ISO curve) to every SAFIR simulation dir
    # run SAFIR simulation using general model and selected fire
    # return 0/1
    for dir in scandir('results'):
        if dir.is_dir() and dir.name.startswith('sim'):
            # copy general model files (already calculated for ISO curve) to every SAFIR simulation dir
            print(dir.path)
            # for i in scandir('config'):

            try:
                [copy(gid.path, '{}\{}'.format(dir.path, gid.name)) for gid in scandir('config')
                 if gid.is_dir() and gid.name.endswith('.gid')]
            except FileExistsError:
                print('SAFIR-GiD files have been already copied ({})'.format(dir.path))
            # run SAFIR simulation using general model and selected fire
            fdsafir.main('LCF', dir.path)


main()
