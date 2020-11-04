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

        def flt_tab(tab):
            for i in range(len(tab)):
                tab[i] = float(tab[i][:-1])
            return tab

        def remove_n(tab):
            for i in range(len(tab)):
                tab[i] = tab[i][:-1]
            return tab

        def data_to_time(time_list, data_list):
            curve = []
            [curve.append('    {}    {}\n'.format(time_list[i], data_list[i])) for i in range(len(data_list))]
            return curve

        time_table = remove_n(ozn[fire_row + 9: fire_row + int(float(ozn[fire_row + 8][:-1])) * 2][::2])

        # import HRR table and convert to floats
        hrr_float = flt_tab(ozn[fire_row+9: fire_row + int(float(ozn[fire_row+8][:-1]))*2][1::2])

        # calculate fire area table
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

        return position, z_ceil, data_to_time(time_table, diam), data_to_time(time_table, hrr_float)

    def gen_lcf(self, position, z_ceil, diam, hrr):
        # add data to LOCAFI.txt core
        lcf = core.copy()
        def ins(arg, index): return lcf[index].split('*')[0] + str(arg) + core[index].split('*')[1]

        lcf[2] = ins(position, 2)
        lcf[3] = ins(z_ceil, 3)

        # add tables of hrr and diameter

        def add(start, tab): [lcf.insert(i+start, tab[i]) for i in range(len(tab))]

        add(lcf.index('DIAMETER\n')+2, diam)
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
