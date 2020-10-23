"""this script will itarate over McOZone results in order to choose simulations for SAFIR simulations

input: number of SAFIR simulations, path to results CSV folder
output: fire scenarios (LOCAFI.txt), coordinates of fire

The script also will run chosen scenarios (basic McSAFIR functionality)"""

import pandas as pd
from os import mkdir, scandir
from core import core
from shutil import copy
import fdsafir

'''open data, import to DataFrame'''
# cwd has to be directory where results CSV is stored
# import to pandas DataFrame

'''iterate over DF to select X simulations or 1% of the worst'''
# sort DF with time of critical temperature
# save 1.percentyl or X simulations (if defined) to another DF
# generate LOCAFI.txt files for each record in chosen DF and save them in different dirs


class Prepare:
    def __init__(self, sim_number=0):
        self.results = pd.read_csv('stoch_rest.csv', sep=',')
        self.sim_number = sim_number

    def percentile(self):
        if self.sim_number == 0:
            len_df = 1000
            perc = int(0.01 * len_df)
        else:
            perc = self.sim_number

        self.results.sort_values(by=['time_crit'])

        for i, r in self.results.iterrows():
            if r['time_crit'] == 0:
                self.results.drop(i)
            else:
                chosen = self.results.loc[0:perc]
                break
            raise TabError("There is an error with results DataFrame")

        return chosen

    def data_import(self, df_row):
        # open PRI file for single sim
        with open('details\{}{}'.format(str(df_row['name']), '.ozn')) as file:
            ozn = file.readlines()

        fire_row = ozn.index('Localised\n')

        # import HRR table
        hrr = ozn[fire_row+9: fire_row + int(ozn[fire_row+8][:-1])*2]
        hrr_float = []
        [hrr_float.append(float(hrr[i][:-1])) for i in range(len(hrr))]

        # calculate fire area table
        max_diam = float(ozn[fire_row + 5][:-1])
        hrrpua = max(hrr[1::2]) / max_diam
        diam = hrr[0::2]
        for i in range(1, len(hrr), step=2):
            diam.insert(i, hrr[i] / hrrpua)

        # import z_ceiling
        z_ceil = float(ozn[fire_row + 1][:-1])

        # import absolute fire position
        position = df_row[['abs_x', 'abs_y', 'abs_z']]

        return position, z_ceil, diam, hrr

    def gen_lcf(self, position, z_ceil, diam, hrr):
        # add data to LOCAFI.txt core
        def ins(arg, index): return core[index].split('*')[0] + str(arg) + core[index].split('*')[1]
        core[2] = ins(position, 2)
        core[3] = ins(z_ceil, 3)

        # add tables of hrr and diameter

        def add(start, tab): [core.insert(i, tab[i]) for i in range(start, len(tab))]
        add(core.index('DIAMETER\n')+2, diam)
        add(core.index('RHR\n') + 2, diam)

        # return LOCAFI.txt as list
        return core

    # create simulations' directories and save LOCAFI.txt in proper dir
    def save_in_dir(self):
        chosen = self.percentile()
        for index, row in chosen.iterrows():
            # create directory
            mkdir('sim{}'.format(str(row['ID'])))

            # generate locafi.txt
            # save locafi.txt to the dir
            with open('{}\LOCAFI.txt'.format('sim{}'.format(str(row['ID'])))) as file:
                file.writelines(self.gen_lcf(*self.data_import(chosen)))

        return 0


'''run chosen scenarios in SAFIR'''


def main():
    # copy general model files (already calculated for ISO curve) to every SAFIR simulation dir
    # run SAFIR simulation using general model and selected fire
    # return 0/1
    for dir in scandir():
        if dir.is_dir() and dir.startswith('sim'):
            # copy general model files (already calculated for ISO curve) to every SAFIR simulation dir
            [copy(gid, dir) for gid in scandir('gids')]

            # run SAFIR simulation using general model and selected fire
            fdsafir.main('LCF')


main()
