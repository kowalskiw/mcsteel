"""this script will itarate over McOZone results in order to choose simulations for SAFIR simulations

input: number of SAFIR simulations, path to results CSV folder
output: fire scenarios (LOCAFI.txt), coordinates of fire

The script also will run chosen scenarios (basic McSAFIR functionality)"""

import pandas as pd
import numpy as np

'''open data, import to DataFrame'''
# cwd has to be directory where results CSV is stored
# import to pandas DataFrame




'''iterate over DF to select X simulations or 1% of the worst'''
# sort DF with time of critical temperature
# save 1.percentyl or X simulations (if defined) to another DF
# generate LOCAFI.txt files for each record in chosen DF and save them in different dirs


class Choose:
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

    def gen_lcf(self, df):
        # open PRI file for this sim
        with open('details\{}{}'.format(str(df['name']), '.ozn')) as file:
            ozn = file.readlines()

        fire_row = ozn.index('Localised\n')

        # import HRR table
        hrr = ozn[fire_row+9: fire_row + int(ozn[fire_row+8][:-1])*2]
        hrr_float = []
        [hrr_float.append(float(hrr[i][:-1])) for i in range(len(hrr))]

        # calculate fire area table
        max_area = 3.1415 * float(ozn[fire_row+5][:-1]) ** 2 / 4
        hrrpua = max(hrr[1::2]) / max_area
        area = hrr[0::2]
        for i in range(1, len(hrr), step=2):
            area.insert(i, hrr[i] / hrrpua)

        # add data to LOCAFI.txt core


        # return LOCAFI.txt as list

        pass


'''run chosen scenarios in SAFIR'''
# copy general model files (already calculated for ISO curve) to every SAFIR simulation dir
# run SAFIR simulation using general model and selected fire
# return 0/1


def runSAFIRs():
    pass