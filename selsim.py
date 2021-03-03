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
from os import makedirs, scandir, chdir
from shutil import copyfile as copy
from math import ceil
from sys import argv

from fdsafir import scripted, user_config


'''iterate over DF to select X simulations or 1% of the worst'''


# sort DF with time of critical temperature
# save 1.percentyl or X simulations (if defined) to another DF
# generate LOCAFI.txt files for each record in chosen DF and save them in different dirs


class Prepare:
    def __init__(self, results_path, title, sim_number=0):
        chdir(results_path)
        self.results = pd.read_csv('{}_results.csv'.format(title), sep=',')
        self.sim_number = sim_number

    def percentile(self):
        # define how many scenarios will be considered in FEM analyses (default 10%)
        if self.sim_number == 0:
            perc = ceil(0.01 * len(self.results.index))
        else:
            perc = self.sim_number

        # sort values by time to exceed temp crit
        copied_df = self.results.sort_values(by=['time_crit'], ascending=True)
        not_exceeded = pd.DataFrame()
        chosen = pd.DataFrame()

        for i, r in copied_df.iterrows():
            if int(r['time_crit']) == 0:     # drop rows where time crit haven't been exceeded
                not_exceeded = not_exceeded.append(copied_df.loc[i])
                copied_df = copied_df.drop([i])
                continue
            elif int(r['time_crit']) != 0:
                chosen = copied_df.head(n=perc)
                break
        try:
            diff = perc - len(chosen.index)
        except TypeError:
            diff = perc

        # when there is no scenario where temp crit have been exceeded
        if diff > 0:
            not_exceeded = not_exceeded.sort_values(by=['temp_max'], ascending=False)
            chosen = chosen.append(not_exceeded.head(n=diff))

        return chosen

    def save_in_dir(self, config_path):
        chosen = self.percentile()

        for index, row in chosen.iterrows():
            # copy locafi.txt and essential SAFIR files
            id = int(row['ID'])

            try:
                makedirs('worst\{}'.format(id))
            except FileExistsError:
                pass

            copy('{}\locafi.txt'.format(id), 'worst\{}\locafi.txt'.format(id))
            for i in scandir(config_path):
                if i.name.endswith('.gid') and i.is_dir():
                    gid = i.name[:-4]
                    copy('{0}\{1}.gid\{1}.in'.format(config_path, gid), 'worst\{}\{}.in'.format(id, gid))

            print('simulation {} prepared'.format(id))

        return 0


'''run chosen scenarios in SAFIR'''


def run(config, sim_number=0):
    # select the worst scenarios and generate input files
    Prepare(config['results_path'], config['case_title'], sim_number=sim_number).save_in_dir(config['config_path'])
    # run SAFIR simulation using global structure model and selected fire
    scripted(config['safir_path'], config['config_path'], config['results_path'])


if __name__ == '__main__':
    run(user_config(argv[1]))
