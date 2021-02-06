from os import listdir, chdir, getcwd
import numpy as np
from selsim import linear_inter
from pandas import read_csv, DataFrame
from fdsafir import run_safir
from sys import argv
import export
from pandas.errors import EmptyDataError


# return user configuration directory
def user_config(user_file):
    user = {}
    with open(user_file) as file:
        for line in file.readlines():
            splited = line.split()
            try:
                value = float(splited[-1])
            except ValueError:
                value = splited[-1]
            user[splited[0]] = value
    return user


'''Run single simulation'''


class RunSim:
    # run SAFIR T2D
    def t2d(self, chid):
        dir_content = listdir()
        [dir_content.remove(i) for i in ['{}.in'.format(chid), 'locafi.txt']]
        run_safir(dir_content[0], mcsteel=True)

    # calculate mean temperature of the profile return mean temp - time table
    def mean_temp(self):
        temp_table = []
        nfiber = 0

        for i in range(1, 3):
            temp = 0
            t = 0
            section_temp = []

            print(getcwd())
            with open('b00001_{}.tem'.format(i)) as file:
                tem = file.readlines()

            for line in tem:
                # set number of elements
                if line.startswith('  NFIBERBEAM'):
                    nfiber = int(line.split()[1])

                # set time step value and reset temperature
                elif line.startswith(' TIME'):
                    temp = 0
                    t = float(line.split()[1])

                # try to save previous step mean cross-section temperature
                elif line.startswith('\n'):
                    try:
                        section_temp.append([t, temp/nfiber])
                    except UnboundLocalError:
                        pass

                # add temperature in element in certain step
                else:
                    try:
                        temp += float(line.split()[-1])
                    except (IndexError, UnboundLocalError, ValueError):
                        pass

            temp_table.append(np.array(section_temp))

        return (temp_table[0] + temp_table[1]) / 2


'''Run node queue'''


class Queue:
    def __init__(self, user_file):
        self.user = user_config(user_file)  # import simulation configuration
        chdir(self.user['results_path'])
        self.set = read_csv('{}_set.csv'.format(self.user['case_title']))       # import queue
        try:
            self.results_df = read_csv('{}_results.csv'.format(self.user['case_title']))
            del self.results_df['Unnamed: 0']   # error in importing via rcsv
        except (FileNotFoundError, EmptyDataError):
            self.results_df = DataFrame(columns=['ID', 'temp_max', 'time_crit'])
        self.rs = RunSim

    # run simulation queue
    def run(self):
        results = {}
        for index, row in self.set.iterrows():
            chid = str(row['ID'])
            if int(chid) in self.results_df.ID.tolist():
                continue  # check if queue element is in results DF

            # run simulation and add section temperature curve to the list
            chdir(chid)
            self.rs().t2d(chid)
            results[chid] = self.rs().mean_temp()
            chdir('..')

            # save results every 4 sim
            if index % 4 == 0:
                self.save_res(results, export.temp_crit(self.user['miu']))

        self.save_res(results, export.temp_crit(self.user['miu']))

    # choose theta_a,max and t_theta,a,cr and add those values to case.res
    def save_res(self, tables, t_crit):
        path_to_res = '{}_results.csv'.format(self.user['case_title'])
        if path_to_res not in listdir('.'):
            with open(path_to_res, 'w+'):
                print('{} created'.format(path_to_res))

        for k, v in tables.items():
            temp_max = 0
            upper_index = 0

            for step in v:
                if step[1] > temp_max: temp_max = step[1]       # find max temperature
                if step[1] > t_crit: upper_index = np.where(v == step)    # find point, where theta_crit were exceeded

            # try to interpolate
            try:
                time_crit = linear_inter(v[upper_index - 1], v[upper_index], t_crit)
            except TypeError:
                time_crit = 0

            self.results_df.loc[int(k)] = [str(k), temp_max, time_crit]

        self.results_df.to_csv(path_to_res)

        # export results to txt summary file
        export.summary(self.results_df, export.temp_crit(self.user['miu']), self.user['RSET'])


'''Run multisimulation on cluster'''


class Cluster:
    def __init__(self):
        # import cluster configuration
        pass

    # send task to the nodes' queues
    def assign(self):
        pass

    # initialize calculations on nodes
    def wake_nodes(self):
        pass

    # check stop condition (sim number or RMSE), stop analysis or generate more scenarios if needed
    def check_stop(self):
        pass


if __name__ == '__main__':
    Queue(argv[1]).run()
