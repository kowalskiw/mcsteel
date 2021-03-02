from os import listdir, chdir, getcwd
import numpy as np
from pandas import read_csv, DataFrame
from sys import argv
import export
from pandas.errors import EmptyDataError

from fdsafir import run_safir, user_config


# linear interpolation between two points with given first coordinate x_i, returns y_i
def linear_inter(point1, point2, x_i):
    return point1[1] + (x_i - point1[0]) / (point2[0] - point1[0]) * (point2[1] - point1[1])


'''Run single simulation'''


class RunSim:
    # run SAFIR T2D
    def t2d(self, chid, safir_dir_path):
        dir_content = listdir()
        for i in dir_content:
            if (not i.endswith('.in')) or i == '{}.in'.format(chid):
                continue
            run_safir(i, safir=safir_dir_path, mcsteel=True)
            break

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
            try:
                del self.results_df['Unnamed: 0']  # error in importing via rcsv
            except:
                print('nie usunieto unnamed')
        except (FileNotFoundError, EmptyDataError):
            self.results_df = DataFrame(columns=['ID', 'temp_max', 'time_crit', 'compared'])
        self.rs = RunSim

    # run simulation queue
    def run(self):
        results = {}
        for index, row in self.set.iterrows():
            chid = str(row['ID'])
            # check if queue element is in results DF
            if (int(chid) in self.results_df.ID.tolist()) or (int(chid) in self.results_df.compared.tolist()):
                continue

            # run simulation and add section temperature curve to the list
            chdir(chid)
            print('Started {} calculations'.format(chid))
            self.rs().t2d(chid, self.user['safir_path'])
            results[chid] = self.rs().mean_temp()
            chdir('..')

            # save results every 5 scenarios (10 sim)
            if (index+1) % 8 == 0:
                self.save_res(results, export.temp_crit(self.user['miu']))
                results.clear()

        return [self.save_res(results, export.temp_crit(self.user['miu'])) if results else -1]

    # choose theta_a,max and t_theta,a,cr and add those values to case.res
    def save_res(self, tables, t_crit):
        path_to_res = '{}_results.csv'.format(self.user['case_title'])
        if path_to_res not in listdir('.'):
            with open(path_to_res, 'w+'):
                print('{} created'.format(path_to_res))

        compared = []
        for k, v in tables.items():
            temp_max = 0
            upper_index = None

            for step in v:
                if step[1] > temp_max: temp_max = step[1]       # find max temperature
                if step[1] > t_crit and not upper_index:  # find point, where theta_crit were exceeded
                    upper_index = np.where(v == step)[0][0]

            # try to interpolate
            try:
                pt1 = list(v[upper_index - 1])
                pt2 = list(v[upper_index])
                [i.reverse() for i in [pt1, pt2]]
                time_crit = linear_inter(pt1, pt2, t_crit)
            except TypeError:
                time_crit = 0

            compared.append([str(k), temp_max, time_crit])

            # compare beam with column scenario
            if len(compared) == 2:
                print(compared)
                # smaller time_crit except 0
                if compared[0][2] + compared[1][2] > 0 and compared[0][2] < compared[1][2]:
                    comp_id = compared.pop(0)[0]
                    print('usunieta belka')

                elif compared[0][1] > compared[1][1]:   # bigger temp_max
                    comp_id = compared.pop(0)[0]
                    print('usunieta belka')

                else:
                    comp_id = compared.pop(1)[0]
                    print('usunieta kolumna')

                self.results_df.loc[len(self.results_df)] = compared[0] + [comp_id]   # save results to the Data Frame
                compared = []

        self.results_df.to_csv(path_to_res)

        # export results to txt summary file
        export.summary(self.results_df, export.temp_crit(self.user['miu']), self.user['RSET'])

        # # clear Data Frame
        # self.results_df = self.results_df.iloc[0:0]

        return 0


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
