from os import scandir, listdir, chdir
import numpy as np
from pandas import read_csv, DataFrame
import sys
import export
from pandas.errors import EmptyDataError

from fdsafir import run_safir, user_config, Logger, progressBar


# linear interpolation between two points with given first coordinate x_i, returns y_i
def linear_inter(point1, point2, x_i):
    return point1[1] + (x_i - point1[0]) / (point2[0] - point1[0]) * (point2[1] - point1[1])


'''Run single simulation'''


# run SAFIR T2D
def t2d(chid, safir_dir_path):
    for i in scandir():
        if (not i.name.endswith('.in')) or i.name == '{}.in'.format(chid):
            continue
        run_safir(i.name, safir=safir_dir_path, mcsteel=True)
        break


def mean_temp(amb_temp):
    temp_table = []
    nfiber = 0
    is_reading_temp = False

    for i in range(1, 3):
        temp = 0
        t = 0
        section_temp = []

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
                is_reading_temp = True

            # try to save previous step mean cross-section temperature
            elif line.startswith('\n'):
                try:
                    section_temp.append([t, temp/nfiber])
                except UnboundLocalError:
                    pass

            # add temperature in element in certain step
            elif is_reading_temp:
                try:
                    fiber_temp = float(line.split()[-1])
                    if fiber_temp >= amb_temp:
                        temp += fiber_temp
                    else:
                        print('[WARNING] SafirError: Fiber temperature is lower than ambient ({} °C < {} °C)'.format(
                            fiber_temp, amb_temp))
                        raise ChildProcessError
                except (IndexError, UnboundLocalError, ValueError):
                    pass

        temp_table.append(np.array(section_temp))

    return (temp_table[0] + temp_table[1]) / 2


class RunSim:

    # calculate mean temperature of the profile return mean temp - time table
    pass


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
                pass
        except (FileNotFoundError, EmptyDataError):
            self.results_df = DataFrame(columns=['ID', 'temp_max', 'time_crit', 'compared'])
        self.rs = RunSim

    # choose theta_a,max and t_theta,a,cr and add those values to case.res
    def save_res(self, tables, t_crit):
        print('Saving results...')
        path_to_res = '{}_results.csv'.format(self.user['case_title'])
        if path_to_res not in listdir('.'):
            with open(path_to_res, 'w+'):
                print('[OK] {} results file created'.format(path_to_res))

        compared = []
        for k, v in tables.items():
            # no element exception
            if type(v) == int:
                temp_max = v
                time_crit = 0

            else:
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
                # smaller time_crit except 0
                if compared[0][2] + compared[1][2] > 0 and compared[0][2] < compared[1][2]:
                    # print('    {} vs {}: Temperature in column scenario is higher -> added to results'.format(
                    #     *[i[0] for i in compared]))
                    comp_id = compared.pop(0)[0]

                elif compared[0][1] > compared[1][1]:   # bigger temp_max
                    # print('   {} vs {}: Temperature in column scenario is higher -> added to results'.format(
                    #     *[i[0] for i in compared]))
                    comp_id = compared.pop(0)[0]

                else:
                    # print('    {} vs {}: Temperature in beam scenario is higher -> added to results'.format(
                    #     *[i[0] for i in compared]))
                    comp_id = compared.pop(1)[0]

                self.results_df.loc[len(self.results_df)] = compared[0] + [comp_id]   # save results to the Data Frame
                compared = []

        self.results_df.to_csv(path_to_res)

        # export results to txt summary file
        export.summary(self.results_df, export.temp_crit(self.user['miu']), self.user['RSET'])

        # # clear Data Frame
        # self.results_df = self.results_df.iloc[0:0]

        return 0

    # run simulation queue
    def run(self):
        ambient_temperature = 20
        results = {}
        errors = ['ID,error_type\n']

        def remove_err(id: str, error_type: str):
            iid = int(id)
            if iid % 2 == 1:
                try:
                    results.pop(str(iid - 1))
                    errors.append('{},{}\n'.format(iid - 1, error_type))
                except KeyError:
                    print('{} not found in results'.format(iid - 1))
                errors.append('{},{}\n'.format(id, error_type))
            else:
                errors.append('{},{}\n'.format(id, error_type))
                errors.append('{},{}\n'.format(iid + 1, error_type))

        x = 0
        l = len(self.set.index)
        for index, row in self.set.iterrows():
            progressBar('Multisimulating in progress ', x, l)
            chid = str(row['ID'])
            # check if queue element is in results DF
            if (int(chid) in self.results_df.ID.tolist()) or (int(chid) in self.results_df.compared.tolist()):
                print('    Simulation {} has been already calculated. Continuing...'.format(chid))
                continue

            print('\nScenario {} calculations started...'.format(chid))

            # check if simulation is a part of error scenario
            try:
                for e in errors:
                    esplt = e.split(',')
                    if esplt[0] == chid:
                        print('[ERROR] Scenario {} is affected by {}'.format(chid, esplt[1][:-1]))
                        raise AttributeError
            except AttributeError:
                continue

            # run simulation and add section temperature curve to the list
            chdir(chid)
            dirs = listdir()

            # ElementNotFoundError
            if '{}.err'.format(chid) in dirs:
                results[chid] = ambient_temperature
                print('[OK] {} calculations finished'.format(chid))

            else:
                t2d(chid, self.user['safir_path'])
                # FatalSafirError
                for i in dirs:
                    if i.endswith('err'):
                        remove_err(chid, 'FatalSafirError')  # remove scenario (both syms & add them to err.csv)
                        print('[ERROR] Scenario {} finished with FatalSafirError'.format(chid))
                try:
                    results[chid] = mean_temp(ambient_temperature)
                    print('[OK] {} calculations finished'.format(chid))
                # FiberTemperatureError
                except ChildProcessError:
                    remove_err(chid, 'FiberTemperatureError')  # remove scenario (both syms & add them to err.csv)
                    print('[ERROR] Scenario {} finished with FiberTemperatureError'.format(chid))

            chdir('..')

            # save results every 5 scenarios (10 sim)
            if (index+1) % 8 == 0:
                self.save_res(results, export.temp_crit(self.user['miu']))
                results.clear()

            # save errors
            if len(errors) > 0:
                with open('{}\err.csv'.format(self.user['results_path']), 'w') as file:
                    file.writelines(errors)

            x += 1

        return [self.save_res(results, export.temp_crit(self.user['miu'])) if results else -1]


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
    sys.stdout = Logger('multi.log')
    Queue(sys.argv[1]).run()
