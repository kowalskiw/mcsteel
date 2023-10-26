from os import scandir, listdir, chdir
import numpy as np
from pandas import read_csv, DataFrame, Series
import sys
from postprocess import temp_crit, summary
from pandas.errors import EmptyDataError
import os
from time import time as sec
from utils import Mechanical, run_safir, out, Config
from datetime import timedelta as dt

global outpth


# linear interpolation between two points with given first coordinate x_i, returns y_i
def linear_inter(point1, point2, x_i):
    return point1[1] + (x_i - point1[0]) / (point2[0] - point1[0]) * (point2[1] - point1[1])


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


class Single:
    def __init__(self, scenario_data: Series, config: Config):
        self.fire_id = scenario_data.fire_id
        self.calc_no = scenario_data.calc_no
        self.chid = f'{self.fire_id}_{self.calc_no}'
        self.dir_path = os.path.join(config.results_path, self.chid)
        self.safir = config.safir_path
        self.data = scenario_data
        self.config_path = config.config_path


class ThermalOnly(Single):
    def calculate(self):
        start_dir = os.getcwd()

        for i in scandir(self.dir_path):
            if (not i.name.endswith('.in')) or i.name == '{}.in'.format(self.chid):
                continue
            run_safir(i.path, safir_exe_path=self.safir, fix_rlx=False)
            break

        os.chdir(start_dir)


class Full(Single):
    def calculate(self):
        start = sec()
        m = Mechanical(self.dir_path)

        # run thermal analyses
        m.make_thermals(self.config_path)

        for t in m.thermals:
            st = sec()
            t.change_in(m.chid)
            t.run(self.safir)
            print(f'Runtime of "{t.chid}" thermal analysis: {dt(seconds=int(sec() - st))}\n')

        # run mechanical analysis
        st = sec()
        m.change_in()
        m.run(self.safir)
        print(f'Runtime of "{m.chid}" mechanical analysis: {dt(seconds=int(sec() - st))}\n')

        print(f'Summary "{m.chid}" runtime: {dt(seconds=int(sec() - start))}\n')


class Queue:
    def __init__(self, cfg: Config, scenarios_data_frame=None, mode=None):
        if scenarios_data_frame:
            self.df = scenarios_data_frame
        else:
            self.df = read_csv(os.path.join(cfg.results_path, f'{cfg.title}_set.csv'))
        self.mode = mode
        self.config = cfg
        self.queue = []
        self.status = -2
        try:
            self.results_df = read_csv('{}_results.csv'.format(self.config.title))
            try:
                del self.results_df['Unnamed: 0']  # error in importing via rcsv
            except:
                pass
        except (FileNotFoundError, EmptyDataError):
            self.results_df = DataFrame(columns=['ID', 'temp_max', 'time_crit'])
        self.to_compare = {}

    def load_set(self):
        for calc in self.df.iterrows():
            if self.mode == 'full':
                self.queue.append(Full(calc[1], self.config))
            elif not self.mode:
                self.queue.append(ThermalOnly(calc[1], self.config))

    # append DataFrame to CSV file
    def writedf2csv(self, iteration_no):
        path = os.path.join(self.config.results_path, f'{self.config.title}_results.csv')
        try:
            with open(path):
                header = False
            to_be_written = self.results_df[iteration_no:]
        except FileNotFoundError:
            header = True
            to_be_written = self.results_df

        to_be_written.to_csv(path_or_buf=path, mode='a', header=header)

    def set_status(self, rmses, iter_number, rmse_limit=0.001):
        if all([e < rmse_limit for e in rmses]):
            self.status = 2     # to be finished due to RMSE limit
        elif iter_number >= self.config.max_iterations:
            self.status = 1     # to be finished due to number of iterations limit
        else:
            self.status = -1     # to be continued

    # find the results of the fire scenario
    def compare(self):
        # find results of calculation from temperature tab
        def results(temp_tab):
            if type(temp_tab) == int:  # no element above the fire - max temp = ambient
                return temp_tab, 0

            # find time of achieving critical temperature
            temp_max = 0
            upper_index = None
            for step in temp_tab:
                temp_max = step[1] if step[1] > temp_max else temp_max  # find max temperature
                if step[1] >= theta_a_cr and not upper_index:  # find point, where theta_a_cr were exceeded
                    upper_index = np.where(temp_tab == step)[0][0]
            # interpolate to find exact time of achieving theta_a_cr
            try:
                pt1 = list(temp_tab[upper_index - 1])
                pt2 = list(temp_tab[upper_index])
                [i.reverse() for i in [pt1, pt2]]
                time_crit = linear_inter(pt1, pt2, theta_a_cr)
            except TypeError:
                time_crit = 0

            return temp_max, time_crit

        theta_a_cr = temp_crit(self.config.miu)

        # create list of data to compare
        fire_res = [[str(chid), *results(temp_tab)] for chid, temp_tab in self.to_compare.items()]

        # find the worst section (the most heated one) to use it as a result of the fire scenario
        chids, temps, times = zip(*fire_res)
        chids_no0, temps_no0, times_no0 = [], [], []
        for i, t in enumerate(list(times).copy()):
            if times != 0:
                [tab1.append(tab2[i]) for tab1, tab2 in [(chids_no0, chids), (temps_no0, temps), (times_no0, times)]]
        if times_no0:
            the_worst = [j[times_no0.index(min(times_no0))] for j in [chids_no0, temps_no0, times_no0]]
        else:
            the_worst = [j[temps.index(max(temps))] for j in [chids, temps, times]]

        return the_worst

    def run(self):
        def err(id: str, error_type: str):
            self.to_compare.clear()
            errors.append(f'{id},{error_type}\n')

        current_fire_id = int
        errors = ['ID,error_type\n']

        for n, i in enumerate(self.queue):
            # check if scenario is in results
            if i.fire_id in self.results_df.ID.tolist():
                print('    Scenario {} has been already calculated. Continuing...'.format(i.fire_id))
                continue

            print('\n {} calculations started...'.format(i.chid))

            # check if simulation is a part of an error scenario
            try:
                for e in errors:
                    esplt = e.split(',')
                    if esplt[0] == i.fire_id:
                        print(f'[ERROR] Scenario {i.fire_id} is affected by {esplt[1][:-1]}. Passing {i.chid}.')
                        raise AttributeError
            except AttributeError:
                continue

            # save results of fire scenario to the Data Frame
            if i.fire_id != current_fire_id and len(self.to_compare) == 2 and not self.mode:
                self.results_df.loc[len(self.results_df)] = self.compare()
                current_fire_id = i.fire_id

            # run calculations
            chdir(i.dir_path)
            self.to_compare[i.chid] = i.calculate()

            # if no element above the fire base
            dirs = listdir()
            if '{}.err'.format(i.chid) in dirs and not self.mode:
                self.to_compare[i.chid] = 20
                print('[OK] {} calculations finished'.format(i.chid))
                continue

            #i_results = i.calculate()
            # FatalSafirError
            for d in dirs:
                if d.endswith('err'):
                    err(i.fire_id, 'FatalSafirError')  # remove scenario (both syms & add them to err.csv)
                    print('[ERROR] Scenario {} finished with FatalSafirError'.format(i.chid))

#            if not i_results:
#                breakpoint()
#                err(i.fire_id, 'FiberTemperatureError')  # remove scenario (both syms & add them to err.csv)
#                print(f'[ERROR] Scenario {i.chid} finished with FiberTemperatureError')
#                continue

            if not self.mode:
                self.to_compare[i.chid] = mean_temp(20)
                print(f'[OK] {i.chid} calculations finished')

                # save results to files every 5% of progress
                print(bool(int((n+1)) % int(len(self.queue) / 50) == 0))
                if int((n+1)) % int(len(self.queue)  / 50) == 0 and n != 0:
                    self.writedf2csv(n)
                    self.set_status(summary(self.results_df, temp_crit(self.config.miu), self.config.RSET,
                                            savepath=self.config.results_path), len(self.results_df.index))

                if self.status >= 0:
                    return self.status

        return self.status


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
    outpth = './multi.log'
    print(out(outpth, '========================================================================================'
                      '\nmulti.py  Copyright (C) 2022  Kowalski W.'
                      '\nThis program comes with ABSOLUTELY NO WARRANTY.'
                      '\nThis is free software, and you are welcome to redistribute it under certain conditions.'
                      '\nSee GPLv3.0 for details (https://www.gnu.org/licenses/gpl-3.0.html).'
                      '\n========================================================================================\n'))

    cfg = Config(sys.argv[1])

    q = Queue(cfg)
    q.load_set()
    q.run()

    print(out(outpth, '========================================================================================'
                      '\nThank you for using mcsteel package :)'
                      '\nVisit project GitHub site: https://github.com/kowalskiw/mcsteel and contribute!'
                      '\n========================================================================================\n'))
