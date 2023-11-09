from os import scandir, chdir, getcwd
from os.path import join
import numpy as np
from sys import argv
from pandas import DataFrame, Series, read_csv, concat
from pandas.errors import EmptyDataError
from time import time as sec
from datetime import timedelta as dt

from utils import Mechanical, run_safir, out, Config
from postprocess import temp_crit, summary

global outpth


# linear interpolation between two points with given first coordinate x_i, returns y_i
def linear_inter(point1, point2, x_i):
    return point1[1] + (x_i - point1[0]) / (point2[0] - point1[0]) * (point2[1] - point1[1])


# calculate timetable of average temperature in element
# adequate for critical temperature (steel only)
def get_temp_table(sim_path, amb_temp=20, ng_beam=2):
    temp_table = []
    nfiber = 0
    is_reading_temp = False

    for i in range(1, ng_beam + 1):
        temp = 0
        t = 0
        section_temp = []

        with open(join(sim_path, 'b00001_{}.tem'.format(i))) as file:
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
                    section_temp.append([t, temp / nfiber])
                except UnboundLocalError:
                    pass

            # add temperature in element in certain step
            elif is_reading_temp:
                try:
                    fiber_temp = float(line.split()[-1])
                    if fiber_temp >= amb_temp:
                        temp += fiber_temp
                    else:
                        out(outpth,
                            '[WARNING] SafirError: Fiber temperature is lower than ambient ({} °C < {} °C)'.format(
                                fiber_temp, amb_temp))
                        raise ChildProcessError
                except (IndexError, UnboundLocalError, ValueError):
                    pass

        temp_table.append(np.array(section_temp))

    return (temp_table[0] + temp_table[1]) / 2


class Single:
    def __init__(self, scenario_data: Series, config: Config):
        self.fire_id = int(scenario_data.fire_id)
        self.calc_no = int(scenario_data.calc_no)
        self.chid = f'{self.fire_id}_{self.calc_no}'
        self.dir_path = join(config.results_path, self.chid)
        self.safir = config.safir_path
        self.data = scenario_data
        self.config_path = config.config_path


class ThermalOnly(Single):
    def calculate(self):
        start_dir = getcwd()
        status = 1
        err = f'No IN file in {self.dir_path}'

        for i in scandir(self.dir_path):
            if (not i.name.endswith('.in')) or i.name == '{}.in'.format(self.chid):
                continue
            status, err = run_safir(i.path, safir_exe_path=self.safir, fix_rlx=False)
            break

        chdir(start_dir)
        return status, err


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
            out(outpth, f'Runtime of "{t.chid}" thermal analysis: {dt(seconds=int(sec() - st))}\n')

        # run mechanical analysis
        st = sec()
        m.change_in()
        m.run(self.safir)
        out(outpth, f'Runtime of "{m.chid}" mechanical analysis: {dt(seconds=int(sec() - st))}\n')

        out(outpth, f'Summary "{m.chid}" runtime: {dt(seconds=int(sec() - start))}\n')


class Queue:
    def __init__(self, cfg: Config, scenarios_data_frame=None, mode=None):
        if scenarios_data_frame:
            self.df = scenarios_data_frame
        else:
            self.df = read_csv(join(cfg.results_path, f'{cfg.title}_set.csv'))
        self.mode = mode
        self.config = cfg
        self.queue = []
        self.status = -2
        self.results_df = self._create_or_load_df('results', ['ID', 'temp_max', 'time_crit'])
        self.errors = self._create_or_load_df('errors', ['ID', 'error_type'])
        self.to_compare = {}

    def _create_or_load_df(self, content, columns):
        try:
            df = read_csv(join(self.config.results_path, f'{self.config.title}_{content}.csv'))
            try:
                del df['Unnamed: 0']  # error in importing via rcsv
            except:
                pass
        except (FileNotFoundError, EmptyDataError):
            df = DataFrame(columns=columns)
        return df

    @property
    def n(self):
        return len(self.queue)

    # append DataFrame to CSV file
    def _writedf2csv(self, scenarios_2_save: DataFrame, ext='results'):
        if scenarios_2_save.empty:
            return False
        path = join(self.config.results_path, f'{self.config.title}_{ext}.csv')
        try:
            with open(path):
                header = False
            to_be_written = scenarios_2_save
        except FileNotFoundError:
            header = True
            to_be_written = scenarios_2_save

        to_be_written.to_csv(path_or_buf=path, mode='a', header=header, index=False)
        return True

    def _save(self, current_i_index: int, last_saved: int):
        if current_i_index == last_saved:
            return False
        self._writedf2csv(self.results_df[int((last_saved - current_i_index) / 2):])
        # consider summary a method
        self._set_status(summary(self.results_df, temp_crit(self.config.miu), self.config.RSET,
                                 savepath=self.config.results_path), len(self.results_df.index))
        return current_i_index

    def _set_status(self, rmses, iter_number, rmse_limit=0.001):
        if all([e < rmse_limit for e in rmses]):
            self.status = 2  # to be finished due to RMSE limit
        else:
            self.status = -1  # to be continued

    # find the results of the fire scenario
    def _compare(self):
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
        fire_res = [[scen.split('_')[0], *results(temp_tab)] for scen, temp_tab in self.to_compare.items()]

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

        self.to_compare.clear()

        return the_worst

    def load_set(self):
        for calc in self.df.iterrows():
            if self.mode == 'full':
                self.queue.append(Full(calc[1], self.config))
            elif not self.mode:
                self.queue.append(ThermalOnly(calc[1], self.config))

    def run(self):
        saves_no = 20
        # save every second iteration = save every scenario (when saves_no = self.queue)
        save_interval = max([int(self.n / saves_no)*2, 2])
        last_saved = 0

        # (I) iterate over tasks in queue
        for n, i in enumerate(self.queue):
            # (0a) check if iteration is a part of an erroneous scenario
            if any(self.errors.ID.isin([i.fire_id])):
                out(outpth, f'[ERROR] Scenario {i.fire_id} is affected by'
                            f' {self.errors.loc[self.errors["ID"] == i.fire_id].error_type[0]}. Passing {i.chid}.')
                last_saved = n
                continue

            # (0b) check if scenario is in results
            if any(self.results_df.ID.isin([i.fire_id])):
                out(outpth, '    Scenario {} has been already calculated. Continuing...'.format(i.fire_id))
                last_saved = n
                continue

            out(outpth, '\n {} calculations started...'.format(i.chid))

            # (1) run iteration calculations
            i_stat, i_err_mess = i.calculate()
            # (2) in case of SAFIR error utils.run_safir saves .err file and returns nonzero
            #    we don't want both iterations of erroneous scenario to be taken as results
            if i_stat != 0:
                self.to_compare.clear()
                err = [i.fire_id, i_err_mess]
                self._writedf2csv(DataFrame([err], columns=['ID', 'error_type']), ext='errors')
                self.results_df = concat([self.results_df, DataFrame([[i.fire_id]], columns=['ID'])], ignore_index=True)
                out(outpth, f'[ERROR] Scenario {i.chid} finished with FatalSafirError: {i_err_mess}')
                continue
            # (3) average temperature functions of time from subsequent iterations within scenario are saved
            self.to_compare[i.chid] = get_temp_table(i.dir_path)
            # (4a) scenario that quicker reached threshold (critical temp for steel) is saved as scenario result
            # (4b) results of fire scenario are saved to the Data Frame
            if len(self.to_compare) == 2 and not self.mode:
                self.results_df.loc[len(self.results_df)] = self._compare()

            # if not self.mode:  # [WK] not sure why this is needed here
            out(outpth, f'[OK] {i.chid} calculations finished')

            # (5) save results (CSV,txt and distributions) at every save_interval
            if (n + 1) % save_interval == 0 and n != 0:
                last_saved = self._save(n, last_saved)

            # (6) check for end condition
            if self.status >= 0:
                return self.status

        # (II) when queue gets emptied save remaining scenarios results
        self._save(self.n - 1, last_saved)
        return 1


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
    out(outpth, '========================================================================================'
                '\nmulti.py  Copyright (C) 2022  Kowalski W.'
                '\nThis program comes with ABSOLUTELY NO WARRANTY.'
                '\nThis is free software, and you are welcome to redistribute it under certain conditions.'
                '\nSee GPLv3.0 for details (https://www.gnu.org/licenses/gpl-3.0.html).'
                '\n========================================================================================\n')

    cfg = Config(argv[1])

    q = Queue(cfg)
    q.load_set()
    endstat = q.run()
    out(outpth, f'[INFO] Finished multi.py with status {endstat}')
    out(outpth, '========================================================================================'
                '\nThank you for using mcsteel package :)'
                '\nVisit project GitHub site: https://github.com/kowalskiw/mcsteel and contribute!'
                '\n========================================================================================\n')
