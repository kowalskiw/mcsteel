import copy
import socket
from os import scandir, listdir, chdir
import numpy as np
from pandas import read_csv, DataFrame
import sys
import export
from pandas.errors import EmptyDataError

from fdsafir import run_safir, user_config, Logger


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
                    print('    {} vs {}: Temperature in column scenario is higher -> added to results'.format(
                        *[i[0] for i in compared]))
                    comp_id = compared.pop(0)[0]

                elif compared[0][1] > compared[1][1]:   # bigger temp_max
                    print('   {} vs {}: Temperature in column scenario is higher -> added to results'.format(
                        *[i[0] for i in compared]))
                    comp_id = compared.pop(0)[0]

                else:
                    print('    {} vs {}: Temperature in beam scenario is higher -> added to results'.format(
                        *[i[0] for i in compared]))
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

        print('Multisimulating started...')

        for index, row in self.set.iterrows():
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

        return [self.save_res(results, export.temp_crit(self.user['miu'])) if results else -1]


'''Run multisimulation on cluster'''


class Cluster:
    def __init__(self, user, nodes_config):
        self.user = user
        self.nodes = nodes_config    # {ip: simulation directory}
        self.s = socket.socket()
        self.host = socket.gethostname()  # Get local machine name
        self.port = 12345  # Reserve a port for your service.
        self.s.bind((self.host, self.port))  # Bind to the port

    # [H] send USER file to the node
    def set_user(self, node):

        pass

    # [H] run mc.py on host
    def run_mc(self):
        pass

    # [H] divide config (chid_set.csv) to the patches
    def divide(self):
        pass

    # >>>loop>>>
    # [H] send patch of simulations config (XXX_set.csv) to the node
    def send_patch(self, node):
        # create the client socket
        print(f"[+] Connecting to {self.host}:{self.port}")
        self.s.send((self.host, self.port))
        print("[+] Connected.")

        # send the filename and filesize
        s.send(f"{filename}{SEPARATOR}{filesize}".encode())

        # start sending the file
        progress = tqdm.tqdm(range(filesize), f"Sending {filename}", unit="B", unit_scale=True, unit_divisor=1024)
        with open(filename, "rb") as f:
            while True:
                # read the bytes from the file
                bytes_read = f.read(BUFFER_SIZE)
                if not bytes_read:
                    # file transmitting is done
                    break
                # we use sendall to assure transimission in
                # busy networks
                s.sendall(bytes_read)
                # update the progress bar
                progress.update(len(bytes_read))

        # close the socket
        s.close()

    # [H] run multi.py on node
    def run_multi(self, node):
        pass

    # [H] download the results (XXX_results.csv) append it to the main DB
    def pull_results(self, ip_address):
        pass

    # [H] check if any further patches are available >BREAK< or >CONTINUE<
    # to calculate or not to calculate
    def hamletize(self):
        pass
    #       - summarize already gathered results (export.summary)
    #       - calculate RMSE and check the finish criteria
    # <<<loop<<<

    def run(self):
        def loop():
            for ip, sim_pth in self.nodes:
                self.set_user(ip)
                self.send_patch(ip+sim_pth)
                self.run_multi(ip)

        self.run_mc()
        self.divide()

        loop()

        while True:
            self.s.listen(5)
            c, addr = self.s.accept()
            self.pull_results(addr)
            if self.hamletize():
                loop()
            c.close()


class Node:
    def __init__(self):
        self.s = socket.socket()
        self.host = socket.gethostname()  # Get local machine name
        self.port = 12345  # Reserve a port for your service

    # [N] call the host when finished
    def finished(self):
        self.s.connect((self.host, self.port))
        self.s.send(b'finished')
        self.s.close()


if __name__ == '__main__':
    n = Node()
    sys.stdout = Logger('multi.log')
    Queue(sys.argv[1]).run()
    n.finished()
