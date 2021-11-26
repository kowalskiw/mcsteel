import subprocess
from time import time as sec
from datetime import timedelta as dt
from os import getcwd, chdir
from os.path import abspath as ap
from os.path import isfile, basename, dirname
from shutil import copy2
from sys import stdout, argv
import argparse as ar


# write stdout to file
class Logger(object):
    def __init__(self, filename: str):
        self.terminal = stdout
        self.log = open('{}_{}'.format(int(sec()), filename), 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass


def out(file, line):
    with open(file, 'a') as f:
        f.write(line + '\n')
    return line


# print progress of process
def progress_bar(title, current, total, bar_length=20):
    percent = float(current) * 100 / total
    arrow = '-' * int(percent / 100 * bar_length - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    print('%s: [%s%s] %d %%' % (title, arrow, spaces, percent), end='\r')


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


# find profile - first element relation in SAFIR S3D file (i.e. HEA180.tem-b_0001.tem)
# creates beam dict {'profile_no':[section_chid, section_line_no, 'bXXXXX.tem'], ...}
# creates shell dict {'profile_no':[section_chid, section_line_no, 'sXXXXX.tem'], ...}
def read_mech_input(path_to_frame):
    with open(path_to_frame) as file:
        frame_lines = file.readlines()

    t_end = 0

    tems = {}
    c = 1
    tshs = {}
    d = 1
    beam_read = False
    shell_read = False

    # read input file line by line
    for no in range(len(frame_lines)):
        lin = frame_lines[no]  # line of frame.in file

        if '.tem' in lin:  # beam TEM files
            tems[str(c)] = [''.join(lin[:-1].split('.tem')), no]
            c += 1

        elif '.tsh' in lin:  # shell TSH files
            tshs[str(d)] = [''.join(lin[:-1].split('.tsh')), no]
            d += 1

        elif 'NODOFBEAM' in lin:
            beam_read = True
        elif 'NODOFSHELL' in lin:
            beam_read = False
            shell_read = True

        elif '      ELEM' in lin:  # add first element name to each profile
            element = lin.split()
            if len(tems[element[-1]]) < 3:
                number = element[1]
                for i in range(5 - len(number)):
                    number = '0{}'.format(number)
                if beam_read:
                    tems[element[-1]].append('b{}_1.tem'.format(number))
                elif shell_read:
                    tshs[element[-1]].append('s{}_1.tem'.format(number))

        elif 'TIME' in lin.split():
            t_end = float(frame_lines[no + 1].split()[1])

    return tems, tshs, t_end


# running SAFIR simulation
def run_safir(in_file_path, safir_exe_path='C:\SAFIR\safir.exe', print_time=True):
    backpath = getcwd()
    chdir(dirname(in_file_path))
    chid = basename(in_file_path)[:-3]

    process = subprocess.Popen(' '.join([safir_exe_path, chid]), shell=False, stdout=subprocess.PIPE)
    print_all = False
    while True:
        if process.poll() is not None:
            break
        try:
            output = process.stdout.readline().strip().decode()
        except UnicodeError:
            continue
        if output:
            if print_all:
                print('    ', output)
            # check for errors
            elif 'ERROR' in output or 'forrtl' in output:
                print('[ERROR] FatalSafirError: ')
                print_all = True
                print('    ', output)
            # check for timestep
            elif 'time' in output and print_time:
                print('%s %s' % ('SAFIR started "{}" calculations: '.format(chid), output), end='\r')

    rc = process.poll()
    chdir(backpath)

    if not rc:
        print('[OK] SAFIR finished "{}" calculations at'.format(chid))


'''New, simpler and more object-oriented code'''


# return the input file path regarding to GiD catalogues or files with no further directory tree
# (make this comment better, please)
def find_paths(config_path, chid, shell=False):
    gid_paths = ['{0}\{1}.gid\{1}{2}'.format(config_path, chid, i) for i in ['.in', '-1.T0R']]
    other_in_path = '{0}\{1}.in'.format(config_path, chid)
    other_tor_paths = ['{0}\{1}{2}'.format(config_path, chid, i) for i in ['-1.T0R', 'T0R', 't0r', 'tor', 'TOR']]
    expected_paths_no = 2

    if shell:
        gid_paths = gid_paths[0]
        other_tor_paths = []
        expected_paths_no = 1

    real_paths = []

    def check(file_path):
        if isfile(file_path):
            real_paths.append(file_path)

    # check if there are GiD directories with thermals
    [check(i) for i in gid_paths]

    if len(real_paths) == expected_paths_no:
        return real_paths

    # check if files are present directly in the config directory
    else:
        check(other_in_path)
        [check(i) for i in other_tor_paths]

    if len(real_paths) == expected_paths_no:
        return real_paths
    else:
        raise FileNotFoundError('[ERROR] It is not possible to locate your {} thermal results. '
                                'Check config path {}.'.format(chid, config_path))


class ThermalTEM:
    def __init__(self, beam_type, from_mech, path_to_config, fire_model, sim_time, sim_dir):
        self.chid = from_mech[0]
        self.config_paths = find_paths(path_to_config, self.chid)  # [input file path, torsion results file path]
        self.model = fire_model.lower()
        self.beam_type = int(beam_type)
        self.t_end = sim_time
        if self.model.lower() in {'iso', 'f20', 'cold'}:
            self.first = self.chid + '.tem'
        else:
            self.first = from_mech[2]
        self.sim_dir = sim_dir
        self.line_no = from_mech[1]

    # changing input file form iso curve to natural fire mode
    def change_in(self, mech_chid):
        # open thermal analysis input file
        with open('{}\{}.in'.format(self.sim_dir, self.chid)) as file:
            init = file.readlines()

        # save backup of input file
        with open('{}\{}.bak'.format(self.sim_dir, self.chid), 'w') as file:
            file.writelines(init)

        # make changes
        for no in range(len(init)):
            line = init[no]
            # type of calculation
            if line == 'MAKE.TEM\n':
                if self.model in {'cfd', 'fds'}:
                    init[no] = 'MAKE.TEMCD\n'
                elif self.model in {'lcf', 'locafi'}:
                    init[no] = 'MAKE.TEMLF\n'
                elif self.model in {'hsm', 'hasemi'}:
                    init[no] = 'MAKE.TEMHA\n'

                # insert beam type
                [init.insert(no + 1, i) for i in ['BEAM_TYPE {}\n'.format(self.beam_type), '{}.in\n'.format(mech_chid)]]
            # change thermal load
            elif line.startswith('   F  ') and 'FISO' in line:  # choose heating boundaries with FISO or FISO0 frontier
                # change FISO0 to FISO
                if 'FISO0' in line:
                    line = 'FISO'.join(line.split('FISO0'))

                if self.model in {'cfd', 'fds'}:
                    if 'F20' not in line:
                        init[no] = 'FLUX {}'.format('CFD'.join(line[4:].split('FISO')))
                    else:
                        init[no] = 'FLUX {}'.format('NO'.join(('CFD'.join(line[4:].split('FISO'))).split('F20')))
                        init.insert(no + 1, 'NO'.join(line.split('FISO')))

                elif self.model in {'lcf', 'locafi'}:
                    if 'F20' not in line:
                        init[no] = 'FLUX {}'.format('LOCAFI'.join(line[4:].split('FISO')))
                    else:
                        init[no] = 'FLUX {}'.format('NO'.join(('LOCAFI'.join(line[4:].split('FISO'))).split('F20')))
                        init.insert(no + 1, 'NO'.join(line.split('FISO')))

                elif self.model in {'hsm', 'hasemi'}:
                    if 'F20' not in line:
                        init[no] = 'FLUX {}'.format('HASEMI'.join(line[4:].split('FISO')))
                    else:
                        init[no] = 'FLUX {}'.format('NO'.join(('HASEMI'.join(line[4:].split('FISO'))).split('F20')))
                        init.insert(no + 1, 'NO'.join(line.split('FISO')))
                elif self.model == 'F20':
                    init[no] = 'F20'.join(line.split('FISO'))

            # change convective heat transfer coefficient of steel to 35 in locafi mode according to EN1991-1-2
            elif self.model in {'lcf', 'locafi', 'hsm', 'hasemi'} and 'STEEL' in line:
                init[no + 1] = '{}'.format('35'.join(init[no + 1].split('25')))

            # change T_END
            elif ('TIME' in line) and ('END' not in line):
                try:
                    init[no + 1] = '    '.join([init[no + 1].split()[0], str(self.t_end), '\n'])
                except IndexError:
                    pass

        # write changed file
        with open('{}\{}.in'.format(self.sim_dir, self.chid), 'w') as file:
            file.writelines(init)

    # insert torsion results to the first TEM file
    def insert_tor(self, tor_lines=False):
        # open torsion analysis results if those are not provided as an argument
        if not tor_lines:
            with open(self.config_paths[1]) as file:
                tor = file.readlines()
        else:
            tor = tor_lines

        # check if torsion results already are in TEM file
        try:
            file_path = '{}\{}'.format(self.sim_dir, self.first)

            with open(file_path) as file:
                tem = file.readlines()
            if '         w\n' in tem:
                print('[OK] Torsion results are already in the TEM')
                return 0

            # looking for start of torsion results regexp in TEM file
            try:
                tor_index = tor.index('         w\n')
            except ValueError:
                raise ValueError('[ERROR] Torsion results not found in the {} file'.format(self.config_paths[1]))

            # find TEM line where torsion results should be passed
            annotation = ''
            if self.model in {'iso', 'f20'}:
                annotation = '       HOT\n'
            elif self.model == 'cfd':
                annotation = '       CFD\n'
            elif self.model in {'lcf', 'locafi'}:
                annotation = '    LOCAFI\n'
            elif self.model in {'hsm', 'hasemi'}:
                annotation = '    HASEMI\n'

            tem_with_tor = []
            try:
                tem_index = tem.index(annotation)
                tem_with_tor = tem[:tem_index] + tor[tor_index:-1] + tem[tem_index:]
            except ValueError:
                print('[WARNING] Flux constraint annotation ("HOT", "CFD", "HASEMI" or "LOCAFI") not found in'
                      '{} file'.format(self.first))
                for line in tem:
                    if 'TIME' in line:
                        tem_index = tem.index(line) - 1
                        tem_with_tor = tem[:tem_index] + tor[tor_index:-1] + [annotation] + tem[tem_index:]
                        break

            # pasting torsion results
            with open(file_path, 'w') as file:
                file.writelines(tem_with_tor)
            print('[OK] Torsion results copied to the TEM')
            return 0

        except FileNotFoundError:
            raise FileNotFoundError('[ERROR] There is no proper TEM file ({}) in {}'.format(self.first, self.sim_dir))

        except UnboundLocalError:
            print('[WARNING] The {} profile is not found in the Structural 3D .IN file'.format(self.chid))
            return -1

    def in2sim_dir(self):
        copy2(self.config_paths[0], self.sim_dir)

    # default calculations (preparations should have already been done)
    def run(self, safir_exe):
        run_safir('{}\{}.in'.format(self.sim_dir, self.chid), safir_exe_path=safir_exe)
        self.insert_tor()


class ThermalTSH:
    def __init__(self, shell_type, from_mech, path_to_config, fire_model, sim_time, sim_dir):
        self.chid = from_mech[0]
        self.config_path = find_paths(path_to_config, self.chid)[0]  # [input file path]
        self.model = fire_model.lower()
        self.shell_type = int(shell_type)
        self.t_end = sim_time
        if self.model.lower() in {'iso', 'f20', 'cold'}:
            self.first = self.chid + '.tsh'
        else:
            self.first = from_mech[2]
        self.sim_dir = sim_dir
        self.line_no = from_mech[1]

    # change input file to natural fire calculations
    def change_in(self, mech_chid):
        # open thermal analysis input file
        with open('{}\{}.in'.format(self.sim_dir, self.chid)) as file:
            init = file.readlines()

        # save backup of input file
        with open('{}.bak'.format(self.sim_dir, self.chid), 'w') as file:
            file.writelines(init)

        # make changes
        for no in range(len(init)):
            line = init[no]
            # type of calculation
            if line == 'MAKE.TEM\n':
                if self.model in {'cfd', 'fds'}:
                    init[no] = 'MAKE.TSHCD\n'
                elif self.model in {'lcf', 'locafi'}:
                    raise ValueError('[ERROR] LOCAFI model is not allowed to be used for SHELL elements.')
                elif self.model in {'hsm', 'hasemi'}:
                    init[no] = 'MAKE.TSHHA\n'

                # insert shell type and mechanical file reference
                [init.insert(no + 1, i) for i in
                 ['{}.in\n'.format(mech_chid), 'SHELL_TYPE {}\n'.format(self.shell_type)]]

            # change thermal load
            elif line.startswith('   F  ') and 'FISO' in line:  # choose heating boundaries with FISO or FISO0 frontier
                # change FISO0 to FISO
                if 'FISO0' in line:
                    line = 'FISO'.join(line.split('FISO0'))

                if self.model in {'cfd', 'fds'}:
                    if 'F20' not in line:
                        init[no] = 'FLUX {}'.format('CFD'.join(line[4:].split('FISO')))
                    else:
                        init[no] = 'FLUX {}'.format('NO'.join(('CFD'.join(line[4:].split('FISO'))).split('F20')))
                        init.insert(no + 1, 'NO'.join(line.split('FISO')))

                elif self.model in {'hsm', 'hasemi'}:
                    if 'F20' not in line:
                        init[no] = 'FLUX {}'.format('HASEMI'.join(line[4:].split('FISO')))
                    else:
                        init[no] = 'FLUX {}'.format('NO'.join(('HASEMI'.join(line[4:].split('FISO'))).split('F20')))
                        init.insert(no + 1, 'NO'.join(line.split('FISO')))
                elif self.model == 'F20':
                    init[no] = 'F20'.join(line.split('FISO'))

            # change convective heat transfer coefficient of steel to 35 in locafi mode according to EN1991-1-2
            elif self.model in {'hsm', 'hasemi'} and 'STEEL' in line:
                init[no + 1] = '{}'.format('35'.join(init[no + 1].split('25')))

            # change T_END
            elif ('TIME' in line) and ('END' not in line):
                try:
                    init[no + 1] = '    '.join([init[no + 1].split()[0], str(self.t_end), '\n'])
                except IndexError:
                    pass

        # write changed file
        with open('{}\{}.in'.format(self.sim_dir, self.chid), 'w') as file:
            file.writelines(init)

    # insert data that were lost in thermal analysis
    # not ready yet - to be developed in the future
    def insert_data(self, data_lines=False):
        print('[WARNING] There is not possible to take rebars into consideration now. These feature will be developed')
        # paste them into tsh file (begining)
        if not data_lines:
            # calculate thickness
            with open(self.first) as file:
                shell_result = file.readlines()

            read = False
            nodes = []
            for line in shell_result:
                if 'NUMBER OF POSITIONS' in line:
                    read = True
                elif 'TIME' in line:
                    break
                if read:
                    nodes += line.split()

            thickness = abs(float(nodes[0]) - float(nodes[-1]))

            #     # find rebars number
            #     with open(self.config_path) as file
            #         mech_input = file.readlines()
            #
            #     for l in shell_input:

            material = 1  # some defaults
            rebars = 0  # no rebars

            with open(self.first, 'w') as file:
                file.writelines([' THICKNESS {}\n MATERIAL {}\n REBARS {}\n\n'.format(thickness, material, rebars)] +
                                shell_result)

    def in2sim_dir(self):
        copy2(self.config_path, self.sim_dir)

    # default calculations (preparations should have already been done)
    def run(self, safir_exe):
        run_safir('{}\{}.in'.format(self.sim_dir, self.chid), safir_exe_path=safir_exe)


class Mechanical:
    def __init__(self, mech_in_path, fire_model='locafi'):
        self.thermals = []  # Thermal2D objects to be included in caculations
        self.model = fire_model
        self.input_file = mech_in_path  # path to Structural 3D input file
        self.chid = basename(mech_in_path)[:-3]
        self.sim_dir = dirname(mech_in_path)
        self.t_end = 0

    def make_thermals(self, path_to_config):
        tems, tshs, t_end = read_mech_input(self.input_file)
        self.t_end = t_end
        for k, v in tems.items():
            self.thermals.append(ThermalTEM(k, v, path_to_config, self.model, t_end, self.sim_dir))

        for k, v, in tshs.items():
            self.thermals.append(ThermalTSH(k, v, path_to_config, self.model, t_end, self.sim_dir))

        # copy input files of thermal analyses to simulation directory
        for t in self.thermals:
            t.in2sim_dir()

    # changing input file form iso curve to natural fire or cold mode
    def change_in(self):
        # open input file
        with open(self.input_file) as file:
            init = file.readlines()

        # save backup of input file
        with open('{}.bak'.format(self.input_file.rsplit('.')[-2]), 'w') as file:
            file.writelines(init)

        if self.model == 'COLD':
            # change to STATIC COLD mode
            for line in init.copy():
                if 'NCORES' in line:
                    init[init.index(line) + 1] = 'STATICCOLD' + init[init.index(line) + 1].split[1:]
                if ('M_BEAM' in line) or ('M_NODE' in line):
                    init.remove(line)
        else:
            # change TEM and TSH names in IN file (i.e. 'hea180.tem' -> 'b00001_1.tem')
            for t in self.thermals:
                init[t.line_no] = '{}\n'.format(t.first)

        # save changes
        with open('{}'.format(self.input_file), 'w') as file:
            file.writelines(init)

    def run(self, safir_exe):
        run_safir(self.input_file, safir_exe_path=safir_exe)


# to be rewritten in the future
class CheckConfig:
    def __init__(self, argums):
        self.a = argums

    # def t0r_vs_in(self, chid):
    #     tor = False
    #     # find profile 'chid' in config files and open T0R file .gid catalogue if possible
    #     for p in ['{0}\{1}.gid\{1}-1.T0R'.format(self.paths[0], chid),
    #               '{0}\{1}.T0R'.format(self.paths[0], chid),
    #               '{0}\{1}-1.T0R'.format(self.paths[0], chid),
    #               '{0}\{1}.t0r'.format(self.paths[0], chid),
    #               '{0}\{1}-1.t0r'.format(self.paths[0], chid)]:
    #         try:
    #             with open(p) as file:
    #                 tor = file.readlines()
    #             break
    #         except FileNotFoundError:
    #             continue
    #
    #     if not tor:
    #         return ['[ERROR] Torsion file of "{}" profile was not found'.format(chid)]
    #
    #     # looking for start of torsion results regexp in TEM file
    #     try:
    #         tor_index = tor.index('         w\n')
    #     except ValueError:
    #         return ['[ERROR] Torsion results were not found in "{}" TOR file'.format(chid)]
    #
    #     # find the number of elements in torsion file
    #     for l in tor:
    #         if 'NFIBERBEAM' in l:
    #             n_tor = int(l.split()[-1])
    #             break
    #
    #     # find the number of elements in initial file of thermal analysis
    #     with open('{}\{}.in'.format(self.paths[1], chid)) as file:
    #         init = file.readlines()
    #     for l in init:
    #         if 'ELEMENTS' in l:
    #             n_in = int(init[init.index(l) + 1].split()[-1])
    #
    #     # ERROR if differences found
    #     if n_in != n_tor:
    #         return ['{0} profile you use does not match {0} you put in config path ({1})'.format(chid, self.paths[0])]
    #
    #     return 0
    #
    # def in_names(self, chid):
    #     info = []
    #     if any(forb in chid for forb in ['.', ' ', ]):
    #         info.append('Filename consists of forbidden characters ("." or " ")')
    #     if chid.lower() != chid:
    #         info.append('Filename need to be lowercase')
    #
    #     return info
    #
    # def nfiber(self):
    #     # check if nfiber is correctly set
    #     # repair if necessary
    #     return []

    def check_all(self):
        pass
        #
        #
        # info = []
        #
        #
        #
        # for f in scandir(self.m.]):
        #     name = f.name
        #     in_splt = name.split('.in')
        #     chid = in_splt[0]
        #
        #     # check thermal
        #     if chid != 'frame':
        #         info.extend(self.t0r_vs_in(chid))
        #     info.extend(self.in_names(chid))
        #
        # info.extend(self.nfiber())
        #
        # if len(info) > 0:
        #     print('[INFO] While checking your input files I found some mistakes:\n', '\n'.join(info))
        #     raise ValueError('[ERROR] {} simulation was improperly set up.'
        #                      'Some mistakes were pointed above'.format(self.paths[1]))


# # to be rewritten
# # when fdsafir.py called by multi.py script of McSteel or when '-s'/'--scripted' flag used
# def scripted(safir_path, config_path, results_path):
#     for case in scandir('{}\worst'.format(results_path)):
#         chdir(case.path)
#
#         CheckConfig(config_path, case.path).check()
#
#         with open('frame.in') as frame:
#             f = frame.readlines()
#             for i in range(len(f)):
#                 if 'TIME' in f[i]:
#                     t_end = f[i + 1].split()[1]
#                     break
#         # Thermal 2D analyses of profiles
#         print('Running {} thermal analysis...'.format(case.name))
#         for i in scandir():
#             f = i.name
#             if f.endswith('.in') and not f == 'frame.in':
#                 chid = f[:-3]
#                 t = ThermalScripted(f, 'LCF', frame_chid='frame', profile_pth=f, time_end=t_end)
#                 t.alias = chid
#                 t.change_in()
#                 run_safir(chid, safir_path)
#                 t.insert_tor(config_path)
#                 del t
#
#         # Structural 3D analysis of the structure
#         print('Running {} mechanical analysis...'.format(case.name))
#         m = MechanicalScripted(frame_pth='frame')
#         m.change_in()
#         run_safir('frame', safir_path)
#
#         print('[OK] {} scenario calculations finished!'.format(case.name))


# run a single simulation with natural fire model
# (have to be already prepared/calculated for FISO)
def run_user_mode(sim_no, arguments):
    start = sec()
    m = Mechanical(arguments.results[sim_no], fire_model=arguments.model)
    m.make_thermals(arguments.config)

    for t in m.thermals:
        t.change_in(m.chid)
        t.run(arguments.safir)

    m.change_in()
    m.run(arguments.safir)

    print('\nRuntime: {}\n'.format(dt(seconds=int(sec() - start))))


def get_arguments(from_argv):
    parser = ar.ArgumentParser(description='Run SAFIR localised fire analysis automatically')

    parser.add_argument('-c', '--config', help='Path to configuration directory', required=True)
    parser.add_argument('-s', '--safir', help='Path to SAFIR executable', default='/safir.exe')
    parser.add_argument('-m', '--model', help='Type of localised fire model to be used (hasemi/hsm or locafi/lcf or'
                                              'cfd/fds)', default='locafi')
    parser.add_argument('-r', '--results', nargs='+', help='Paths to mechanical analysis IN files (one scenario ->'
                                                           'one IN file)', required=True)
    argums = parser.parse_args(args=from_argv)

    # change paths to absolute
    for k in argums.__dict__:
        if k == 'model':
            continue
        try:
            argums.__dict__[k] = ap(argums.__dict__[k])
        except TypeError:
            l = []
            for p in argums.__dict__[k]:
                l.append(ap(p))
            argums.__dict__[k] = l

    return argums


if __name__ == '__main__':
    first = argv[1]
    # if first == ('-d' or '--dev'):
    #     print('=== ATTENTION! ===\nYou have entered developer mode. It is not design for regular user - be careful,'
    #           'things don\'t work here well ;-)\n==================\n')
    #     paths = [input('safir_path = '), input('config_path = '), input('results_path = ')]
    #     scripted(*paths)

    print('\n====== Welcome to fdsafir2 ======\n fdsafir.py is one of the components of McSteel package.\n\n'
          'I am trying to run your case now. I will keep you updated on the progress. \n==================\n')

    args = get_arguments(argv[1:])

    for n in range(len(args.results)):
        CheckConfig(args).check_all()
        run_user_mode(n, args)

    print('\n==================\nThank you for using our tools. We will appreciate your feedback and bug reports'
          ' on github: https://www.github.com/kowalskiw/mcsteel\n'
          'Have a good day!\n==================\n')
