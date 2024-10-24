from datetime import datetime as dt1
from os.path import isfile, basename, dirname, abspath, join
from os import chdir, getcwd
import subprocess
from shutil import copy2
from numpy import random


# running SAFIR simulation
def run_safir(in_file_path, safir_exe_path='C:\SAFIR\safir.exe', print_time=True, fix_rlx=True, verbose=False):
    start = dt1.now()
    print(f'[INFO] Calculations started at {start}') if print_time else None
    backpath = getcwd()
    dirpath = dirname(in_file_path)
    chdir(dirpath)
    chid = basename(in_file_path)[:-3]

    print(f'Reading {chid}.in file...')
    process = subprocess.Popen(' '.join([safir_exe_path, chid]), shell=False,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print_all = verbose
    success = True
    count = 0
    err_mess = ''

    # clear output
    storage = []
    n = 0
    while True:
        n+=1
        if process.poll() is not None:
            break
        try:
            output = process.stdout.readline().strip().decode()
        except UnicodeError:
            continue
        # else print SAFIR output
        if output:
            step = output[7:]
            if print_all:
                print('    ', output)
            # check for errors in stdout (handled errors from SAFIR are redirected to stdout)
            elif 'ERROR' in output or 'forrtl' in output:
                print('[ERROR] FatalSafirError: ')
                print_all = True
                success = False
                print('    ', output)
                with open(f'{chid}.err', 'w') as errfile:
                    errfile.write(output)
                err_mess = output
            # look for characters at the beginning of the simulation
            if '======================' in output:
                count += 1
            # check for timestep
            elif 'time' in output:
                if print_time:
                    print(f'SAFIR started "{chid}" (sim #{count}) calculations: {step}', end='\r')
        else:
            step = ''

    rc = process.poll()
    chdir(backpath)

    if not rc:
        if success:
            print(f'[OK] SAFIR finished {count} "{chid}" calculations at{step}')
            print(f'[INFO] Computing time: {dt1.now() - start}')
            repair_relax(f'{dirpath}\\{chid}.XML') if fix_rlx else None
            return 0, err_mess
        else:
            print(f'[WARNING] SAFIR finished "{chid}" calculations with error!')
            return -1, err_mess
    else:
        err = process.stderr.read().decode()
        # check for errors
        if err:
            print('[ERROR] FatalSafirError: ')
            print('    ', err)
            with open(f'{chid}.err', 'w') as errfile:
                errfile.write(err)
            return -1, err
        raise ChildProcessError('Something went wrong with SAFIR subprocess')


# repair invalid relaxations results in XML file
def repair_relax(path_to_xml, copyxml=True):
    rlx_lines = []
    index = 0
    fixed = 0

    with open(path_to_xml) as xmlfile:
        for line in xmlfile:
            if 'RLX' in line:
                rlx_lines.append('0'.join('-1'.join(line.split('-0.100E+01')).split('0.000E+00')))
                fixed += 1
            else:
                rlx_lines.append(line)

            index += 1

    with open(f'{path_to_xml[:-4]}_fixed.XML' if copyxml else path_to_xml, 'w') as newxml:
        newxml.writelines(rlx_lines)

    print(f'[OK] {fixed} XML file lines fixed (relaxations bug)')

    return 0


# print progress of the process
def progress_bar(title, current, total, bar_length=20):
    percent = float(current) * 100 / total
    arrow = '-' * int(percent / 100 * bar_length - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    print('%s: [%s%s] %d %%' % (title, arrow, spaces, percent), end='\r')


# append to the output file
def out(outfile, line, silent=False):
    with open(outfile, 'a') as f:
        f.write(line + '\n')
    if not silent:
        print(line)
    return line


# return user configuration dict
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


# triangular distribution sampler
def triangular(left, right, mode=False):
    if not mode:  # default mode in 1/3 of the left-right distance
        mode = (right - left) / 3 + left
    return random.triangular(left, mode, right)


class Config:
    def __init__(self, user_filepath):
        self.config = user_config(user_filepath)
        self.title = self.config['title']
        self.config_path = self.config['config_path']
        self.results_path = self.config['results_path']
        self.fire_type = self.config['fire_type']
        self.safir_path = self.optional('safir_path')
        self.miu = self.optional('miu')
        self.RSET = self.config['RSET']
        self.max_iterations = self.optional('max_iterations')
        self.stop_condition = self.config['stop_condition']
        self.time_end = self.config['time_end']
        self.occupancy = self.optional('occupancy')
        if self.fire_type.lower() != 'cfast':
            self.fuel = self.config['fuel']
        print('[OK] Configuration file loaded')

    def optional(self, key):
        try:
            if key == 'max_iterations':
                return int(self.config[key])
            return self.config[key]
        except KeyError:
           return None

    def section_path(self, section_chid):
        return find_paths(self.config_path, section_chid, shell=False)[0]


# find profile - first element relation in SAFIR S3D file (ent.e. HEA180.tem-b_0001.tem)
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
            if len(element) > 7:
                continue
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


# return the input file path regarding to GiD catalogues or files with no further directory tree
# (make this comment better, please)
def find_paths(config_path, chid, shell=False):
    in_paths = [join(config_path, f'{chid}.{i}') for i in ['in', 'IN', 'In', 'iN']]
    tor_paths = [join(config_path, f'{chid}{i}') for i in ['-t.TOR', '-t.T0R', '.TOR', '.T0R']]
    encaps_in_paths = [join(dirname(i), chid, basename(i)) for i in in_paths]
    encaps_tor_paths = [join(dirname(i), chid, basename(i)) for i in tor_paths]
    expected_paths_no = 2

    if shell:
        expected_paths_no = 1

    real_paths = []

    def check(file_path):
        if isfile(file_path):
            real_paths.append(abspath(file_path))
            return True

    for p in in_paths + encaps_in_paths:
        if check(p):
            break

    for p in tor_paths + encaps_tor_paths:
        if check(p):
            break

    if len(real_paths) == expected_paths_no:
        return real_paths
    else:
        raise FileNotFoundError('[ERROR] It is not possible to locate your {} thermal or torsion results. '
                                'Check config path {}.'.format(chid, config_path))

def find_paths2019(config_path, chid, shell=False):
    gid_paths = [join(config_path, f'{chid}.gid', f'{chid}.{i}') for i in ['.in', '-1.T0R']]
    other_in_path = [join(config_path, f'{chid}.{i}') for i in ['in', 'IN', 'In', 'iN']]
    other_tor_paths = [join(config_path, f'{chid}{i}') for i in ['-1.T0R', '.T0R', '.t0r', '.tor', '.TOR']]
    expected_paths_no = 2

    if shell:
        gid_paths = gid_paths[0]
        other_tor_paths = []
        expected_paths_no = 1

    real_paths = []

    def check(file_path):
        if isfile(file_path):
            real_paths.append(abspath(file_path))

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
        raise FileNotFoundError('[ERROR] It is not possible to locate your {} thermal or torsion results. '
                                'Check config path {}.'.format(chid, config_path))


class ThermalTEM:
    def __init__(self, beam_type, from_mech, path_to_config, fire_model, sim_time, sim_dir):
        self.chid = from_mech[0]
        self.config_paths = find_paths(path_to_config, self.chid)  # [input file path, torsion results file path]
        self.model = fire_model.lower()
        self.beam_type = int(beam_type)
        self.t_end = sim_time
        if self.model in {'iso', 'standard', 'fiso', 'f20', 'cold'}:
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
            spltd_l = line.split()
            # type of calculation
            if line == 'MAKE.TEM\n' and self.model not in {'cold', 'f20', 'iso', 'fiso', 'standard'}:
                if self.model in {'cfd', 'fds'}:
                    init[no] = 'MAKE.TEMCD\n'
                elif self.model in {'lcf', 'locafi'}:
                    init[no] = 'MAKE.TEMLF\n'
                elif self.model in {'hsm', 'hasemi'}:
                    init[no] = 'MAKE.TEMHA\n'

                # insert beam type
                [init.insert(no + 1, i) for i in ['BEAM_TYPE {}\n'.format(self.beam_type), '{}.in\n'.format(mech_chid)]]

            # change thermal attack functions
            elif ('F' in spltd_l) and ('FISO' in spltd_l):  # choose heating boundaries with FISO or FISO0 frontier
                # change FISO0 to FISO
                if 'FISO0' in spltd_l:
                    line = 'FISO'.join(line.split('FISO0'))

                # choose function to be changed with
                thermal_attack = 'F20'
                if self.model in {'cfd', 'fds'}:
                    thermal_attack = 'CFD'
                elif self.model in {'lcf', 'locafi'}:
                    thermal_attack = 'LOCAFI'
                elif self.model in {'hsm', 'hasemi'}:
                    thermal_attack = 'HASEMI'
                elif self.model in {'iso', 'fiso', 'standard'}:
                    break

                if thermal_attack == 'F20':
                    init[no] = 'F20'.join(line.split('FISO'))

                elif 'F20' not in line:
                    init[no] = f'FLUX {" ".join(map(lambda x: x.replace("FISO", thermal_attack), line.split()[1:]))}\n'
                else:
                    # replace FISO frontiers with thermal_attack FLUX
                    init[no] = f'FLUX {" ".join(map(lambda x: x.replace("FISO", thermal_attack), line.split()[1:]))}\n'
                    init[no] = " ".join(map(lambda x: x.replace("F20", "NO"), init[no].split()))
                    # append remaining F20 frontiers
                    init.insert(no+1, 'NO'.join(line.split('FISO')) + '\n')

            # change convective heat transfer coefficient of steel to 35 in locafi mode according to EN1991-1-2
            elif self.model in {'lcf', 'locafi', 'hsm', 'hasemi'} and 'STEEL' in line:
                init[no + 1] = '{}'.format('35'.join(init[no + 1].split('25')))

            # change T_END
            elif ('TIME' in line) and ('END' not in line):
                try:
                    timeline = init[no+1].split()
                    init[no+1] = '    '.join([timeline[0], str(self.t_end), timeline[2], '\n'])
                except IndexError:
                    pass

        # write changed file
        with open('{}\{}.in'.format(self.sim_dir, self.chid), 'w') as file:
            file.writelines(init)

    # insert torsion results to the first TEM file
    def insert_tor(self):
        tem_with_tor = []

        # check if torsion results already are in TEM file
        try:
            tem_file_path = '{}\{}'.format(self.sim_dir, self.first)
            with open(tem_file_path) as file:
                tem = file.read()

        except FileNotFoundError:
            raise FileNotFoundError('[ERROR] There is no proper TEM file ({}) in {}'.format(self.first, self.sim_dir))

        with open(self.config_paths[1]) as file:
            tor = file.read()

        # try:
        # looking for torsion results regexp in TEM file
        if all(t in tem for t in ['GJ', 'w\n']):
            print('[OK] Torsion results are already in the TEM')
            tem_with_tor = tem

        else:
            # check for torsion results in T0R file
            if all(t in tor for t in ['GJ', 'w\n']):
                tor_indexes = [tor.index(i) for i in ['w\n', 'COLD']]
            else:
                raise ValueError('[ERROR] Torsion results not found in the {} file'.format(self.config_paths[1]))

            # find TEM line where torsion results should be passed
            # annotation = ''
            # if self.model in {'cold', 'f20', 'iso', 'standard', 'fiso'}:
            #     annotation = '       HOT\n'
            # elif self.model == 'cfd':
            #     annotation = '       CFD\n'
            # elif self.model in {'lcf', 'locafi'}:
            #     annotation = '    LOCAFI\n'
            # elif self.model in {'hsm', 'hasemi'}:
            #     annotation = '    HASEMI\n'

            # insert torsion results to thermal results file
            tem_parts = []
            for i in ['HOT', 'CFD', 'HASEMI', 'LOCAFI']:
                if i in tem:
                    tem_parts = tem.split(i)
                    tem_with_tor = i.join([tem_parts[0] + tor[tor_indexes[0]:tor_indexes[1]], tem_parts[1]])
                    break

            if not tem_parts:
                raise ValueError('[ERROR] Flux constraint annotation ("HOT", "CFD", "HASEMI" or "LOCAFI") not'
                                 'found in {} file'.format(self.first))

        # change for cold if necessary
        if self.model in {'cold', 'f20'}:
            tem_with_tor = 'COLD'.join(tem_with_tor.split('HOT'))

        # pasting torsion results
        with open(tem_file_path, 'w') as file:
            file.writelines(tem_with_tor)
        print('[OK] Torsion results copied to the TEM')
        return 0

        # except UnboundLocalError:
        #     print('[WARNING] The {} profile is not found in the Structural 3D .IN file'.format(self.chid))
        #     return -1

    def in2sim_dir(self):
        copy2(self.config_paths[0], self.sim_dir)

    # default calculations (preparations should have already been done)
    def run(self, safir_exe):
        run_safir('{}\{}.in'.format(self.sim_dir, self.chid), safir_exe_path=safir_exe, fix_rlx=False)
        self.insert_tor()


class ThermalTSH:
    def __init__(self, shell_type, from_mech, path_to_config, fire_model, sim_time, sim_dir):
        self.chid = from_mech[0]
        self.config_path = find_paths(path_to_config, self.chid)[0]  # [input file path]
        self.model = fire_model.lower()
        self.shell_type = int(shell_type)
        self.t_end = sim_time
        if self.model.lower() in {'iso', 'fiso', 'standard', 'f20', 'cold'}:
            self.first = self.chid + '.tsh'
        else:
            self.first = from_mech[2]
        self.sim_dir = sim_dir
        self.line_no = from_mech[1]

    # change input file to natural fire calculations
    def change_in(self, mech_chid):
        in_file_path = '{}\{}.in'.format(self.sim_dir, self.chid)

        # open thermal analysis input file
        with open(in_file_path) as file:
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

            # change thermal attack functions
            elif line.startswith('   F  ') and 'FISO' in line:  # choose heating boundaries with FISO or FISO0 frontier

                # choose function to be changed with
                thermal_attack = 'F20'
                if self.model in {'cfd', 'fds'}:
                    thermal_attack = 'CFD'
                elif self.model in {'lcf', 'locafi'}:
                    raise ValueError('[ERROR] LOCAFI model is not allowed to be used for SHELL elements.')
                elif self.model in {'hsm', 'hasemi'}:
                    thermal_attack = 'HASEMI'
                elif self.model in {'iso', 'fiso', 'standard'}:
                    break

                # replace FISO0 with FISO
                if 'FISO0' in line:
                    line = 'FISO'.join(line.split('FISO0'))

                # change thermal attack functions
                if thermal_attack == 'F20':
                    init[no] = 'F20'.join(line.split('FISO'))
                elif 'F20' not in line:
                    init[no] = 'FLUX {}'.format(thermal_attack.join(line[4:].split('FISO')))
                else:
                    init[no] = 'FLUX {}'.format('NO'.join((thermal_attack.join(line[4:].split('FISO'))).split('F20')))
                    init.insert(no + 1, 'NO'.join(line.split('FISO')))

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
        with open(in_file_path, 'w') as file:
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
        run_safir('{}\{}.in'.format(self.sim_dir, self.chid), safir_exe_path=safir_exe, fix_rlx=False)


class Mechanical:
    def __init__(self, mech_in_path, fire_model='locafi'):
        self.thermals = []  # Thermal2D objects to be included in caculations
        self.model = fire_model.lower()
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
        with open('{}\{}.bak'.format(self.sim_dir, self.chid), 'w') as file:
            file.writelines(init)

        # change to STATIC COLD mode
        if self.model in {'cold', 'f20'}:
            for line in init.copy():
                if 'NCORES' in line:
                    init[init.index(line) + 1] = 'STATICCOLD  {}\n'.format(init[init.index(line) + 1].split()[-1])
                if any(m in line for m in ['M_BEAM', 'M_NODE', 'MASS', 'END_MASS']):
                    init.remove(line)

        # change TEM and TSH names in IN file (ent.e. 'hea180.tem' -> 'b00001_1.tem')
        else:
            for t in self.thermals:
                init[t.line_no] = '{}\n'.format(t.first)

        # save changes
        with open('{}'.format(self.input_file), 'w') as file:
            file.writelines(init)

    def run(self, safir_exe):
        run_safir(self.input_file, safir_exe_path=safir_exe)


# to be rewritten in the future
class Check:
    def __init__(self, mechanical: Mechanical):
        self.mech = mechanical

    def t0r_vs_in(self, thermal: ThermalTEM):
        # check if torsion results are compatible with thermal input file
        info = []

        try:
            with open(thermal.config_paths[1]) as file:
                tor = file.read()
        except FileNotFoundError:
            return ['Torsion file of "{}" profile was not found'.format(thermal.chid)]

        # looking for start of torsion results regexp in T0R file
        try:
            [tor.index(i) for i in ['w\n', 'GJ']]
        except ValueError:
            info.append(['Torsion results were not found in "{}" TOR file'.format(thermal.chid)])

        # find the number of elements in torsion file
        n_tor = int(tor[:50].split('NFIBERBEAM')[1].split()[0])

        # find the number of elements in initial file of thermal analysis
        try:
            with open(thermal.config_paths[0]) as file:
                n_in = int(file.read()[:500].split('SOLID')[1].split()[0])
        except FileNotFoundError:
            return ['Thermal analysis input file of "{}" profile was not found'.format(thermal.chid)]

        # ERROR if differences found
        if n_in != n_tor:
            info.append(['Numbers of fibers in torsion ({}) and thermal ({}) analyses do not match.'.format(
                n_tor, n_in)])

        return info

    def name(self, file_name):
        info = []

        if any(forb in file_name for forb in ['.', ' ', ]):
            info.append('There is forbidden character in filename "{}"'.format(file_name))

        if file_name.lower() != file_name:
            info.append('Filename "{}" has to be lowercase'.format(file_name))

        return info

    def nfiber(self):
        # check if nfiber is correctly set
        with open(self.mech.input_file) as file:
            in_mech_file = file.readlines()
            nfiber = int(''.join(in_mech_file)[:500].split('NFIBER')[1].split()[0])

        nsolid_max = nfiber
        for t in self.mech.thermals:
            with open(t.config_paths[0]) as file:
                n_in = int(file.read()[:500].split('SOLID')[1].split()[0])
            if n_in > nsolid_max:
                nsolid_max = n_in

        # write with updated nfiber if necessary
        if nsolid_max > nfiber:
            print('[WARNING] NFIBER is too low - is {} and should be {}. I will fix it.'.format(nfiber, nsolid_max))

            for line in in_mech_file:
                if 'NFIBER' in line:
                    in_mech_file[in_mech_file.index(line)] = '    NFIBER    {}\n'.format(nsolid_max)
                    break

            with open(self.mech.input_file, 'w') as file:
                file.writelines(in_mech_file)

        return []

    def full_mech(self):
        info = []

        # check mechanical
        [info.extend(check) for check in [self.name(self.mech.chid), self.nfiber()]]

        # check thermals
        for t in self.mech.thermals:
            [info.extend(check) for check in [self.name(t.chid), self.t0r_vs_in(t)]]

        if len(info) > 0:
            print('[INFO] While checking your input files I found some mistakes:\n', '\n'.join(info))
            raise ValueError('[ERROR] {} simulation was improperly set up.'.format(self.mech.chid))
        else:
            print('[OK] Config files seems to be OK')
