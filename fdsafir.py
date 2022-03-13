import subprocess
from time import time as sec
from os import getcwd, listdir, chdir, scandir
from shutil import copyfile
from sys import stdout, argv


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
def progressBar(title, current, total, bar_length=20):
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
# creates beam dict {'profile_no':[section_name, section_line_no, 'bXXXXX.tem'], ...}
# creates shell dict {'profile_no':[section_name, section_line_no, 'sXXXXX.tem'], ...}

def sections(frame):
    with open('{}.in'.format(frame)) as file:
        frame_lines = file.readlines()

    tems = {}
    c = 1
    tshs = {}
    d = 1
    beam_read = False
    shell_read = False

    for n in range(len(frame_lines)):  # find sections files
        l = frame_lines[n]  # line of frame.in file

        if '.tem' in l:  # beam TEM files
            tems[str(c)] = [l[:-1], n]
            c += 1

        elif '.tsh' in l:  # shell TSH files
            tshs[str(d)] = [l[:-1], n]
            d += 1

        elif 'NDOFBEAM' in l:
            beam_read = True
        elif 'NDOFSHELL' in l:
            beam_read = False
            shell_read = True

        elif '      ELEM' in l:  # add first element name to each profile
            element = l.split()
            if len(tems[element[-1]]) < 3:
                number = element[1]
                for i in range(5 - len(number)):
                    number = '0{}'.format(number)
                if beam_read:
                    tems[element[-1]].append('b{}_1.tem'.format(number))
                elif shell_read:
                    tshs[element[-1]].append('s{}_1.tem'.format(number))

    return tems, tshs


'''Safir_thermal2D analyses'''


class Thermal:
    def __init__(self, chid: str, model: str, time_end: str = 1800):
        self.chid = ''.join(chid.split('.')[:-1]) if chid.endswith('.in') or chid.endswith('.gid') else chid
        self.model = model
        self.path = getcwd()
        self.t_end = time_end

        info = '[INFO] Default mode - single fdsafir.py simulation'
        self.frame = '{}'.format('frame')
        self.profile_pth = '{0}\{1}.gid\{1}.in'.format(self.path, self.chid)
        self.alias = self.chid
        self.tems, self.tshs = sections('{0}\\{1}.gid\\{1}'.format(self.path, self.frame))
        self.path_prefix = '{}\{}\\'.format(self.path, self.chid)

    def is_in_structural(self):
        if self.chid in self.tems:
            return True

        return False

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open(self.profile_pth) as file:
            init = file.readlines()

        with open('{}.bak'.format(self.alias), 'w') as file:
            file.writelines(init)

        # make changes
        for n in range(len(init)):
            l = init[n]
            # type of calculation
            if l == 'MAKE.TEM\n':
                if self.model == 'CFD':
                    init[n] = 'MAKE.TEMCD\n'
                elif self.model in {'LCF', 'LOCAFI'}:
                    init[n] = 'MAKE.TEMLF\n'
                elif self.model in {'HSM', 'HASEMI'}:
                    init[n] = 'MAKE.TEMHA\n'

                # insert beam type
                for k, v in self.tems.items():
                    if v[0][:-4] == self.alias:
                        [init.insert(n + 1, i) for i in ['BEAM_TYPE {}\n'.format(k), '{}.in\n'.format(self.frame)]]

            elif l == 'MAKE.TSH\n':
                if self.model == 'CFD':
                    init[n] = 'MAKE.TSHCD\n'
                elif self.model in {'LCF', 'LOCAFI'}:
                    raise ValueError('[ERROR] It is not possible to calculate shell elements with LOCAFI fire model')
                elif self.model in {'HSM', 'HASEMI'}:
                    init[n] = 'MAKE.TSHHA\n'

                # insert shell type
                for k, v in self.tems.items():
                    if v[0][:-4] == self.alias:
                        [init.insert(n + 1, i) for i in ['SHELL_TYPE {}\n'.format(k), '{}.in\n'.format(self.frame)]]

            # change thermal load
            elif l.startswith('   F  ') and 'FISO' in l:  # choose heating boundaries with FISO or FISO0frontier
                # change FISO0 to FISO
                if 'FISO0' in l:
                    l = 'FISO'.join(l.split('FISO0'))

                if self.model == 'CFD':
                    if 'F20' not in l:
                        init[n] = 'FLUX {}'.format('CFD'.join(l[4:].split('FISO')))
                    else:
                        init[n] = 'FLUX {}'.format('NO'.join(('CFD'.join(l[4:].split('FISO'))).split('F20')))
                        init.insert(n + 1, 'NO'.join(l.split('FISO')))
                elif self.model in {'LCF', 'LOCAFI'}:
                    if 'F20' not in l:
                        init[n] = 'FLUX {}'.format('LOCAFI'.join(l[4:].split('FISO')))
                    else:
                        init[n] = 'FLUX {}'.format('NO'.join(('LOCAFI'.join(l[4:].split('FISO'))).split('F20')))
                        init.insert(n + 1, 'NO'.join(l.split('FISO')))
                elif self.model in {'HSM', 'HASEMI'}:
                    if 'F20' not in l:
                        init[n] = 'FLUX {}'.format('HASEMI'.join(l[4:].split('FISO')))
                    else:
                        init[n] = 'FLUX {}'.format('NO'.join(('HASEMI'.join(l[4:].split('FISO'))).split('F20')))
                        init.insert(n + 1, 'NO'.join(l.split('FISO')))
                elif self.model == 'F20':
                    init[n] = 'F20'.join(l.split('FISO'))

            # change convective heat transfer coefficient of steel to 35 in locafi mode
            elif self.model in {'LCF', 'LOCAFI', 'HSM', 'HASEMI'} and 'STEEL' in l:
                init[n + 1] = '{}'.format('35'.join(init[n + 1].split('25')))

            # change T_END
            elif ('TIME' in l) and ('END' not in l):
                try:
                    init[n + 1] = '    '.join([init[n + 1].split()[0], str(self.t_end), '\n'])
                except IndexError:
                    pass

        # write changed file
        with open(self.profile_pth, 'w') as file:
            file.writelines(init)

    # copying fire and frame file to section catalogue (natural fire only)
    def copy_ess(self):
        copyfile('{0}.gid\{0}.in'.format('frame'), '{}\{}.gid\{}.in'.format(self.path, self.chid, 'frame'))

        if self.model == 'CFD':
            copyfile('cfd.txt', '{}.gid\cfd.txt'.format(self.chid))
        elif self.model in {'LCF' or 'LOCAFI'}:
            copyfile('locafi.txt', '{}.gid\locafi.txt'.format(self.chid))
        elif self.model in {'HSM' or 'HASEMI'}:
            copyfile('hasemi.txt', '{}.gid\hasemi.txt'.format(self.chid))

    # adding torsion analysis results to first TEM file
    def insert_tor(self, config_path='.'):
        # chdir('{}\{}.gid'.format(self.path, self.chid))
        with open('{0}\{1}.gid\{1}-1.T0R'.format(config_path, self.chid)) as file:
            tor = file.readlines()

        # picking TEM file to insert torsion results to

        if self.model == ('ISO' or 'F20'):
            first_b = self.chid + '.tem'
        else:
            for v in self.tems.values():
                if v[0][:-4] == self.chid:
                    first_b = v[-1]
        # check if torsion results already are in TEM file
        try:
            file_path = '{}\{}'.format(self.path_prefix, first_b)

            with open(file_path) as file:
                tem = file.readlines()
            if '         w\n' in tem:
                # chdir('..')
                print('[OK] Torsion results are already copied to the TEM')
                return 0

            # looking for start of torsion results regexp in TEM file
            try:
                tor_index = tor.index('         w\n')
            except ValueError:
                raise ValueError('[ERROR] Torsion results not found in the TOR')

            # find TEM line where torsion results should be passed
            annotation = ''
            if self.model in {'ISO', 'F20'}:
                annotation = '       HOT\n'
            elif self.model == 'CFD':
                annotation = '       CFD\n'
            elif self.model in {'LCF', 'LOCAFI'}:
                annotation = '    LOCAFI\n'
            elif self.model in {'HSM', 'HASEMI'}:
                annotation = '    HASEMI\n'

            tem_index = int
            tem_with_tor = []
            try:
                tem_index = tem.index(annotation)
                tem_with_tor = tem[:tem_index] + tor[tor_index:-1] + tem[tem_index:]
            except ValueError:
                print('[WARNING] Flux constraint annotation ("HOT", "CFD", "HASEMI" or "LOCAFI") not found in the'
                      '{} file'.format(first_b))
                for l in tem:
                    if 'TIME' in l:
                        tem_index = tem.index(l) - 1
                        tem_with_tor = tem[:tem_index] + tor[tor_index:-1] + [annotation] + tem[tem_index:]
                        break

            # pasting torsion results
            with open(file_path, 'w') as file:
                file.writelines(tem_with_tor)
            # chdir('..')
            print('[OK] Torsion results copied to the TEM')

            return 0

        except FileNotFoundError:
            raise FileNotFoundError('[ERROR] There is no proper TEM file ({}) in {}'.format(first_b, getcwd()))

        except UnboundLocalError:
            print('[WARNING] The {} profile is not found in the Structural 3D .IN file'.format(self.chid))
            return -1

    # running single SAFIR simulation
    def run(self):
        # iso fire curve
        if self.model == 'ISO':
            run_safir(self.chid)
            self.insert_tor()
            print('[OK] ISO-curve {} thermal 2D analysis finished\n\n'.format(self.chid))

        # natural fire or F20
        else:
            if not self.model == 'F20':
                self.copy_ess()
            try:
                if not self.is_in_structural():
                    print('[WARNING] The {} profile is not found in the Structural 3D .IN file. Passing the scenario'
                          ''.format(self.chid))
                    return False
                self.change_in()
                run_safir(self.chid)
                self.insert_tor()
                print('[OK] {}-data {} Thermal 2D analysis finished\n\n'.format(self.model, self.chid))
            except ValueError:
                raise ValueError('[ERROR] change_in not possible')


class ThermalScripted(Thermal):
    def __init__(self, chid: str, model: str, frame_chid: str = 'default', profile_pth: str = 'default',
                 time_end: str = 1800):
        super().__init__(chid, model, time_end)

        info = '[INFO] scripted mode - multisimulations'
        self.frame = frame_chid
        self.profile_pth = profile_pth
        self.alias = 'safir'
        self.tems, self.tshs = sections(self.frame)
        self.path_prefix = self.path


'''Safir_structural3D analyses'''


class Mechanical:
    def __init__(self, t2Ds, model, frame_pth='frame.gid/frame'):
        self.model = model  # ISO/other
        self.t2Ds = t2Ds
        self.frame = frame_pth

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open('{}.in'.format(self.frame)) as file:
            init = file.readlines()

        with open('{}.bak'.format(self.frame.rsplit('.')[-1]), 'w') as file:
            file.writelines(init)

        if self.model == 'COLD':
            # change to STATIC COLD mode
            for l in init.copy():
                if 'NCORES' in l:
                    init[init.index(l) + 1] = 'STATICCOLD' + init[init.index(l) + 1].split[1:]
                if 'M_BEAM' in l:
                    init.remove(l)
        else:
            # change TEM names in IN file
            tems, tshs = sections(self.frame)
            for v in tems.values():
                init[v[1]] = '{}\n'.format(v[-1])
            # change TSH names in IN file
            for v in tshs.values():
                init[v[1]] = '{}\n'.format(v[-1])

        with open('{}.in'.format(self.frame), 'w') as file:
            file.writelines(init)

    # copying natural fire TEM files
    def copy_tems(self):
        for prof in self.t2Ds:
            for f in listdir('{}\{}'.format(getcwd(), prof)):
                if f.endswith('.tem'):
                    copyfile('{}\{}\{}'.format(getcwd(), prof, f), 'frame.gid\{}'.format(f))

    # running single SAFIR simulation
    def run(self):
        # iso fire curve
        if self.model in {'ISO', 'F20'}:
            self.copy_tems()
            run_safir('frame')
            print('\n[OK] {}}-curve Structural 3D analysis finished\n\n'.format(self.model))

        # natural fire
        else:
            self.change_in()
            self.copy_tems()
            run_safir('frame')
            print('\n [OK] Natural fire Structural 3D analysis finished\n\n')

#
# class CheckConfig:
#     def __init__(self, config_path: str, scenario_path: str):
#         self.paths = [config_path, scenario_path]
#
#     def t0r_vs_in(self, chid):
#         tor = False
#         # find profile 'chid' in config files and open T0R file .gid catalogue if possible
#         for p in ['{0}\{1}.gid\{1}-1.T0R'.format(self.paths[0], chid),
#                   '{0}\{1}.T0R'.format(self.paths[0], chid),
#                   '{0}\{1}-1.T0R'.format(self.paths[0], chid),
#                   '{0}\{1}.t0r'.format(self.paths[0], chid),
#                   '{0}\{1}-1.t0r'.format(self.paths[0], chid)]:
#             try:
#                 with open(p) as file:
#                     tor = file.readlines()
#                 break
#             except FileNotFoundError:
#                 continue
#
#         if not tor:
#             return ['[ERROR] Torsion file of "{}" profile was not found'.format(chid)]
#
#         # looking for start of torsion results regexp in TEM file
#         try:
#             tor_index = tor.index('         w\n')
#         except ValueError:
#             return ['[ERROR] Torsion results were not found in "{}" TOR file'.format(chid)]
#
#         # find the number of elements in torsion file
#         for l in tor:
#             if 'NFIBERBEAM' in l:
#                 n_tor = int(l.split()[-1])
#                 break
#
#         # find the number of elements in initial file of thermal analysis
#         with open('{}\{}.in'.format(self.paths[1], chid)) as file:
#             init = file.readlines()
#         for l in init:
#             if 'ELEMENTS' in l:
#                 n_in = int(init[init.index(l) + 1].split()[-1])
#
#         # ERROR if differences found
#         if n_in != n_tor:
#             return ['{0} profile you use does not match {0} you put in config path ({1})'.format(chid, self.paths[0])]
#
#         return 0
#
#     def in_names(self, chid):
#         info = []
#         if any(forb in chid for forb in ['.', ' ', ]):
#             info.append('Filename consists of forbidden characters ("." or " ")')
#         if chid.lower() != chid:
#             info.append('Filename need to be lowercase')
#
#         return info
#
#     def nfiber(self):
#         # check if nfiber is correctly set
#         # repair if necessary
#         return []
#
#     def check(self):
#         for f in scandir(self.paths[1]):
#             name = f.name
#             in_splt = name.split('.in')
#             chid = in_splt[0]
#
#             # check thermal
#             if chid != 'frame':
#                 info.extend(self.t0r_vs_in(chid))
#                 info.extend(self.nfiber())
#             info.extend(self.in_names(chid))
#
#             if len(info) > 0:
#                 raise ValueError('[ERROR] While checking your input files I found some mistakes:\n', '\n'.join(info))


# running SAFIR simulation in shell
def run_safir(in_file, safir='C:\SAFIR', gid_style=False):
    safir_path = '{}\safir.exe'.format(safir)

    if not gid_style:
        chid = in_file.split('.')[0]
    else:
        chid = in_file
        chdir('.\{}.gid'.format(chid))

    print('SAFIR started {} calculations...'.format(chid))
    output = subprocess.check_output(' '.join([safir_path, '"{}"'.format(chid)]), universal_newlines=True)

    # check for errors in output:
    try:
        err_text = output.split('ERROR')[1]
        with open('{}.err'.format(chid), 'w') as err:
            mess = '[ERROR] FatalSafirError: Simulation {} \n{}\nScenario passed, simulations saved to \'err.csv\',' \
                   ' continuing...'.format(chid, err_text)
            err.write(mess)
        print(mess)

    except IndexError:
        print('[OK] SAFIR finished {} calculations'.format(chid))

    chdir('..') if not gid_style else 0


def main(model, calc_type='s3dt2d', path='.'):
    model = model.upper()
    calc_type = calc_type.lower()

    wd = getcwd()
    if path != wd:
        chdir(path)
    folders = listdir(getcwd())
    folders.remove('frame.gid')
    for f in reversed(folders):
        if not f.endswith('.gid'):
            folders.remove(f)

    # these commented lines below are for ISO analysis, probably not necessary in common use
    # for prof in folders:
    #     Thermal(prof, 'ISO').run()

    # Mechanical(folders, 'ISO').run()       # frame structural analysis - ISO curve

    if 't2d' in calc_type:
        for prof in folders:
            Thermal(prof, model).run()  # normal temperatre mode

    if 's3d' in calc_type:
        Mechanical(folders, model).run()  # frame structural analysis

    print('\n[OK] All SAFIR calculations finished, well done engineer!\n\n')

    if path != wd:
        chdir(wd)


# when fdsafir.py called by multi.py script of McSteel or when '-s'/'--scripted' flag used
def scripted(safir_path, config_path, results_path):
    for case in scandir('{}\worst'.format(results_path)):
        chdir(case.path)

        # CheckConfig(config_path, case.path).check()

        with open('frame.in') as frame:
            f = frame.readlines()
            for i in range(len(f)):
                if 'TIME' in f[i]:
                    t_end = f[i + 1].split()[1]
                    break
        # Thermal 2D analyses of profiles
        print('Running {} thermal analysis...'.format(case.name))
        for i in scandir():
            f = i.name
            if f.endswith('.in') and not f == 'frame.in':
                chid = f[:-3]
                t = ThermalScripted(f, 'LCF', frame_chid='frame', profile_pth=f, time_end=t_end)
                t.alias = chid
                t.change_in()
                run_safir(chid, safir_path)
                t.insert_tor(config_path)
                del t

        # Structural 3D analysis of the structure
        print('Running {} mechanical analysis...'.format(case.name))
        m = Mechanical(case.name, 'NF', frame_pth='frame')
        m.change_in()
        run_safir('frame', safir_path)

        print('[OK] {} scenario calculations finished!'.format(case.name))


if __name__ == '__main__':
    model = argv[1]
    if model == ('-s' or '--scripted'):
        print('=== ATTENTION! ===\nYou have entered scripted mode. It is designed to support developers - be careful,'
              'things don\'t work here well ;-)\n==================\n')
        paths = [input('safir_path = '), input('config_path = '), input('results_path = ')]
        scripted(*paths)
    else:
        calc_type = argv[2]
        path = argv[3]
        main(model, calc_type, path)
