from os import getcwd, listdir, chdir, scandir
import subprocess
from shutil import copyfile
from sys import argv


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


def sections(frame):
    with open('{}.in'.format(frame)) as file:
        frame_lines = file.readlines()

    tems = {}
    c = 1
    for n in range(len(frame_lines)):  # sections TEM files
        l = frame_lines[n]  # line of frame.in file
        # create dict {'profile_no':[section_name, section_line_no, 'bXXXXX.tem'], ...}
        if '.tem' in l:
            tems[str(c)] = [l[:-1], n]
            c += 1
        if '      ELEM' in l:
            element = l.split()
            if len(tems[element[-1]]) < 3:
                number = element[1]
                for i in range(5 - len(number)):
                    number = '0{}'.format(number)
                tems[element[-1]].append('b{}_1.tem'.format(number))

    return tems


'''Safir_thermal2D analyses'''


class Thermal:
    def __init__(self, chid, model, frame_chid='default', profile_pth='default'):
        self.chid = chid
        self.model = model
        self.path = getcwd()

        if frame_chid == 'default':
            self.frame = 'frame'
        else:
            self.frame = frame_chid

        if profile_pth == 'default':
            self.profile_pth = '{0}\{1}.gid\{1}.in'.format(getcwd(), self.chid)
            self.alias = self.chid
        else:
            self.profile_pth = profile_pth
            self.alias = 'safir'

        self.sections = sections(self.frame)

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open(self.profile_pth) as file:
            init = file.readlines()

        with open('{}.backup'.format(self.alias), 'w') as file:
            file.writelines(init)

        # make changes
        for n in range(len(init)):
            l = init[n]

            # type of calculation
            if l == 'MAKE.TEM\n':
                if self.model == 'CFD':
                    init[n] = 'MAKE.TEMCD\n'
                elif self.model == 'LCF':
                    init[n] = 'MAKE.TEMLF\n'

                # insert beam type
                for k, v in self.sections.items():
                    if v[0][:-4] == self.alias:
                        [init.insert(n + 1, i) for i in ['BEAM_TYPE {}\n'.format(k), '{}.in\n'.format(self.frame)]]

            # change thermal load
            elif l.startswith('   F  ') and 'FISO' in l:  # heating boundaries with FISO
                if self.model == 'CFD':
                    if 'F20' not in l:
                        init[n] = 'FLUX {}'.format('CFD'.join(l[4:].split('FISO')))
                    else:
                        init[n] = 'FLUX {}'.format('NO'.join(('CFD'.join(l[4:].split('FISO'))).split('F20')))
                        init.insert(n + 1, 'NO'.join(l.split('FISO')))
                elif self.model == 'LCF':
                    if 'F20' not in l:
                        init[n] = 'FLUX {}'.format('LOCAFI'.join(l[4:].split('FISO')))
                    else:
                        init[n] = 'FLUX {}'.format('NO'.join(('LOCAFI'.join(l[4:].split('FISO'))).split('F20')))
                        init.insert(n + 1, 'NO'.join(l.split('FISO')))
                # elif self.model == 'HSM':
                #     init[n] = 'FLUX {}'.format('HSM'.join(l.split('FISO')[1:]))

        # write changed file
        with open(self.profile_pth, 'w') as file:
            file.writelines(init)

    # searching information about number of element
    # def beam_type(self):
    #     with open('{}.in'.format(self.frame)) as file:
    #         frame = file.readlines()
    #
    #     # check if the profile is present in frame.in file
    #     try:
    #         frame.index('{}.tem\n'.format(self.alias))
    #     except ValueError:
    #         print('{} profile is not present in frame.in'.format(self.alias))
    #
    #     while not frame.pop(0).startswith(' NODOFBEAM'):
    #         pass
    #
    #     return round((frame.index('{}.tem\n'.format(self.alias)) - 1) / 3) + 1

    # copying fire and frame file to section catalogue
    def copy_ess(self):
        copyfile('{0}.gid\{0}.in'.format('frame'), '{}\{}.gid\{}.in'.format(self.path, self.chid, 'frame'))

        if self.model == 'CFD':
            copyfile('cfd.txt', '{}.gid\locafi.txt'.format(self.chid))
        elif self.model == 'LCF':
            copyfile('locafi.txt', '{}.gid\locafi.txt'.format(self.chid))

# dodaj first_b dla scripted
    # adding torsion analysis results to first TEM file
    def insert_tor(self, config_path='.'):
        # chdir('{}\{}.gid'.format(self.path, self.chid))
        with open('{0}\{1}.gid\{1}-1.T0R'.format(config_path, self.chid)) as file:
            tor = file.readlines()

        # picking TEM file to insert torsion results to
        if self.model == 'ISO':
            first_b = self.chid + '.tem'
        else:
            for v in self.sections.values():
                print(v[0][:-4], self.chid, v[0][:-4] == self.chid)
                if v[0][:-4] == self.chid:
                    first_b = v[-1]
                    print(first_b)
        # check if torsion results already are in TEM file
        try:
            with open(first_b) as file:
                tem = file.readlines()
            if '         w\n' in tem:
                # chdir('..')
                print('Torsion results are already copied to the TEM')
                return 0

            # looking for start of torsion results regexp in TEM file
            try:
                tor_index = tor.index('         w\n')
            except ValueError:
                raise ValueError("Torsion results not found in the TOR")

            # find TEM line where torsion results should be passed
            try:
                if self.model == 'ISO':
                    tem_index = tem.index('       HOT\n')  # if model == ISO
                elif self.model == 'CFD':
                    tem_index = tem.index('       CFD\n')  # if model == CFD
                elif self.model == 'LCF':
                    tem_index = tem.index('    LOCAFI\n')  # if model == LCF
                # elif self.model == 'HSM':
                #     tem_index = tem.tem_index('       HSM\n')   # if model == HSM
            except ValueError:
                raise ValueError("Flux constraint information not found in the TEM")

            # pasting torsion results
            with open(first_b, 'w') as file:
                file.writelines(tem[:tem_index] + tor[tor_index:-1] + tem[tem_index:])
            # chdir('..')
            print('Torsion results copied to the TEM')

            return 0

        except FileNotFoundError:
            raise FileNotFoundError('There is no proper TEM file in {}'.format(getcwd()))

        except UnboundLocalError:
            print('The {} profile is not used in the structure'.format(self.chid))
            return -1

    # running single SAFIR simulation
    def run(self):
        # iso fire curve
        if self.model == 'ISO':
            run_safir(self.chid)
            self.insert_tor()
            print('\nISO-curve {} thermal 2D analysis finished\n\n'.format(self.chid))

        # natural fire
        else:
            self.copy_ess()
            try:
                self.change_in()
                run_safir(self.chid)
                self.insert_tor()
                print('\n{}-data {} Thermal 2D analysis finished\n\n'.format(self.model, self.chid))
            except ValueError:
                raise ValueError("change_in not possible")


'''Safir_structural3D analyses'''


class Mechanical:
    def __init__(self, t2Ds, model, frame_pth='frame.gid/frame.in'):
        self.model = model  # ISO/other
        self.t2Ds = t2Ds
        self.frame = frame_pth

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open('{}.in'.format(self.frame)) as file:
            init = file.readlines()

        with open('frame.backup', 'w') as file:
            file.writelines(init)

        # change TEM names in IN file
        tems = sections(self.frame)
        for v in tems.values():
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
        if self.model == 'ISO':
            self.copy_tems()
            run_safir('frame')
            print('\nISO-curve Structural 3D analysis finished\n\n')

        # natural fire
        else:
            self.change_in()
            self.copy_tems()
            run_safir('frame')
            print('\nNatural fire Structural 3D analysis finished\n\n')


# running SAFIR simulation in shell
def run_safir(in_file, safir='C:\SAFIR', mcsteel=False):
    safir_path = '{}\safir.exe'.format(safir)

    if mcsteel:
        chid = in_file.split('.')[0]
    else:
        chdir('.\{}.gid'.format(in_file))

    subprocess.call(' '.join([safir_path, '"{}"'.format(in_file)]), shell=True)

    chdir('..') if not mcsteel else 0


def main(model, calc_type='s3dt2d', path='.'):
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

    if model == 'LCF' and 't2d' in calc_type:
        for prof in folders:
            Thermal(prof, model).run()  # natural fire mode

    if model == 'LCF' and 's3d' in calc_type:
        Mechanical(folders, 'NF').run()  # frame structural analysis - natural fire

    print('\nAll SAFIR calculations finished, well done engineer!\n\n')

    if path != wd:
        chdir(wd)


# cwd = worst\case
def scripted(safir_path, config_path, results_path):
    for case in scandir('{}\worst'.format(results_path)):
        chdir(case.path)

        # Thermal 2D analyses of profiles
        for i in scandir():
            f = i.name
            if f.endswith('.in') and not f == 'frame.in':
                chid = f[:-3]
                t = Thermal(chid, 'LCF', profile_pth=f)
                t.alias = chid
                t.change_in()
                run_safir(chid, safir_path, mcsteel=True)
                t.insert_tor(config_path)
                del t

        # Structural 3D analysis of the structure
        m = Mechanical(case.name, 'NF', frame_pth='frame')
        m.change_in()
        run_safir('frame', safir_path, mcsteel=True)

        print('\n{} case calculations finished!\n\n'.format(case.name))


if __name__ == '__main__':
    main(*argv[1:])

