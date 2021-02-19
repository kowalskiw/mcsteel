from os import getcwd, listdir, chdir
import subprocess
from shutil import copyfile
from sys import argv

'''Safir_thermal2D analyses'''


class Thermal:
    def __init__(self, chid, model, frame_chid='default', profile_pth='default'):
        self.chid = chid[:-4]
        self.mode = model.startswith('ISO')
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

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open(self.profile_pth) as file:
            init = file.readlines()

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
                [init.insert(n + 1, i) for i in ['BEAM_TYPE {}\n'.format(self.beam_type()),
                                                 '{}.in\n'.format(self.frame)]]

            # change thermal load
            elif l.startswith('   F  '):  # heating boundaries
                if self.model == 'CFD':
                    init[n] = 'FLUX {}'.format('CFD'.join(l[4:].split('FISO')))
                elif self.model == 'LCF':
                    init[n] = 'FLUX {}'.format('LOCAFI'.join(l[4:].split('FISO')))
                # elif self.model == 'HSM':
                #     init[n] = 'FLUX {}'.format('HSM'.join(l.split('FISO')[1:]))

        # write changed file
        with open(self.profile_pth, 'w') as file:
            file.writelines(init)

    # searching information about number of element
    def beam_type(self):
        with open('{}.in'.format(self.frame)) as file:
            frame = file.readlines()

        # check if the profile is present in frame.in file
        try:
            frame.index('{}.tem\n'.format(self.alias))
        except ValueError:
            raise ValueError('{} profile is not present in frame.in'.format(self.alias))

        while not frame.pop(0).startswith(' NODOFBEAM'):
            pass

        return round((frame.index('{}.tem\n'.format(self.alias)) - 1) / 3) + 1

    # copying fire and frame file to section catalogue
    def copy_ess(self):
        copyfile('{0}.gid\{0}.in'.format('frame'), '{}\{}.gid\{}.in'.format(self.path, self.chid, 'frame'))

        if self.model == 'CFD':
            copyfile('cfd.txt', '{}.gid\locafi.txt'.format(self.chid))
        elif self.model == 'LCF':
            copyfile('locafi.txt', '{}.gid\locafi.txt'.format(self.chid))

    # adding torsion analysis results to first TEM file
    def insert_tor(self):
        chdir('{}\{}.gid'.format(self.path, self.chid))
        with open('{}-1.T0R'.format(self.chid)) as file:
            tor = file.readlines()

        # nfiber = int(tor[0].split(' ')[-1])

        # picking TEM file to insert torsion results to
        if self.mode:
            first_b = self.chid + '.tem'
        else:
            for f in listdir(getcwd()):
                if f.endswith('.tem') and f.startswith('b'):
                    first_b = f
                    break

        # check if torsion results already are in TEM file
        with open(first_b) as file:
            tem = file.readlines()
        if '         w\n' in tem:
            chdir('..')
            print('Torsion results are already copied to the TEM')
            return 0

        # looking for start of torsion results regexp in TEM file
        try:
            tor_index = tor.index('         w\n')
        except ValueError:
            raise ValueError("Torsion results not found in the TOR")

        # find TEM line where torsion results should be passed
        try:
            if self.mode:
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
        chdir('..')
        print('Torsion results copied to the TEM')
        return 0

    # running single SAFIR simulation
    def run(self):
        # iso fire curve
        if self.mode:
            run_safir(self.chid)
            self.insert_tor()
            print('\nPrimary {} thermal 2D analysis finished\n\n'.format(self.chid))

        # natural fire
        else:
            self.copy_ess()
            try:
                self.change_in()
                run_safir(self.chid)
                self.insert_tor()
                print('\n{}-data {} thermal 2D analysis finished\n\n'.format(self.model, self.chid))
            except ValueError:
                raise ValueError("change_in not possible")


'''Safir_structural3D analyses'''


class Mechanical:
    def __init__(self, t2Ds, mode='ISO'):
        self.mode = mode.startswith('ISO')
        self.t2Ds = t2Ds

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open('frame.gid\{}.in'.format('frame')) as file:
            init = file.readlines()

        for n in range(len(init)):  # sections TEM files
            l = init[n]  # line of frame.in file
            if l.endswith('.tem\n'):
                for f in listdir('{}.gid'.format(l[:-5])):  # iterate over files in section dir
                    if f.endswith('.tem') and f.startswith('b'):
                        init[n] = '{}\n'.format(f)
                        break

        with open('{0}\{1}.gid\{1}.in'.format(getcwd(), 'frame'), 'w') as file:
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
        if self.mode:
            self.copy_tems()
            run_safir('frame')
            print('\nPrimary structural 3D analysis finished\n\n')

        # natural fire
        else:
            self.change_in()
            self.copy_tems()
            run_safir('frame')
            print('\nNatural fire structural 3D analysis finished\n\n')


# running SAFIR simulation in shell
def run_safir(chid, safir='C:\SAFIR', mcsteel=False):
    safir_path = '{}\safir.exe'.format(safir)

    if mcsteel:
        chid = chid.split('.')[0]
    else:
        chdir('.\{}.gid'.format(chid))

    subprocess.call(' '.join([safir_path, '"{}"'.format(chid)]), shell=True)

    if not mcsteel: chdir('..')


def main(model, type='s3dt2d', path=getcwd()):
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

    # Mechanical(folders).run()       # frame structural analysis - ISO curve

    if model == 'LCF' and 't2d' in type:
        for prof in folders:
            Thermal(prof, model).run()  # natural fire mode

    if model == 'LCF' and 's3d' in type:
        Mechanical(folders, mode='NF').run()  # frame structural analysis - natural fire

    print('\nAll SAFIR calculations finished, well done engineer!\n\n')

    if path != wd:
        chdir(wd)


if __name__ == '__main__':
    main(*argv[1:])
