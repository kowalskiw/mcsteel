from os import getcwd, listdir, chdir
import subprocess
from shutil import copyfile
from sys import argv

'''Safir_thermal2D analyses'''


class Thermal:
    def __init__(self, chid, model, mode='ISO'):
        self.chid = chid[:-4]
        self.mode = mode.startswith('ISO')
        self.model = model

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open('{0}\{1}.gid\{1}.in'.format(getcwd(), self.chid)) as file:
            init = file.readlines()

        for n in range(len(init)):
            l = init[n]

            if l == 'MAKE.TEM\n':       # type of calculation
                if self.model == 'CFD':
                    init[n] = 'MAKE.TEMCD\n'
                elif self.model == 'LCF':
                    init[n] = 'MAKE.TEMLF\n'
                [init.insert(n+1, i) for i in ['BEAM_TYPE {}\n'.format(self.beam_type()), 'frame.in\n']]

            elif l.startswith('   F  '):        # heating boundaries
                print(l.split('FISO'))
                if self.model == 'CFD':
                    init[n] = 'FLUX {}'.format('CFD'.join(l[4:].split('FISO')))
                elif self.model == 'LCF':
                    init[n] = 'FLUX {}'.format('LOCAFI'.join(l[4:].split('FISO')))
                    print(init[n])
                # elif self.model == 'HSM':
                #     init[n] = 'FLUX {}'.format('HSM'.join(l.split('FISO')[1:]))

        with open('{0}\{1}.gid\{1}.in'.format(getcwd(), self.chid), 'w') as file:
            file.writelines(init)

    # searching information about number of element
    def beam_type(self):
        with open('{}\{}.gid\{}'.format(getcwd(), self.chid, 'frame.in')) as file:
            frame = file.readlines()

        while not frame.pop(0).startswith(' NODOFBEAM'):
            pass

        return round((frame.index('{}.tem\n'.format(self.chid))-1)/3) + 1

    # copying fire and frame file to section catalogue
    def copy_ess(self):
        copyfile('{0}.gid\{0}.in'.format('frame'), '{}\{}.gid\{}.in'.format(getcwd(), self.chid, 'frame'))

        if self.model == 'CFD':
            copyfile('cfd.txt', '{}.gid\locafi.txt'.format(self.chid))
        elif self.model == 'LCF':
            copyfile('locafi.txt', '{}.gid\locafi.txt'.format(self.chid))

    # adding torsion analysis results to first TEM file
    def insert_tor(self):
        chdir('.\{}.gid'.format(self.chid))
        with open('{}-1.T0R'.format(self.chid)) as file:
            tor = file.readlines()

        nfiber = int(tor[0].split(' ')[-1])

        # looking for results regexp
        while not tor.pop(0).startswith(' ==='):
            pass

        # picking first TEM file
        for f in listdir(getcwd()):
            if f.endswith('.tem'):
                first_b = f
                break

        with open(first_b) as file:
            tem = file.readlines()
        if '         w\n' in tem:       # check if torsion results already are in TEM file
            chdir('..')
            return 0
        # looking for start of torsion results regexp in TEM file
        try:
            if self.model == 'CFD':
                index = tem.index('       CFD\n')  # if model == CFD
            elif self.model == 'LCF':
                index = tem.index('       LCF\n')  # if model == LCF
            # elif self.model == 'HSM':
            #     index = tem.index('       HSM\n')   # if model == HSM
        except ValueError:
            index = tem.index('       HOT\n')

        # pasting torsion results
        with open(first_b, 'w') as file:
            file.writelines(tem[:index] + tor[:nfiber+4] + tem[index:])
        chdir('..')

    # running single SAFIR simulation
    def run(self):
        # iso fire curve
        if self.mode:
            run_safir(self.chid)
            self.insert_tor()
            print('Primary {} thermal 2D analysis finished'.format(self.chid))

        # natural fire
        else:
            self.copy_ess()
            self.change_in()
            run_safir(self.chid)
            self.insert_tor()
            print('{}-data {} thermal 2D analysis finished'.format(self.model, self.chid))


'''Safir_structural3D analyses'''


class Mechanical:
    def __init__(self, t2Ds, mode='ISO'):
        self.mode = mode.startswith('ISO')
        self.t2Ds = t2Ds

    # changing input file form iso curve to natural fire mode
    def change_in(self):
        with open('frame.gid\{}.in'.format('frame')) as file:
            init = file.readlines()

        for n in range(len(init)):      # sections TEM files
            l = init[n]
            if l.endswith('.tem\n'):
                for f in listdir('{}.gid'.format(l[:-5])):
                    if f.endswith('.tem') and f.startswith('b'):
                        init[n] = '{}\n'.format(f)

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
            print('Primary structural 3D analysis finished')

        # natural fire
        else:
            self.change_in()
            self.copy_tems()
            run_safir('frame')
            print('{}-data structural 3D analysis finished'.format(self.model))
        

# running SAFIR simulation in shell
def run_safir(chid):
    safir_path = 'C:\SAFIR\safir.exe'
    chdir('.\{}.gid'.format(chid))

    subprocess.call(' '.join([safir_path, chid, ]), shell=True)
    chdir('..')


if __name__ == '__main__':
    model = argv[1]
    folders = listdir(getcwd())
    folders.remove('frame.gid')
    for f in reversed(folders):
        if not f.endswith('.gid'):
            folders.remove(f)

    ### for prof in folders:
    ###     Thermal(prof, model).run()

    Mechanical(folders).run()       # frame structural analysis - ISO curve

    for prof in folders:                        # calculating meshed sections
        Thermal(prof, model, mode='NF').run()   # natural fire mode

    Mechanical(folders, mode='NF').run()        # frame structural analysis - natural fire

    print('All SAFIR calculations finished, well done engineer!')
