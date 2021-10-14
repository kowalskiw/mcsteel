import dxfgrabber as dxf
from sys import argv
import numpy as np
from math import ceil
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

'''Script genarates locafi_script.txt file which contains set of template fires
burned one by one to achieve certain t-squared curve

You have to specify template fire [TXT], possible fire locations (centeres of fire base circles) [DXF]
and fire growth factor [W/s^2] you tend to'''


class TendToTSquared:
    def __init__(self, ins: dict, time_step=1):
        self.coordinates_pth = ins['dxf']
        self.alpha = ins['alpha']  # [W/s^2]
        self.t_end, self.z_ceiling, self.diameter, self.rhr, self.origin, self.q_plateau = \
            self.basic_locafi(ins['lcf'])
        self.n_of_fires = []
        self.estimated_fc = np.array([])
        self.lim_fc = np.array([])
        self.on_fc = np.array([])
        self.time_step = time_step

    def basic_locafi(self, template_pth):
        origin = []
        z_ceiling = 0
        read = [False, False]
        diameter = ['DIAMETER\n\n', 'END_DIAM\n\n']
        rhr = ['RHR\n\n', 'END_RHR\n\n']

        with open(template_pth) as file:
            template = file.readlines()

        for l in template:
            if 'FIRE_POS' in l:
                origin = [float(i) for i in l.split()[1:]]
            elif 'Z_CEILING' in l:
                z_ceiling = float(l.split()[-1])
            elif 'DIAMETER' in l:
                read[0] = True
            elif 'END_DIAM' in l:
                read[0] = False
            elif 'RHR' in l:
                read[1] = True
            elif 'END_RHR' in l:
                read[1] = False

            elif read[0]:
                diameter.insert(-1, l)

            elif read[1]:
                rhr.insert(-1, l)

        q_tab = []
        for i in rhr[1:-1]:
            try:
                q_tab.append(float(i.split()[-1]))
            except IndexError:
                pass

        print('[OK] locafi template imported')

        for i in reversed(diameter):
            try:
                tend = int(i.split()[0])
                break
            except (IndexError, ValueError):
                pass

        return tend, z_ceiling, diameter, rhr, origin, max(*q_tab)

    def locations(self):
        coord_list = []

        dxffile = dxf.readfile(self.coordinates_pth)

        for ent in dxffile.entities:
            if ent.dxftype == 'POINT':
                coord_list.append(list(ent.point))

        print('[OK] fires coordinates read')

        return coord_list

    def sort_coords(self):
        norigin = np.array(self.origin)
        sortd = []
        latest = 1e90
        raw = self.locations()
        for c in raw:
            nc = np.array(c)
            dist = np.sqrt(np.sum((norigin - nc) ** 2, axis=0))
            if dist < latest:
                latest = dist
                sortd.insert(0, c)
            else:
                sortd.append(c)

        return sortd

    def plateau_method(self, corrector: float = 0.0):
        # find how many single_fc fires needs to be burning to achive lim_fc
        self.n_of_fires.clear()
        for t in range(self.t_end):
            x = (self.alpha / self.q_plateau) ** 0.5 * t
            self.n_of_fires.append(ceil(x * (x + 1 + corrector)))

    def actual(self, tab, dt):
        changed = [tab[0], '\t0.0 0.0']
        for i in list(tab[1:-1]):
            itab = i.split()
            try:
                if float(itab[0]) + dt > self.t_end:
                    last = changed[-1].split()
                    t_step = float(itab[0]) - float(last[0])
                    v_step = float(itab[1]) - float(last[1])
                    if t_step <= 0:
                        v_end = 0.0
                    else:
                        v_end = float(last[1]) + (self.t_end - float(last[0])) / t_step * v_step
                    changed.append('\t{} {}\n'.format(self.t_end, v_end))
                    break
                else:
                    changed.append('\t{} {}\n'.format(float(itab[0]) + dt, itab[1]))
            except IndexError:
                changed.append('\n')
        changed.append(tab[-1])

        return changed

    def make_lcfs(self):
        fires_coords = self.sort_coords()
        self.t_end = len(self.n_of_fires)

        new_lcf = ['NFIRE {}\n\n'.format(self.n_of_fires[-1])]

        l = 0
        for t in range(1, self.t_end):
            for i in range(self.n_of_fires[t] - self.n_of_fires[t - 1]):
                try:
                    new_lcf.extend(['FIRE_POS  {} {} {}\nZ_CEILING  {}\nPLUME_TYPE CONIC\n'.format(*fires_coords[l],
                                                                                                   self.z_ceiling),
                                    *self.actual(self.diameter, t - 1), *self.actual(self.rhr, t - 1)])
                except IndexError:
                    raise ValueError('Number of fires locations is lower than required. Check DXF file and fire growth'
                                     'factor [W/s^2]')
                l += 1

        with open('locafi_script.txt', 'w') as file:
            file.writelines(new_lcf)

        print('[OK] "locafi_script.txt" written')

    def gimme_fc(self):
        one_fc = lcf2array(self.rhr, one_d=True)[1:]
        lim_fc = np.array([self.alpha * t ** 2 for t in range(self.t_end)])
        estimated_fc = np.array([0.0] * self.t_end)
        actual = 0
        t = 0
        for i in self.n_of_fires[1:]:
            t += 1
            new = i - actual
            if new > 0:
                for j in range(new):
                    actual += 1
                    estimated_fc = np.add(estimated_fc, np.array([0.0] * t + one_fc)[:self.t_end])
        return one_fc, lim_fc, estimated_fc

    def rhr_charts(self):
        one_fc, lim_fc, estimated_fc = self.gimme_fc()

        # print(estimated_fc)
        # print('\n\n', lim_fc)

        # Create 2x2 sub plots
        gs = gridspec.GridSpec(2, 2)
        fig = plt.figure()
        fig.suptitle('travelling fire estimation')
        
        err = fig.add_subplot(gs[1, 0])  # row 0, col 0
        err.plot(range(self.time_step, self.t_end), [estimated_fc[t] - lim_fc[t] for t in range(self.time_step, self.t_end)])
        err.set(ylabel='RHR error, [W]', xlabel='Time, [s]')

        rel_err = fig.add_subplot(gs[1, 1])  # row 0, col 1
        rel_err.plot(range(self.time_step, self.t_end), [(estimated_fc[t] - lim_fc[t]) / lim_fc[t] for t in range(self.time_step, self.t_end)])
        rel_err.set(ylabel='Relative RHR error, [-]', xlabel='Time, [s]')

        fc = fig.add_subplot(gs[0, :])  # row 1, span all columns
        fc.plot(0.95 * lim_fc, color='#c30404ff', linestyle='--')
        glf, = fc.plot(1.05 * lim_fc, color='#c30404ff', linestyle='--', label='Lim Fire Curve +/-5%')
        lf, = fc.plot(lim_fc, color='#c30404ff', linestyle='dotted', label='Lim Fire Curve')
        ef, = fc.plot(estimated_fc, label='Estimated Fire Curve', color='#0080d7')
        fc.legend(handles=[lf, ef, glf])
        fc.set(ylabel='Rate of Heat Release, [W]', xlabel='Time, [s]')
        [i.grid(True) for i in [err, rel_err, fc]]

        # plt.show()
        plt.savefig('locafi_script.png')
        print('[OK] figure saved')

    def optimize_corr(self, precision, relative=True):
        delta_mean = [1]
        corr = 0
        stop = False
        multi = 1
        while abs(delta_mean[-1]) > precision and not stop:
            print('%s: %f' % ('delta_mean', delta_mean[-1]), end='\r')
            if delta_mean[-1] > 0:  # rhr too high
                corr -= multi * precision
            else:  # rhr too low
                corr += multi * precision
            self.plateau_method(corrector=corr)
            one_fc, lim_fc, estimated_fc = self.gimme_fc()
            if relative:
                delta_mean.append(sum([(estimated_fc[t] - lim_fc[t]) / lim_fc[t] for t in range(self.time_step, self.t_end)])
                                  / (self.t_end - 1))
            else:
                delta_mean.append(sum([estimated_fc[t] - lim_fc[t] for t in range(self.time_step, self.t_end)])
                                  / (self.t_end - 1))
            if True:
                multi *= 2

            try:
                if delta_mean[-1] * delta_mean[-2] < 0:
                    if multi >= -1:
                        multi /= 2
                else:
                    pass

                if multi < 1:
                    break

                c = 0
                for t in range(-1, -21, -1):
                    if delta_mean[t] * delta_mean[t-1] < 0:
                        c += 1
                    if c == 20:
                        stop = True
            except IndexError:
                pass
        print('[OK] optimization with modifier finished')


    def optimize_iter(self, bottom=-0.01, top=0.05, time_step=10):

        # funkcja iterująca sekunda po sekundzie, jeśli dQ'/Q'_lim < -0.01: przesuwa pożar sekundę wcześniej
        # aktualizuje krzywą i tabelę zmiany pożarów
        # rusza od początku
        # jeśli kryterium OK: print('{}//{}'.fomrat(t, self.t_end))
        np_of_fires = np.array(self.n_of_fires)
        one_fc, lim_fc, estimated_fc = self.gimme_fc()
        conv_time = 0
        t = time_step
        comeback = 0

        # def rem(comeback):
        #     if np_of_fires[t - comeback * time_step] <= 0:
        #         comeback += 1
        #     np_of_fires[t - comeback * time_step:] -= 1
        # def add():
        #     np_of_fires[t - time_step:] += 1
        def modify(back=None):
            print('cb ', comeback)
            newt = t + back * time_step
            if 0 > newt:
                do_move = False  # remove fire
                newt = 0
            elif newt > self.t_end:
                do_move = False  # remove fire
                newt = self.t_end
            else:
                do_move = True  # move fire

            if back < 0:  # rhr too high
                if not do_move:
                    np_of_fires[newt:] -= 1
                    print('removed one fire at {} s'.format(newt))
                else:
                    for i in range(newt, t):
                        back -= 1
                        np_of_fires[i] -= 1
                    print('one fire moved from {} s to {} s'.format(t, newt))

            elif back > 0:  # rhr too low
                if not do_move:
                    np_of_fires[newt:] += 1
                    print('added one fire at {} s'.format(newt))
                else:
                    for i in range(t, newt):
                        np_of_fires[i] += 1
                        back += 1
                    print('one fire moved from {} s to {} s'.format(newt, t))

            else:
                raise ValueError('"comeback" parameter cannot have value "0"')

            return back

        while True:
            delta_rhr = (estimated_fc[t] - lim_fc[t]) / lim_fc[t]
            if delta_rhr < bottom:  # rhr too low

                print('{}/{}'.format(t, self.t_end))
                print(lim_fc[t], estimated_fc[t])
                print('RHR difference: {}'.format(delta_rhr))

                comeback += 1
                comeback = modify(back=comeback)
                self.n_of_fires = list(np_of_fires)
                one_fc, lim_fc, estimated_fc = self.gimme_fc()
                t = time_step

            elif delta_rhr > top:  # rhr too high

                print('{}/{}'.format(t, self.t_end))
                print(lim_fc[t], estimated_fc[t])
                print('RHR difference: {}'.format(delta_rhr))

                comeback -= 1
                comeback = modify(back=comeback)
                self.n_of_fires = list(np_of_fires)
                one_fc, lim_fc, estimated_fc = self.gimme_fc()
                t = time_step

            else:  # rhr ok
                conv_time = t
                t += time_step
                comeback = 0

            if conv_time == self.t_end - 1:
                break


def lcf2array(tab, one_d=False):
    new = []
    for l in tab:
        ls = l.split()
        if len(ls) > 1:
            new.append([float(ls[0]), float(ls[1])])

    if one_d:
        newer = []
        for t in range(int(new[-1][0])):
            newer.append(np.interp(t, *zip(*new)))
        return newer

    else:
        return new


def array2lcf(tab, type: str):
    if type.lower() in ['rhr', 'hrr', 'q', 'q\'']:
        new = ['RHR\n', 'END_RHR\n']
    elif type.lower() in ['diam', 'diameter', 'd']:
        new = ['DIAMETER\n', 'END_DIAM\n']
    else:
        raise ValueError('"type" variable not valid')

    for row in tab:
        new.append('\t{} {}'.format(*row))

    return new


class TendToFire(TendToTSquared):
    def __init__(self, files: dict, fire_curve: list, time_step: int = 5):
        super().__init__(files)
        self.lim_fire = np.array(fire_curve)
        self.diameter = lcf2array(self.diameter)
        self.rhr = lcf2array(self.rhr)
        self.dt = time_step

    def sum_n_check(self):
        def check(t):
            fc_t = np.interp(t, zip(*fc))
            req_t = np.interp(t, zip(*self.lim_fire))
            frac = fc_t / req_t
            if abs(frac - 1) > 0.05:
                return False
            else:
                return True

        def act(t):
            fc_t = np.interp(t, zip(*fc))
            req_t = np.interp(t, zip(*self.lim_fire))
            frac = fc_t / req_t
            if frac < 1:
                pass
                # add
            else:
                pass
                # remove

        calc = True
        fc = self.rhr.copy()
        step = 0
        dt = self.dt

        while calc:
            if check(step):
                step += dt
            else:
                act(step - dt)
                check(step - dt / 2)

        for t in range(self.t_end, self.dt):
            pass


if __name__ == '__main__':
    args = {}

    for a in argv[1:]:
        if a.split('.')[-1].lower() == 'dxf':
            args['dxf'] = a

        elif a.split('.')[-1].lower() == 'txt':
            args['lcf'] = a
        else:
            try:
                args['alpha'] = float(a)
            except:
                raise ValueError('Invalid alpha value')

    try:
        with open('config.txt') as file:
            foo = [i.split() for i in file.readlines()]
        config = {}
        for i in range(len(foo)):
            config[foo[i][0]] = foo[i][1]

        print('[OK] config file imported')

        a = TendToTSquared(args, time_step=int(config['time_step']))
        a.plateau_method()
        if config['optima'].lower() in ['1', 'true', 'yes', 't', 'y']:
            if config['relative'].lower() in ['0', 'false', 'no', 'f', 'n']:
                a.optimize_corr(float(config['precision']), relative=False)
                print('false')
            else:
                a.optimize_corr(float(config['precision']))

    except FileNotFoundError:
        a = TendToTSquared(args)
        a.plateau_method()

    a.make_lcfs()
    a.rhr_charts()
