from numpy import random, pi
from pandas import read_csv, merge, DataFrame
from math import exp
from steputils import p21


# triangular distribution sampler
def triangular(left, right, mode=False):
    if not mode:
        mode = (right - left) / 3 + left
    return random.triangular(left, mode, right)


# drawn fire site among fuel area
def mc_rand(csv):
    ases = []   # list with partial factors A of each fuel area
    probs = []  # list with probabilities of ignition in each fuel area

    # calculate partial factor A (volume * probability) of each fuel area
    for i, r in csv.iterrows():
        a = (r['XB'] - r['XA']) * (r['YB'] - r['YA']) * (r['ZB'] - r['ZA']) * r['MC']
        ases.append(a)

    # calculate probability of ignition in each fuel area
    for a in ases:
        probs.append(a/sum(ases))

    # return drawn fuel area
    return random.choice(len(probs), p=probs), ases


# drawn fire localization
# cwd == confg_path
def f_localization(title, config='fuel&stp'):
    if config == 'fuel&stp':
        ffile = Fuel(title).read_fuel()  # read FUL file
    else:
        ffile = read_csv('{}.ful'.format(title))
    fire_site = mc_rand(ffile)[0]  # generate fire coordinates from MC function
    config = ffile.iloc[fire_site]  # information about chosen fuel site

    def random_position(xes, yes, zes):
        coordinates = []
        [coordinates.append(random.randint(int(10 * i[0]), int(10 * i[1])) / 10) for i in [xes, yes, zes]]
        return coordinates

    return config, random_position((config.XA, config.XB), (config.YA, config.YB), zes=(config.ZA, config.ZB))


'''Import fuel configuration from STEP (geometrical distribution) and FUEL (data)
- to be merged with Properties config

cwd == config_path'''


class Fuel:
    def __init__(self, title):
        self.title = title
        self.step = p21.readfile('{}.stp'.format(title))

    # return points of the solid
    def find_points(self, volume):
        def params(ref, index=1):
            par = self.step.get(ref).entity.params

            while '*' in par[index]:
                index += 1

            if type(par[index]) == p21.Reference:
                return [par[index]]
            else:
                return [*par[index]]

        def add_point(new_point):
            if new_point not in points:
                points.append(new_point)

        points = []

        for i in params(volume.ref):    # manifold
            for j in params(i):  # closed shells
                for k in params(j):  # advanced face
                    for l in params(k):  # face outer bound
                        for m in params(l):  # edge loop
                            for n in params(m):  # oriented edge
                                for o in params(n, index=1):  # edge curve -- first point
                                    for p in params(o):  # vertex_point
                                        add_point(params(p))  # cartesian_point
                                for o in params(n, index=2):  # edge curve -- first point
                                    for p in params(o):  # vertex_point
                                        add_point(params(p))  # cartesian_point

        return points

    # return all volumes present in step file
    def read_step(self):

        vols = []
        layers = []
        for i in self.step:
            if type(i) == p21.SimpleEntityInstance:
                name = i.entity.name
                if name == 'MANIFOLD_SOLID_BREP':    # find volumes
                    vols.append(i)
                elif name == 'PRESENTATION_LAYER_ASSIGNMENT':            # find layers
                    layers.append(i)

        return vols, layers

    # return solid point entities in [[XA, XB], [YA, YB], [ZA, ZB]]format
    def pts2fds(self, pts):
        ptset = []
        xyz = [[], [], []]
        [[xyz[i].append(p[i]) for i in range(3)] for p in pts]
        [ptset.extend([min(i), max(i)]) for i in xyz]

        return ptset

    # return layer name: fuelX, where X is an integer correspondent to .FUL index
    def layer(self, ref, layers):
        # find layers names
        for l in layers:
            params = l.entity.params
            if ref in params[-1]:
                return params[0]
        return False

    def merge_data(self, vols):
        merged = DataFrame(columns=['XA', 'XB', 'YA', 'YB', 'ZA', 'ZB', 'MC', 'hrrpua_min', 'hrrpua_max', 'hrrpua_mode',
                                    't_sprink'])

        fuel = read_csv('{}.fuel'.format(self.title))

        for v in vols:
            data = fuel[fuel.name == v[0]]
            # convert list to DF
            coords = DataFrame([v[0], *v[1]], index=['name', 'XA', 'XB', 'YA', 'YB', 'ZA', 'ZB']).T
            # merge DF1 and DF2
            merged = merged.append([merge(coords, data).drop('name', axis=1)])

        return merged

    def read_fuel(self):
        fuel = []
        vols, lays = self.read_step()
        for v in vols:
            pts = self.find_points(v)
            fuel.append([self.layer(v.ref, lays), self.pts2fds(pts)])

        # merge with .FUL config type
        return self.merge_data(fuel)


'''Draw fire properties from input distributions
operates on the fire types included in the class - to be changed'''


class Properties:
    def __init__(self, t_end):
        self.t_end = t_end  # duration of simulation
        self.hrr_max = 3e8  # W model limitation of HRR

    # NOT WORKING - prepare to use

    # def newzealand(self, name):
    #     fuel_height = (0.32, 34.1)
    #     fuel_xes = (0.3, 23.1)
    #     fuel_yes = (10.3, 101.7)
    #
    #     H = fuel_height[1] - fuel_height[0]
    #     A_max = (fuel_xes[1] - fuel_xes[0]) ** 2 * 3.1415 / 4
    #
    #     config = rcsv('{}.ful'.format(name), sep=',')
    #     alpha = triangular(*config.alpha_min, *config.alpha_max, mode=float(config.alpha_mode))
    #     area = triangular(0, A_max)
    #
    #     print('alpha:{}, radius: {}'.format(alpha, (area / 3.1415) ** 0.5))
    #     hrr = []
    #     for i in range(0, int(self.t_end/99)):
    #         hrr.extend([i, round(H * alpha * (i ** 3) * 1000, 4)])  # [time /s/, HRR /W/]
    #         if hrr[-1] > self.hrr_max:
    #             hrr[-1] = self.hrr_max
    #
    #     return hrr, area, fuel_height, fuel_xes, fuel_yes
    #
    # def pool_fire(self, title, only_mass=False):
    #     with open('{}.ful'.format(title)) as file:
    #         fuel_prop = file.readlines()[1].split(',')
    #
    #     # random mass of fuel
    #     try:
    #         mass = triangular(int(fuel_prop[5]), int(fuel_prop[6]))
    #     except ValueError:
    #         mass = int(fuel_prop[5])
    #
    #     # random area of leakage
    #     if only_mass:
    #         area_ = mass * 0.03  # 0.019 # glycerol # 0.03 methanol leakage
    #         area = triangular(area_ * 0.9 * 100, area_ * 1.1 * 100) / 100
    #     else:
    #         try:
    #             area = triangular(int(fuel_prop[3]), int(fuel_prop[4]))
    #         except ValueError:
    #             area = int(fuel_prop[3])
    #
    #     if area < 0.28:
    #         ml_rate = triangular(0.015 * .9, 0.015 * 1.1)
    #     elif area < 7.07:
    #         ml_rate = triangular(0.022 * .9, 0.022 * 1.1)
    #     else:
    #         ml_rate = triangular(0.029 * .9, 0.029 * 1.1)
    #
    #     if self.a_max < area:
    #         area = self.a_max
    #
    #     print('mass loss rate = {}'.format(ml_rate))
    #     hrr_ = float(fuel_prop[1]) * ml_rate * area  # [MW] - heat release rate
    #     hrr = triangular(hrr_ * .8, hrr_ * 1.2)
    #
    #     time_end = mass / ml_rate / area
    #     if time_end > self.t_end:
    #         time_end = self.t_end
    #         hrr_list = [0, hrr, time_end / 60, hrr]
    #     else:
    #         if time_end < 60:
    #             time_end = 60
    #         hrr_list = [0, hrr, time_end / 60, hrr]
    #         hrr_list.extend([hrr_list[-2] + 1 / 6, 0, self.t_end / 60, 0])
    #
    #     print('HRR = {}MW'.format(hrr))
    #
    #     fuel_h = round(1 / float(fuel_prop[2]) / float(fuel_prop[5]), 2)
    #
    #     return hrr_list, area, fuel_h

    def test_fire(self):
        hrr = [0, 0, 15, 40]
        area = 10
        height = 0

        return hrr, area, height

    # t-squared fire
    def alfa_t2(self, name, property=None):
        config = f_localization(name)[0]

        # calculate HRRPUA according to triangular distribution specified by user
        hrrpua = triangular(config.hrrpua_min, config.hrrpua_max, mode=config.hrrpua_mode) * 1000     # kW/m2

        # calculate ALPHA according to the experimental log-norm or user's triangular distribution
        if not property:
            alpha = triangular(config.alpha_min, config.alpha_max, mode=config.alpha_mode)      # kW/s2
        elif property == 'store':
            alpha = hrrpua * random.lognormal(-9.72, 0.97)       # kW/s2

        hrr_tab = []
        diam_tab = []
        for i in range(0, 99):
            t = int(self.t_end * i/98)

            # calculate HRR, append
            hrr_tab.append([t, round(alpha * (t ** 2) * 1000, 4)])  # [time /s/, HRR /W/]
            if hrr_tab[-1][-1] > self.hrr_max:    # check if hrr_tab does not exceed model limitation
                hrr_tab[-1][-1] = self.hrr_max

            # calculate diameter, append
            diam_tab.append([t, 2 * (hrr_tab[-1][-1] / (hrrpua*1000*pi))**0.5])  # [time /s/, diameter /m/]

        return hrr_tab, diam_tab, hrrpua, alpha

    # curve taking non-effective sprinklers into account
    def sprink_noeff(self, name, property=None):
        config = f_localization(name)[0]

        # calculate HRRPUA according to triangular distribution specified by user
        hrrpua = triangular(config.hrrpua_min, config.hrrpua_max, mode=config.hrrpua_mode) * 1000  # [kW]

        # calculate ALPHA according to the experimental log-norm or user's triangular distribution
        if not property:
            alpha = triangular(config.alpha_min, config.alpha_max, mode=config.alpha_mode)  # [kW/s2]
        elif property == 'store':
            alpha = hrrpua * random.lognormal(-9.72, 0.97)   # [kW/s2]

        hrr_tab = []
        diam_tab = []
        for i in range(0, 99):
            t = int(self.t_end * i/98)

            # calculate HRR (steady-state after sprinklers activation), append
            if t < config.t_sprink:
                hrr_tab.append([t, round(alpha * (t ** 2) * 1000, 4)])  # [time /s/, HRR /W/]
            else:
                hrr_tab.append([t, round(alpha * (config.t_sprink ** 2) * 1000, 4)])  # [time /s/, HRR /W/]

            if hrr_tab[-1][-1] > self.hrr_max:  # check if hrr_tab does not exceed model limitation
                hrr_tab[-1][-1] = self.hrr_max

            # calculate diameter, append
            diam_tab.append([t, 2 * (hrr_tab[-1][-1] / (hrrpua*1000*pi))**0.5])  # [time /s/, diameter /m/]

        return hrr_tab, diam_tab, hrrpua, alpha

    # curve taking effective sprinklers into account
    def sprink_eff(self, name, property=None):
        config = f_localization(name)[0]

        # calculate HRRPUA according to triangular distribution specified by user
        hrrpua = triangular(config.hrrpua_min, config.hrrpua_max, mode=config.hrrpua_mode) * 1000  # [kW]

        # calculate ALPHA according to the experimental log-norm or user's triangular distribution
        if not property:
            alpha = triangular(config.alpha_min, config.alpha_max, mode=config.alpha_mode)  # [kW/s2]
        elif property == 'store':
            alpha = hrrpua * random.lognormal(-9.72, 0.97)   # [kW/s2]

        hrr_tab = []
        diam_tab = []
        q_0 = alpha * config.t_sprink ** 2 * 1000  # [W]
        q_limit = round(0.15 * q_0)
        for i in range(0, 99):
            t = int(self.t_end * i/98)

            # calculate HRR (extinguishing phase after sprinklers activation), append
            if t < config.t_sprink:
                hrr_tab.append([t, round(alpha * (t ** 2) * 1000, 4)])  # [time /s/, HRR /W/]
            else:
                q = round(q_0 * exp(-0.0024339414 * (t - config.t_sprink)), 4)  # W
                if q >= q_limit:
                    hrr_tab.append([t, q])  # [time /s/, HRR /W/]
                else:
                    hrr_tab.append([t, q_limit])

            if hrr_tab[-1][-1] > self.hrr_max:  # check if hrr_tab does not exceed model limitation
                hrr_tab[-1][-1] = self.hrr_max

            # calculate diameter, append
            diam_tab.append([t, 2 * (hrr_tab[-1][-1] / (hrrpua*1000*pi))**0.5])  # [time /s/, diameter /m/]

        return hrr_tab, diam_tab, hrrpua, alpha
