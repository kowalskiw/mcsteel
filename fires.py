from numpy import random, pi
from pandas import read_csv, merge, DataFrame
from math import exp
from steputils import p21
from os import scandir


# triangular distribution sampler
def triangular(left, right, mode=False):
    if not mode:
        mode = (right - left) / 3 + left
    return random.triangular(left, mode, right)


# drawn fire site among fuel area
def mc_rand(csv):
    ases = []   # list with partial factors A of each fuel area
    probs = []  # list with probabilities of ignition in each fuel area

    # calculate partial factor A (area * probability) of each fuel area
    for i, r in csv.iterrows():
        a = (r['XB'] - r['XA']) * (r['YB'] - r['YA']) * r['MC']
        ases.append(a)

    # calculate probability of ignition in each fuel area
    for a in ases:
        probs.append(a/sum(ases))

    # return drawn fuel area
    return random.choice(len(probs), p=probs)


# drawn fire localization
# cwd == confg_path
def f_localization(ffile):
    def random_position(xes, yes, zes):
        coordinates = []
        [coordinates.append(random.randint(int(10 * i[0]), int(10 * i[1])) / 10) for i in [xes, yes, zes]]
        return coordinates

    fire_site = mc_rand(ffile)  # generate fire coordinates from MC function
    config = ffile.iloc[fire_site]  # information about chosen fuel site

    return config, random_position((config.XA, config.XB), (config.YA, config.YB), zes=(config.ZA, config.ZB))


'''Import fuel configuration from STEP (geometrical distribution) and FUEL (data)
- to be merged with Properties config

cwd == config_path'''


class Fuel:
    def __init__(self, title):
        self.title = title
        self.layers = []

    # return points of the solid
    def find_points(self, volume, step):
        def params(ref, index=1):
            par = step.get(ref).entity.params

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
    def read_step(self, step):
        vols = []

        for i in step:
            if type(i) == p21.SimpleEntityInstance:
                name = i.entity.name
                if name == 'MANIFOLD_SOLID_BREP':    # find volumes
                    vols.append(i)

        return vols

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
        for lay in scandir():
            splt = lay.name.split('.')
            if splt[-1] == ('step' or 'stp'):
                step = p21.readfile(lay.name)
                vols = self.read_step(step)
                for v in vols:
                    pts = self.find_points(v, step)
                    fuel.append([splt[0], self.pts2fds(pts)])
        # merge with .FUL config type
        return self.merge_data(fuel)


'''Draw fire config from input distributions
operates on the fire types included in the class - to be changed'''


class Properties:
    def __init__(self, t_end, properties):
        self.t_end = t_end  # duration of simulation
        self.hrr_max = 3e8  # W model limitation of HRR
        self.config = properties

    # def test_fire(self):
    #     hrr = [0, 0, 15, 40]
    #     area = 10
    #     height = 0
    #
    #     return hrr, area, height

    # t-squared fire
    def alfa_t2(self, property=None):

        # calculate HRRPUA according to triangular distribution specified by user
        hrrpua = triangular(self.config.hrrpua_min, self.config.hrrpua_max, mode=self.config.hrrpua_mode)     # [kW/m2]

        # calculate ALPHA according to the experimental log-norm or user's triangular distribution
        if not property:
            alpha = triangular(self.config.alpha_min, self.config.alpha_max, mode=self.config.alpha_mode)      # [kW/s2]
        elif property == 'store':
            alpha = hrrpua * random.lognormal(-9.72, 0.97)       # [kW/s2]

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
    def sprink_noeff(self, property=None):

        # calculate HRRPUA according to triangular distribution specified by user
        hrrpua = triangular(self.config.hrrpua_min, self.config.hrrpua_max, mode=self.config.hrrpua_mode)  # [kW/m2]

        # calculate ALPHA according to the experimental log-norm or user's triangular distribution
        if not property:
            alpha = triangular(self.config.alpha_min, self.config.alpha_max, mode=self.config.alpha_mode)  # [kW/s2]
        elif property == 'store':
            alpha = hrrpua * random.lognormal(-9.72, 0.97)   # [kW/s2]

        hrr_tab = []
        diam_tab = []
        for i in range(0, 99):
            t = int(self.t_end * i/98)

            # calculate HRR (steady-state after sprinklers activation), append
            if t < self.config.t_sprink:
                hrr_tab.append([t, round(alpha * (t ** 2) * 1000, 4)])  # [time /s/, HRR /W/]
            else:
                hrr_tab.append([t, round(alpha * (self.config.t_sprink ** 2) * 1000, 4)])  # [time /s/, HRR /W/]

            if hrr_tab[-1][-1] > self.hrr_max:  # check if hrr_tab does not exceed model limitation
                hrr_tab[-1][-1] = self.hrr_max

            # calculate diameter, append
            diam_tab.append([t, 2 * (hrr_tab[-1][-1] / (hrrpua*1000*pi))**0.5])  # [time /s/, diameter /m/]

        return hrr_tab, diam_tab, hrrpua, alpha

    # curve taking effective sprinklers into account
    def sprink_eff(self, property=None):

        # calculate HRRPUA according to triangular distribution specified by user
        hrrpua = triangular(self.config.hrrpua_min, self.config.hrrpua_max, mode=self.config.hrrpua_mode)  # [kW]

        # calculate ALPHA according to the experimental log-norm or user's triangular distribution
        if not property:
            alpha = triangular(self.config.alpha_min, self.config.alpha_max, mode=self.config.alpha_mode)  # [kW/s2]
        elif property == 'store':
            alpha = hrrpua * random.lognormal(-9.72, 0.97)   # [kW/s2]

        hrr_tab = []
        diam_tab = []
        q_0 = alpha * self.config.t_sprink ** 2 * 1000  # [W]
        q_limit = round(0.15 * q_0)
        for i in range(0, 99):
            t = int(self.t_end * i/98)

            # calculate HRR (extinguishing phase after sprinklers activation), append
            if t < self.config.t_sprink:
                hrr_tab.append([t, round(alpha * (t ** 2) * 1000, 4)])  # [time /s/, HRR /W/]
            else:
                q = round(q_0 * exp(-0.0024339414 * (t - self.config.t_sprink)), 4)  # W
                if q >= q_limit:
                    hrr_tab.append([t, q])  # [time /s/, HRR /W/]
                else:
                    hrr_tab.append([t, q_limit])

            if hrr_tab[-1][-1] > self.hrr_max:  # check if hrr_tab does not exceed model limitation
                hrr_tab[-1][-1] = self.hrr_max

            # calculate diameter, append
            diam_tab.append([t, 2 * (hrr_tab[-1][-1] / (hrrpua*1000*pi))**0.5])  # [time /s/, diameter /m/]

        return hrr_tab, diam_tab, hrrpua, alpha
