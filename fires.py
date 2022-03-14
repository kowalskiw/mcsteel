from numpy import random, pi
from pandas import read_csv, merge, DataFrame, Series
from math import exp
from steputils import p21
from os import scandir

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


class FuelOBJ(Fuel):

    def find_points(self, **kwargs):
        file_lines = kwargs['obj_file']
        volume = []
        volumes = []
        vertices = []
        is_grouped = False
        def save_verts(volume): volumes.append([vertices[int(v_no)-1] for v_no in volume])

        for l in file_lines:
            if l[0] == 'v' and l[1] != 'n':
                vertices.append([float(v)/kwargs['scale'] for v in l.split()[1:]])  # add vertice coords to the list
            if l.startswith('g'):       # find volume
                is_grouped = True
                save_verts(volume)

                volume.clear()
            if l.startswith('f'):       # find face
                for v in l.split()[1:]:
                    vertice = v.split('//')[0]
                    if vertice not in volume:   # add face vertice to volume list if not already present
                        volume.append(vertice)
        if not is_grouped:
            raise RuntimeError('[ERROR] There is no group in OBJ file. It is required to group all volume boxes using'
                               '"g" at the line beginning. For further information about Wavefront OBJ format and '
                               'grouping see: http://fegemo.github.io/cefet-cg/attachments/obj-spec.pdf [2021-09-16]')
        save_verts(volume)

        return volumes[1:]

    def read_fuel(self):
        fuel = []
        for lay in scandir():
            splt = lay.name.split('.')
            if splt[-1] == ('obj'):
                with open(lay.name) as file:
                    obj = file.readlines()
                    for volume in self.find_points(obj_file=obj, scale=1000):
                        fuel.append([splt[0],  self.pts2fds(volume)])
                if len(fuel) == 0:
                    raise RuntimeError('[ERROR] No fuel data imported from {}'.format(''.join(splt)))
                elif fuel[-1][0] != splt[0]:
                    raise RuntimeError('[ERROR] No fuel data imported from {}'.format(''.join(splt)))

        return self.merge_data(fuel)  # merge with .FUL config type


class OldFuel(Fuel):
    def read_fuel(self):
        print('[ERROR] This fuel format is not working yet')
        exit(-1)


'''Draw fire config from input distributions
operates on the fire types included in the class - to be changed'''


class AlfaT2:
    def __init__(self, t_end, alpha, hrrpua):
        self.t_end = t_end  # duration of simulation
        self.hrr_max = 5e7  # W model limitation of HRR
        self.d_max = 10     # [m] model limitation of diameter
        self.alpha = alpha
        self.hrrpua = hrrpua

    # t-squared fire
    def t_squared(self):
        hrr_tab = []
        diam_tab = []
        for i in range(0, 99):
            t = int(self.t_end * i/98)

            # calculate HRR, append
            hrr_tab.append([t, round(self.alpha * (t ** 2) * 1000, 4)])  # [time /s/, HRR /W/]
            if hrr_tab[-1][-1] > self.hrr_max:    # check if hrr_tab does not exceed model limitation
                hrr_tab[-1][-1] = self.hrr_max

            # calculate diameter, append
            diam_tab.append([t, 2 * (hrr_tab[-1][-1] / (self.hrrpua*1000*pi))**0.5])  # [time /s/, diameter /m/]
            if diam_tab[-1][-1] > self.d_max:    # check if diam_tab does not exceed model limitation
                diam_tab[-1][-1] = self.d_max

        return hrr_tab, diam_tab

    def change(self, hrr_tab, d_tab):
        # as is
        return hrr_tab, d_tab

    def burn(self):
        hrr_tab, diam_tab = self.t_squared()

        changed_tabs = [self.change(*tab) for tab in [hrr_tab, diam_tab]]

        return changed_tabs[0], changed_tabs[1]


class SprinkNoEff(AlfaT2):
    def __init__(self, t_end, alpha, hrrpua, sprink_time):
        super().__init__(t_end, alpha, hrrpua)
        self.sprink_time = sprink_time

    # modify t-squared curve to taking non-effective sprinklers into account
    def change(self, hrr_tab, d_tab):
        q_0 = round(self.alpha * (self.sprink_time ** 2) * 1000, 4)
        d_0 = (q_0 / (self.hrrpua * 1000 * pi)) ** 0.5
        for tab, lim in [(hrr_tab, q_0), (d_tab, d_0)]:
            for i in range(len(tab)):
                if tab[i][0] >= self.sprink_time:
                    tab[i] = [tab[i][0], lim]
        return hrr_tab, d_tab


class SprinkEff(AlfaT2):
    def __init__(self, t_end, alpha, hrrpua, sprink_time):
        super().__init__(t_end, alpha, hrrpua)
        self.sprink_time = sprink_time

    # modify t-squared curve to taking non-effective sprinklers into account
    def change(self, hrr_tab, d_tab):
        q_0 = round(self.alpha * (self.sprink_time ** 2) * 1000, 4)
        d_0 = (q_0 / (self.hrrpua * 1000 * pi)) ** 0.5
        for tab, lim in [(hrr_tab, q_0), (d_tab, d_0)]:
            lower_limit = round(0.15 * lim)
            for i in range(len(tab)):
                if tab[i][0] >= self.sprink_time:
                    value_red = round(lim * exp(-0.0024339414 * (tab[i][0] - self.sprink_time)), 4)
                    if value_red > lower_limit:
                        tab[i] = [tab[i][0], value_red]
                    else:
                        tab[i] = [tab[i][0], lower_limit]
        return hrr_tab, d_tab
