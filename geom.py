import gmsh
import numpy as np


class Geometry:
    def __init__(self):
        gmsh.initialize()
        self.mesh_size = 0.1    # element size [m]
        self.nodes = list()
        self.lines = list()
        self.shells = list()
        self.volumes = list()

        self.sections = list()

    def get_node_id(self, idx: int): return self.nodes[idx]

    def get_node_tag(self, tag: int): return self.nodes[tag-1]

    def get_line_id(self, idx: int): return self.lines[idx]

    def get_line_tag(self, tag: int): return self.lines[tag-1]

    def get_line_nodes(self, tag: int): return [self.nodes[node_tag] for node_tag in self.lines[tag]]

    def get_line_startpoint(self, tag: int): return self.nodes[self.lines[tag][0]]

    def get_line_endpoint(self, tag: int): return self.nodes[self.lines[tag][0]]

    def get_shell_id(self, idx: int): return self.shells[idx]

    def get_shell_tag(self, tag: int): return self.shells[tag-1]

    # maybe we don't need to mesh at all
    def _mesh(self):
        # exclude lines that are used for surfaces from meshing
        for surface in gmsh.model.getEntities(2):
            loop_tag, loop_lines = gmsh.model.occ.getCurveLoops(1)
            print(loop_lines)
            print(surface)
        # mesh lines only
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), self.mesh_size)
        gmsh.model.mesh.generate(1)

    def _import_nodes(self):
        # imoport nodes as sets of coordinates
        _, node_tags = zip(*gmsh.model.getEntities(0))
        self.nodes = [Node(*gmsh.model.getValue(0,tag, [0]), tag=tag) for tag in node_tags]

    def add_node(self, x, y, z, safe=False):
        if safe:
            for i, n in enumerate(self.nodes):
                if x==n.x and y==n.y and z==n.z:
                    return i

        self.nodes.append(Node(x, y, z))
        return len(self.nodes) - 1

    def _import_lines(self):
        _, line_tags = zip(*gmsh.model.getEntities(1))
        for tag in line_tags:
            # find start and endpoint nodes
            _, line_nodes_tags = zip(*gmsh.model.getBoundary([(1, tag)]))

            # find group name
            group_tag = gmsh.model.getPhysicalGroupsForEntity(1, tag)
            if group_tag.size == 0:
                group_name = None    # not physical lines
            else:
                group_name = gmsh.model.getPhysicalName(1, group_tag[0])
                if group_name not in self.sections:
                    self.sections.append(group_name)

            # add line
            current_nodes = [self.get_node_tag(n) for n in line_nodes_tags]    # nodes are already imported
            current_line = Line(*current_nodes, group=group_name, tag=tag)
            self.lines.append(current_line)

    def _import_shells(self):
        _, shell_tags = zip(*gmsh.model.getEntities(2))
        for tag in shell_tags:
            # find edges
            _, shell_lines_tags = zip(*gmsh.model.getBoundary([(2, tag)]))

            # add shell
            current_lines = [self.get_line_tag(t) for t in shell_lines_tags]
            current_shell = Shell(current_lines)
            self.shells.append(current_shell)

    def _import_volumes(self):
        _, vols_tags = zip(*gmsh.model.getEntities(3))
        for tag in vols_tags:
            # find edges
            _, vols_shells_tags = zip(*gmsh.model.getBoundary([(3, tag)]))

            # add shell
            current_shells = [self.get_shell_tag(t) for t in vols_shells_tags]
            current_vol = Volume(current_shells)
            self.shells.append(current_vol)


    # read geometry
    def read_geo_geom(self, geo_path: str):
        gmsh.initialize()
        gmsh.open(geo_path)
        gmsh.model.occ.synchronize()
        # self._mesh()
        self._import_nodes()
        self._import_lines()
        self._import_shells()
        self._import_volumes()

    def close_gmsh(self): gmsh.clear()

    def find_closest_section(self, fire_coords: list[float], physical_only=True):
        section_coords = {k: [(None, None, None), 0] for k in self.sections}

        # iterate over elements
        for l in self.lines:
            if physical_only and not l.group:
                continue
            projection = gmsh.model.getClosestPoint(1, l.tag, fire_coords)
            if not self.points_between_shells(fire_coords, projection):
                continue
            distance = np.linalg.norm(np.array(fire_coords) - np.array(projection))
            if section_coords[l.group][1] > distance:
                section_coords[l.group] = [projection, distance]

        # return {k:v[0] for k,v in section_coords.items()}
        return section_coords

    def points_between_shells(self, p1, p2):
        if all([(sh.level - p1[2]) * (sh.level - p2[2]) > 0 for sh in self.shells]):
            return True
        else:
            return False


class Node:
    def __init__(self, x: float, y: float, z: float, tag: int = None):
        self.x = x
        self.y = y
        self.z = z
        self.tag = tag

    def __repr__(self): return f'{self.x} {self.y} {self.z}'

    def __str__(self): return f'Node ({self.x} {self.y} {self.z})'


class Line:
    def __init__(self, start: Node, end: Node, group: str = None, tag: int = None):
        self.nodes = [start, end]
        self.group = group
        self.tag = tag

    def __repr__(self): return f'{self.nodes}'

    def __str__(self): return f'Line ({self.nodes})'


class Shell:
    def __init__(self, lines: list[Line], tag: int = None):
        self.lines = lines
        self.level = self._get_level()

    def __repr__(self): return f'{self.level}'

    def __str__(self): return f'Shell ({self.level})'

    def _get_level(self):
        level = 1e9     # super-high
        for l in self.lines:
            level = min([level] + [n.z for n in l.nodes])

        return level


class Volume:
    def __init__(self, shells: list[Shell], group: str = None, tag: int = None):
        self.shells = shells
        self.group = group
        self.tag = tag

    def __repr__(self): return f'{self.group} {self.shells}'

    def __str__(self): return f'Volume ({self.group} {self.shells})'


def read_geo_geom(geo_path):
    if not geo_path.endswith('.geo'):
        geo_path += '.geo'
    geometry = Geometry()
    geometry.read_geo_geom(geo_path)
    return geometry


def find_closest_section(fire_coords: list[float], structure: Geometry):
    return structure.find_closest_section(fire_coords)

# if __name__ == '__main__':
#     geom = Geometry()
#     geom.read_geo_geom('/Users/wk/Library/CloudStorage/OneDrive-AkademiaPożarnicza/01_instytut/2023/05_aamks_for_mezzanines/_analiza/test_gmsh.geo')
#     # geom = read_step_geom('/Users/wk/Library/CloudStorage/OneDrive-AkademiaPożarnicza/01_instytut/2023/05_aamks_for_mezzanines/_analiza/part_lines_legacy.step')
#     # geom = read_iges_geom('/Users/wk/Library/CloudStorage/OneDrive-AkademiaPożarnicza/01_instytut/2023/05_aamks_for_mezzanines/_analiza/part_lines.iges')
#     # lines = select_lines_in_floor(geom, None)
#     pass