import numpy as np
import math
import matplotlib.pyplot as plt

from stl import mesh

COS_30 = math.sqrt(3)/2
SIN_30 = 1/2


def getVerticesAndFacesFromMesh(mesh):
    vertices = []
    faces = []

    for triangle in mesh:
        face = []
        for vertex in triangle:
            if not vertex in vertices:
                vertices.append(vertex)

            face.append(vertices.index(vertex))
        faces.append(face)

    return np.array(vertices), np.array(faces)


def getXYFromPoints(points):
    x = np.array([p[0] for p in points])
    y = np.array([p[1] for p in points])

    return x, y


def createNutHousingPolygon(margin, s):
    # e/2 * cos(30°) = s/2 => e = s / cos(30°)
    # https://www.fastenerdata.co.uk/media/wysiwyg/technical/FULLNUT.jpg
    e = s/COS_30
    e_2 = e/2
    return [[0, s/2],
            [margin+e_2*(1 + SIN_30), s/2],
            [e+margin, 0],
            [margin+e_2*(1 + SIN_30), -s/2],
            [0, -s/2]]


def saveTriangleMeshToSTL(trianglesMesh, filename):
    vertices, faces = getVerticesAndFacesFromMesh(trianglesMesh)

    _mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            _mesh.vectors[i][j] = vertices[f[j], :]

    # Write the mesh to file "cube.stl"
    _mesh.save(filename + ".stl")


class Configuration(dict):
    def __init__(self, dictionary: dict) -> None:
        for k in dictionary:
            v = dictionary[k]
            self[k] = Configuration(v) if isinstance(
                v, dict) else parseConfigurationValue(v)

    def __getattr__(self, attr: str) -> object:
        return self.get(attr)


def parseConfigurationValue(value):
    if not isinstance(value, str):
        return value

    def isfloat(x):
        try:
            float(x)
        except (TypeError, ValueError):
            return False
        else:
            return True

    if isfloat(value[:-1]):
        if value.endswith("m"):
            return float(value[:-1]) * 1000
        else:
            raise NotImplementedError
    else:
        return value
