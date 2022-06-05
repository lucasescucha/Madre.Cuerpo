import numpy as np
import math

from scipy.spatial.transform.rotation import Rotation

COS_30 = math.sqrt(3)/2
SIN_30 = 1/2

def getUnityVector(v):
    return v/np.linalg.norm(v)

def rotateAroundAxis(v, axis, angle):
    return Rotation.from_rotvec(angle * getUnityVector(axis)).apply(v)

def getRotationAngleAndAxis(initialVector, directionVector):
    initialVector = getUnityVector(initialVector)
    directionVector = getUnityVector(directionVector)

    axis = np.cross(initialVector, directionVector)
    
    return axis, np.arccos(np.clip(np.dot(initialVector, directionVector), -1.0, 1.0))
  
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
