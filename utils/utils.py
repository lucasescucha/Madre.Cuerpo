import numpy as np
import math

import warnings

warnings.filterwarnings("error")

from scipy.spatial.transform.rotation import Rotation

COS_30 = math.sqrt(3)/2
SIN_30 = 1/2

def planeLineIntersection(planePoint, planeNormal, linePoint, lineDirection):
    if np.dot(planeNormal, getUnityVector(lineDirection)) == 0:
        return None
    
    t = (np.dot(planeNormal, planePoint) - np.dot(planeNormal, linePoint)) / np.dot(planeNormal, getUnityVector(lineDirection))
    return (getUnityVector(lineDirection) * t) + linePoint

def getUnityVector(v):
    try:
        return v/np.linalg.norm(v)
    except RuntimeWarning:
        print("Error")

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

def generateNutPolygon(s):
    COS_30 = math.sqrt(3)/2
    e = s/COS_30
    radius = e/2 

    return generateRegularPolygon(6, radius)

def generateRegularPolygon(sides, radius, rotation=0, z=0):
    one_segment = math.pi * 2 / sides

    return [
        [math.sin(one_segment * i + rotation) * radius,
         math.cos(one_segment * i + rotation) * radius,
         z] for i in range(sides + 1)]

class Configuration(dict):
    def __init__(self, dictionary: dict) -> None:
        for k in dictionary:
            v = dictionary[k]
            self[k] = Configuration(v) if isinstance(
                v, dict) else parseConfigurationValue(v)

    def __getattr__(self, attr: str) -> object:
        return self.get(attr)

    def check(self) -> bool:
        return True

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

