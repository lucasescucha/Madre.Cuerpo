import numpy as np

from surfaces.iSurface import ISurface

class MeshSurface(ISurface):
    def __init__(self, verticesGrid: np.array, triangleMesh: np.array) -> None:
        self.verticesGrid = verticesGrid
        self.triangleMesh = triangleMesh

    def F(self, P: np.array) -> float:
        x, y = P
        
        v1, v2, v3 = np.where(checkPointInTriangle(P, self.triangleMesh))[0]
        a, b, c = np.cross(v2 - v1, v3 - v1)
        x0, y0, z0 = v1
        
        return (-1/c * (a * (x-x0) + b * (y-y0))) + z0
    
    def gradF(self, P: np.array) -> np.array:
        raise NotImplementedError
    
    def arcLenght(self, direction: str, start: float, end: float) -> float:
        raise NotImplementedError

    def FOffset(self, P: np.array, r) -> np.array:
        raise NotImplementedError

    def getMeshData(self):
        return self.verticesGrid, self.triangleMesh

def checkPointInTriangle(point, triangleVertexs):

    def sign (p1, p2, p3):
        return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

    v1, v2, v3 = triangleVertexs

    d1 = sign(point, v1, v2)
    d2 = sign(point, v2, v3)
    d3 = sign(point, v3, v1)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not (has_neg and has_pos)