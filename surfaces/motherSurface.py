from surfaces.iSurface import ISurface
import numpy as np


class MotherSurface(ISurface):
    def __init__(self, a: float, b: float) -> None:
        self.a = a
        self.b = b

    def F(self, P: np.array) -> float:
        x, y = P
        return self.a*(x ** 2) - self.b*(y ** 2)

    def gradF(self, P: np.array) -> np.array:
        x, y = P
        return np.array([2*self.a*x, -2*self.b*y])
