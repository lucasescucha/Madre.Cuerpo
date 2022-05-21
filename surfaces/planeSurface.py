from surfaces.iSurface import ISurface
import numpy as np


class PlaneSurface(ISurface):
    def __init__(self, a: float, b: float):
        self.a = a
        self.b = b

    def F(self, P: np.array) -> float:
        x, y = P
        return self.a*x + self.b*y

    def gradF(self, P: np.array) -> np.array:
        return np.array([self.a, self.b])
