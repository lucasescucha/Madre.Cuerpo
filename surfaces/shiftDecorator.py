from surfaces.iSurface import ISurface
import numpy as np


class ShiftDecorator(ISurface):

    _surface: ISurface = None
    _shift: np.array = None

    def __init__(self, surface: ISurface, shift: np.array) -> None:
        self._surface = surface
        self._shift = shift

    @property
    def surface(self) -> str:
        return self._surface

    def F(self, P: np.array) -> float:
        shift = self._shift
        return self._surface.F(P + shift[0:2]) + shift[2]

    def gradF(self, P: np.array) -> np.array:
        shift = self._shift
        return self._surface.gradF(P + shift[0:2])

    def arcLenght(self, direction: str, start: float, end: float) -> float:
        shift = self._shift
        offset =  shift[0] if direction == "x" else shift[1]
        return self._surface.arcLenght(direction, start + offset, end + offset)
    
    def FOffset(self, P: np.array, r: float) -> np.array:
        normal = np.append((self.gradF(P) * -1), [1])
        hat = normal/np.linalg.norm(normal)
        return np.append(P, [self.F(P)]) + hat*r
