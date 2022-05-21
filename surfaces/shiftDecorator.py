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
