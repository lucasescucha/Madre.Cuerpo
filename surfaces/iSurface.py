from core.iFunction import IFunction
import numpy as np


class ISurface(IFunction):
    def arcLenght(self, direction: str, start: float, end: float) -> float:
        pass

    def FOffset(self, P: np.array, r: float) -> np.array:
        pass
