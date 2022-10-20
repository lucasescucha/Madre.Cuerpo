import numpy as np


class ISurface:
    def F(self, P: np.array) -> float:
        pass

    def gradF(self, P: np.array) -> np.array:
        pass
    
    def arcLenght(self, direction: str, start: float, end: float) -> float:
        pass

    def FOffset(self, P: np.array, r: float) -> np.array:
        pass
