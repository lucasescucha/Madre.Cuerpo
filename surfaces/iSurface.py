import numpy as np


class ISurface:
    def F(self, P: np.array) -> float:
        pass

    def gradF(self, P: np.array) -> np.array:
        pass
