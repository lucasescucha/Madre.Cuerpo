import math
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
    
    #\frac{1}{2}x\sqrt{4b^2x^2+1}+\frac{1}{4b}\ln \left|2bx+\sqrt{1+4b^2x^2}\right|+C
    def arcLenght(self, direction: str, start: float, end: float) -> float:
        def integral(k, v):
            return 0.5*v*math.sqrt(4*(k**2)*(v**2) + 1)+(1/(4*k))*math.log(math.abs(2*k*v + math.sqrt(1+4*(k**2)*(v**2))))

        k = self.a if direction == "x" else self.b 
        
        return integral(k, end) - integral(k, start)  
