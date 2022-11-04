from core.iFunction import IFunction
import numpy as np
import scipy.optimize as opt

SOLVER_TOLERANCE = 1e-5

def solveF_Factory(Po, gradF, r):
    def solveF(P):
        n_hat = getNormalHat(gradF(P))
        remainder = Po - P - (n_hat[0:-1] * r)
        
        return remainder if np.linalg.norm(remainder) > SOLVER_TOLERANCE \
                    else np.zeros_like(remainder)

    return solveF

def getNormalVector(gradF_P):
    return np.append(-gradF_P, [1])

def getNormalHat(gradF_P):
    normal = getNormalVector(gradF_P)
    return normal/np.linalg.norm(normal)

def getOffsetFunctionZ(function: IFunction, Po, r):
    solveF = solveF_Factory(Po, function.gradF, r)
    Pr = opt.excitingmixing(solveF, Po, f_tol=SOLVER_TOLERANCE)
    n_hat = getNormalHat(function.gradF(Pr))
    return function.F(Pr) + (n_hat[-1] * r)