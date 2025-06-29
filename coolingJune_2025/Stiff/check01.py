import numpy as np


def f1(t):
    return -2. * np.exp(-t) + 1000. * np.exp(-1000.*t)


print('Result = ', f1(1e-3))

