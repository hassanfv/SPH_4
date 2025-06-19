import numpy as np


SAFETY = 0.9
PGROW = -0.2
PSHRNK = -0.25
ERRCON = 1.89e-4

h = 0.1
errmax = 0.09

htemp = SAFETY * h * errmax**PSHRNK

print(htemp)
