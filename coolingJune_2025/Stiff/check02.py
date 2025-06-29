import numpy as np

def f1(t):
    return -2. * np.exp(-t) + 1000. * np.exp(-1000.*t)


#print('Result = ', f1(1e-3))

h = 0.0001

t = 0.0
u = 1.0 # initial value of u at t = 0 is 1.0 ! 

for i in range(1000):

    t += h
    u = u + h * f1(t)
    if t == 0.01:
        print(t, u, f1(t))

print()
#print(t, u)