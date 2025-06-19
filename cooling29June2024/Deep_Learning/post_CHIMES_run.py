
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import pickle
from numba import njit
import struct
import glob

'''
AbundanceEvolution
TableBins
TemperatureEvolution
TimeArray_seconds
'''


f = h5py.File('grid_noneq_evolution_NeuralNet.hdf5', 'r')

# Print the attributes of HDF5 objects
for name, obj in f.items():
  print(name)
  for key, val in obj.attrs.items():
    print("    %s: %s" % (key, val))
    
TemperatureEvolution = f['TemperatureEvolution'][:]
print(TemperatureEvolution.shape)

AbundanceEvolution = f['AbundanceEvolution']
print('AbundanceEvolution.shape = ', AbundanceEvolution.shape)

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

N_nH = n_densities = f['TableBins/N_Densities'][()]
print("N_Densities:", n_densities)

N_Z = n_metallicities = f['TableBins/N_Metallicities'][()]
print("N_Metallicities:", n_metallicities)

N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
print("N_Temperatures:", n_temperatures)
print()
densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]

print('temperatures = ', temperatures)

inH = 0
print('nH = ', 10**densities[inH])
iTemp = 20
print('T = ', 10**temperatures[iTemp])
iZ = 0
print('Z = ', metallicities[iZ])


TEvol = TemperatureEvolution[iTemp, inH, iZ, :]

nelec = AbundanceEvolution[iTemp, inH, iZ, 0, :] * 10**densities[inH]

nH0 = AbundanceEvolution[iTemp, inH, iZ, 1, :] * 10**densities[inH]
nHp = AbundanceEvolution[iTemp, inH, iZ, 2, :] * 10**densities[inH]

nHe0 = AbundanceEvolution[iTemp, inH, iZ, 4, :] * 10**densities[inH]
nHep = AbundanceEvolution[iTemp, inH, iZ, 5, :] * 10**densities[inH]
nHepp= AbundanceEvolution[iTemp, inH, iZ, 6, :] * 10**densities[inH]


chimesAbundEvol = AbundanceEvolution[iTemp, inH, iZ, :, :] * 10**densities[inH]
print(chimesAbundEvol.shape)

dictZ = {'TempEvol': TEvol, 'chimesAbundEvol': chimesAbundEvol, 't_Arr_in_yrs': t_Arr_in_yrs}
with open('chimesData.pkl', 'wb') as f:
  pickle.dump(dictZ, f)


T = TEvol
nHeTot = nHe0 + nHep + nHepp


plt.figure(figsize = (12, 6))

plt.subplot(1, 2, 1)
plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 5, color = 'k')
plt.xlim(0, 300)
plt.ylim(3, 8)


plt.subplot(1, 2, 2)
plt.scatter(t_Arr_in_yrs, nHe0, s = 5, color = 'orange', label = 'nHe0')
plt.scatter(t_Arr_in_yrs, nHep, s = 5, color = 'lime', label = 'nHep')
plt.scatter(t_Arr_in_yrs, nHepp, s = 5, color = 'yellow', label = 'nHepp')
plt.yscale('log')
plt.legend()

plt.savefig('primordial.png')

plt.show()






