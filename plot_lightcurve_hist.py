import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import sys

inputfile = sys.argv[1]
try:
    bin_num= int(sys.argv[2])
except IndexError:
    bin_num = 100
all_array= np.genfromtxt(inputfile, names= True)

fluxes= np.copy(all_array['flux'])
plt.figure(figsize= (20,9))
plt.hist(fluxes, bins= bin_num)
plt.axvline(x= 0,linestyle = '-', color = 'k' , ymin = 0, ymax = 100000, linewidth = 1, alpha = 0.2)
plt.xlabel("Residual Flux Values")
plt.ylabel("N Bins with Flux Value")
plt.title(inputfile)
plt.show()
