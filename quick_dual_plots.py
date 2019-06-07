from __future__ import print_function


import numpy as np
import sys
import run_dual_plot_grid as rdpg

#target_directories = np.genfromtxt('target_dirs.txt', 'str')
target_directories = np.genfromtxt('target_for_lc.txt', 'str')


stepsize = int(sys.argv[1])
for target in target_directories:
    target= str(target)
    if 'grid_lightcurves' in target:
        pass
    else:
        try:
            rdpg.run_dual_plot_grid(target[:-1], stepsize) #removes the slash at the end of the directory name
        except IOError as error:
            print(error)
