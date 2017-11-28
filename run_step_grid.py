import numpy as np
from glob import glob
import os
import sys
import config

import try_lightcurve as tlc


"""
calls try_lightcurve.py and then runs it over a grid of step sizes, and potentially will iterate over the different target directories.

"""

target_dir = sys.argv[1] + '/'
#step_grid= np.arange(1,121, 1) #by going one above the desired value, it is included in the range. For example to include 120 make it 121.
#step_grid= np.arange(15,121, 15) #by going one above the desired value, it is included in the range. For example to include 120 make it 121.
#limit_list= [[1130,1430],[1130,1800],[1130,1850],[1130,1900]]
step_grid= np.arange(5,31,5)
app_array= np.arange(45,61,15)
step_grid= np.append(step_grid, app_array)
step_grid= np.append(step_grid, [90,120,180,240,300])

def run_step_grid(target_dir):
    for wave_limit in config.wave_limit_list:
        outstring = tlc.make_lightcurve(target_dir, 1, wave_limit, plotall = False) #one second binning before running the other steps in the grid
        for stepsize in step_grid:
            outstring = tlc.make_lightcurve(target_dir, stepsize, wave_limit, plotall = False)
    return ''

if __name__ == '__main__':
    run_step_grid(target_dir)

