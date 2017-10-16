import numpy as np
from glob import glob
import os
import sys

import try_lightcurve as tlc


"""
calls try_lightcurve.py and then runs it over a grid of step sizes, and potentially will iterate over the different target directories.

"""

target_dir = sys.argv[1] + '/'
step_grid= np.arange(1,121, 1) #by going one above the desired value, it is included in the range. For example to include 120 make it 121.
limit_list= [[1130,1800],[1130,1850],[1130,1900]]

def run_step_grid(target_dir):
    for wave_limit in limit_list:
        for stepsize in step_grid:
            outstring = tlc.make_lightcurve(target_dir, stepsize, wave_limit, plotall = False)
    return ''

if __name__ == '__main__':
    run_step_grid(target_dir)

