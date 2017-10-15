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


for stepsize in step_grid:
    outstring = tlc.make_lightcurve(target_dir, stepsize, plotall = False)
