import numpy as np
from glob import glob
import os
import sys
from run_step_grid import limit_list
from make_dual_plots import make_dual_plots

#import try_lightcurve as tlc
target_dir= sys.argv[1]
#lcfile= sys.argv[2]
stepsize = int(sys.argv[2])
#period= float(sys.argv[3])
#unit_arg = sys.argv[4]

for limit in limit_list:
    #nullstring = make_dual_plots(target_dir, stepsize, period, unit_arg, wave_limits= limit)
    nullstring = make_dual_plots(target_dir, stepsize,  wave_limits= limit)



