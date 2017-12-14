import numpy as np
from glob import glob
import os
import sys
import config
#from run_step_grid import limit_list
from make_dual_plots import make_dual_plots

def run_dual_plot_grid(target_dir, stepsize):
    for limit in config.wave_limit_list:
        #nullstring = make_dual_plots(target_dir, stepsize, period, unit_arg, wave_limits= limit)
        try:
            nullstring = make_dual_plots(target_dir, stepsize,  wave_limits= limit)
        except IndexError as error:
            print "Wavelengths: ", limit
            print error



if __name__ == "__main__":
    #import try_lightcurve as tlc
    target_dir= sys.argv[1]
    #lcfile= sys.argv[2]
    stepsize = int(sys.argv[2])
    #period= float(sys.argv[3])
    #unit_arg = sys.argv[4]
    run_dual_plot_grid(target_dir, stepsize)
