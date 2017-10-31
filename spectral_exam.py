import numpy as np
import os
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

#from calcos import calcos
from costools import timefilter
import local_lightcurve as lc
from plot_all_lightcurve import gain_change_list_mjd, lpos_list
