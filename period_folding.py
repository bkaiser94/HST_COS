import numpy as np
import os
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.time import Time
from astropy.stats import LombScargle
import sys

from calcos import calcos
from costools import timefilter
import lightcurve as lc


inputfile= sys.argv[1]
period= float(sys.argv[2])
unit_arg = sys.argv[3]
second_per_mjd= 1./1.15741e-5     #SECOND_PER_MJD value from lightcurve.cos.extract, but I couldn't import it for whatever reason... so I just copied and pasted

if unit_arg.startswith('s'):
    print "period in seconds"
    time_converter = second_per_mjd
    time_string = 's'
    
elif unit_arg.startswith('h'):
    print "period in hours"
    time_converter= second_per_mjd/3600.
    time_string = 'hrs'

elif unit_arg.startswith('d'):
    print 'period in days'
    time_converter = 1.
    time_string = 'days'

elif unit_arg.startswith('m'):
    print 'period in minutes'
    time_converter = second_per_mjd/60.
    time_string = 'min'
    

all_array = np.genfromtxt(inputfile, names=True)
#times= Time(all_array['mjd'], format='mjd')
times= np.copy(all_array['mjd'])
fluxes= np.copy(all_array['flux'])

times= (times - times[0]) *time_converter
fold_times = times%period


plt.figure(figsize= (20,9))
#plt.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
plt.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
plt.scatter(fold_times, fluxes)
plt.xlabel("Time ("+ time_string + ")")
plt.xlim(0, period)
plt.title(inputfile + ' Period fold'+ str(period) + ' ' + time_string)
plt.show()


#iterative_periods= np.linspace(10.,20,20)
#for ped in iterative_periods:
    ##plt.figure(figsize= (20,9))
    #fold_times= times%ped
    #plt.scatter(fold_times, fluxes)
    #plt.xlabel("Time ("+ time_string + ")")
    #plt.title(inputfile + ' Period fold'+ str(ped) + ' ' + time_string)
    #plt.show()
