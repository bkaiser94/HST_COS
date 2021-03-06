from __future__ import print_function


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
import config
#from calcos import calcos
#from costools import timefilter
import local_lightcurve as lc



#lpos_list= [56132, 57063] #mjd dates of position changes
#gain_change_list=Time( ['2009-05-11', '2009-08-12', '2011-03-08', '2012-03-26', '2012-07-23', '2013-06-24', '2014-07-21', '2014-11-03', '2015-02-09'], scale='utc')
gain_change_list_mjd = Time(config.gain_change_list, scale='utc').mjd
print(gain_change_list_mjd)
second_per_mjd= 1./1.15741e-5 
def plot_events(times):
    min_time= np.min(times)
    max_time = np.max(times)
    #for gain_change in gain_change_list_mjd:
        ##if ((fppos > times.min) & (fppos < times.max)):
        #if ((gain_change > min_time) & (gain_change< max_time)):
            #plt.axvline(x= gain_change ,linestyle = '-', color = 'r' , ymin = -10, ymax = 10, linewidth = 1, alpha = 0.2)
    #for lpos in config.lpos_list:
        ##if ((fppos > times.min) & (fppos < times.max)):
        #if ((lpos > min_time) & (lpos < max_time)):
            #plt.axvline(x= lpos,linestyle = '-', color = 'k' , ymin = -10, ymax = 10, linewidth = 1, alpha = 0.8)
    gain_change_list_bmjd = Time(config.gain_change_list, scale= 'utc')
    gain_change_list_bmjd= gain_change_list_bmjd.tdb.mjd
    for gain_change in gain_change_list_bmjd:
        #if ((fppos > times.min) & (fppos < times.max)):
        if ((gain_change > min_time) & (gain_change< max_time)):
            plt.axvline(x= gain_change ,linestyle = '-', color = 'r' , ymin = -10, ymax = 10, linewidth = 1, alpha = 0.2)
    for lpos in Time(config.lpos_list, format = 'mjd', scale = 'utc').tdb.mjd:
        #if ((fppos > times.min) & (fppos < times.max)):
        if ((lpos > min_time) & (lpos < max_time)):
            plt.axvline(x= lpos,linestyle = '-', color = 'k' , ymin = -10, ymax = 10, linewidth = 1, alpha = 0.8)

def plot_stuff(inputfile):
    try:
        all_array = np.genfromtxt(inputfile, names=True, skip_header=1)
        times= np.copy(all_array['bmjd_tdb'])
    except ValueError as error:
        print("ValueError:",error)
        print("probably due to this lightcurve being before the addition of comments,\nso now we won't skip any lines.")
        all_array = np.genfromtxt(inputfile, names=True)
        times= np.copy(all_array['bmjd_tdb'])
        
    #all_array = np.genfromtxt(inputfile, names=True)
    #times= Time(all_array['mjd'], format='mjd')
    #times= np.copy(all_array['bmjd_tdb'])
    fluxes= np.copy(all_array['flux'])
    #fluxes= np.copy(all_array['counts'])
    try:
        gross = np.copy(all_array['gross']) #approximate poisson noise of the plot.
        flux_err = np.sqrt(gross)/gross
        #flux_err = np.sqrt(gross)
    except ValueError as error:
        print(error)
        print("Flux uncertainties wrong")
        flux_err= 0
    stepsize = inputfile.split('step')[-1]
    stepsize= float(stepsize.split('_')[0]) #hopefully the stepsize
    time_err = stepsize/2./second_per_mjd


    plt.figure(figsize= (20,9))
    #plt.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
    plt.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
    plot_events(times)
    #plt.scatter(times, fluxes)
    plt.errorbar(times, fluxes, flux_err, time_err, fmt= 'o')
    plt.xlabel("Time (BMJD_TDB)")
    plt.ylabel("Flux Residuals about the Mean")
    plt.title(inputfile)
    plt.show()

if __name__ == '__main__':
    inputfile= sys.argv[1]
    plot_stuff(inputfile)
