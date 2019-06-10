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

from calcos import calcos
from costools import timefilter
import lightcurve as lc


inputfile= sys.argv[1]
period= np.float_(sys.argv[2])
unit_arg = sys.argv[3]
second_per_mjd= 1./1.15741e-5     #SECOND_PER_MJD value from lightcurve.cos.extract, but I couldn't import it for whatever reason... so I just copied and pasted
phase_bins = 120. #number of bins to do with in-phase binning

if unit_arg.startswith('s'):
    print("period in seconds")
    time_converter = second_per_mjd
    time_string = 's'
    
elif unit_arg.startswith('h'):
    print("period in hours")
    time_converter= second_per_mjd/3600.
    time_string = 'hrs'

elif unit_arg.startswith('d'):
    print('period in days')
    time_converter = 1.
    time_string = 'days'

elif unit_arg.startswith('m'):
    print('period in minutes')
    time_converter = second_per_mjd/60.
    time_string = 'min'
    
def bin_in_phase(fold_times, flux):
    #sorted_indices = np.argsort(fold_times)
    #sorted_fold_times = fold_times[sorted_indices]
    #sorted_flux = flux[sorted_indices]
    bin_width = period/phase_bins
    bin_edges = np.linspace(0,period+bin_width, phase_bins)
    binned_times = []
    binned_flux = []
    binned_std = []
    for i in range(0,int(phase_bins)-1):
        condition  = (fold_times >= bin_edges[i]).astype(int) * (fold_times <= bin_edges[i+1]).astype(int) #array of zeros and ones for values that are in bounds
        good_vals = np.where(condition> 0)
        good_flux = flux[good_vals]
        print(good_flux.shape[0])
        std_err = np.std(good_flux)
        bin_flux = np.mean(good_flux)
        binned_flux.append(bin_flux)
        binned_std.append(std_err)
        binned_times.append(np.mean([bin_edges[i], bin_edges[i+1]]))
    return np.array(binned_times), np.array(binned_flux), np.array(binned_std)
        
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
#times= np.copy(all_array['mjd'])
#times= np.copy(all_array['bmjd_tdb'])
fluxes= np.copy(all_array['flux'])

times= (times - times[0]) *time_converter
#fold_times = times%period
fold_times= np.mod(times, period)

print("Base time unit:", times[1]-times[0], times[2]-times[1], times[3]-times[2])


plt.figure(figsize= (20,9))
#plt.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
plt.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
plt.scatter(fold_times, fluxes)
plt.xlabel("Time ("+ time_string + ")")
plt.xlim(0, period)
plt.title(inputfile + ' Period fold'+ str(period) + ' ' + time_string)
plt.show()

binned_times, binned_flux, binned_std = bin_in_phase(fold_times, fluxes)
print(binned_std)
plt.figure(figsize= (20,9))
#plt.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
plt.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
#plt.scatter(binned_times, binned_flux)
#plt.plot(binned_times, binned_flux, marker = 'o')
plt.errorbar(binned_times, binned_flux, yerr= binned_std,  marker = 'o')

plt.xlabel("Time ("+ time_string + ")")
plt.ylabel("Mean normalized flux")
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
