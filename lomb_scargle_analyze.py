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
num_freq = 1000

all_array = np.genfromtxt(inputfile, names=True)
#times= Time(all_array['mjd'], format='mjd')
times= np.copy(all_array['mjd'])
times = (times- times[0])*24.*3600. #rezeroing the time of observation.
#times= times * u.day
#times= times.to(u.s)
#print times
#period_range = np.array([10,20]) #hours
period_range= [100, 1400] #seconds
#period_range= period_range *3600. #seconds
frequency_range= np.linspace(1./period_range[1],1./period_range[0], 100000)
frequency = frequency_range
print "max freq: " , np.nanmax(frequency)
print "freq min: ", np.nanmin(frequency)
#frequency, power = LombScargle(times, all_array['flux']).autopower()
power = LombScargle(times, all_array['flux']).power(frequency_range)


plt.figure(figsize= (20,9))
plt.plot(frequency, power)
plt.xlabel('frequency (Hz)')
plt.title(inputfile + ' LombScargle')
plt.show()



plt.figure(figsize= (20,9))
plt.plot(1./frequency, power)
plt.xlabel("Period (seconds)")
plt.title(inputfile + ' LombScargle')
plt.show()


plt.figure(figsize= (20,9))
plt.plot(np.log10(1./frequency), power)
plt.xlabel("log10(Period (seconds))")
plt.title(inputfile + ' LombScargle')
plt.show()


#plt.figure(figsize= (20,9))
#plt.plot((1./frequency)/3600., power)
#plt.xlabel("Period (hrs)")
#plt.title(inputfile + ' LombScargle')
#plt.show()
