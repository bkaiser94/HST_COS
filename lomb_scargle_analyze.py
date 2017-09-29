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
num_freq = 10000

all_array = np.genfromtxt(inputfile, names=True)
#times= Time(all_array['mjd'], format='mjd')
times= np.copy(all_array['mjd'])
#times = (times- times[0])*24.*3600. #rezeroing the time of observation.
times = (times- times[0]) #rezeroing the time of observation.

#times= times * u.day
#times= times.to(u.s)
#print times
#period_range = np.array([10,20]) #hours
#period_range= [100, 1400] #seconds
period_range= [1/400., 1/1e-8]
#period_range= period_range *3600. #seconds
frequency_range= np.linspace(1./period_range[1],1./period_range[0], num_freq)
frequency = frequency_range
print "max freq: " , np.nanmax(frequency)
print "freq min: ", np.nanmin(frequency)
print "freq step: ", frequency[1]-frequency[0]
#frequency, power = LombScargle(times, all_array['flux']).autopower()
#power = LombScargle(times, all_array['flux']).power(frequency_range)
power = LombScargle(times, all_array['flux']-np.nanmean(all_array['flux'])).power(frequency_range)

clean_indices= ~np.isnan(power)
clean_power= power[clean_indices]
clean_frequency= frequency[clean_indices]

best_index= np.argmax(clean_power)
best_power= clean_power[best_index]
best_freq= clean_frequency[best_index]

print "best_freq: ", best_freq, best_freq/86400.
print "period: ", 1./best_freq*86400., "s"
print "best_power: ", best_power



plt.figure(figsize= (20,9))
plt.plot(frequency, power)
#plt.xlabel('frequency (Hz)')
plt.xlabel('frequency')
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
