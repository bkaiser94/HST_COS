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

second_per_mjd= 1./1.15741e-5     #SECOND_PER_MJD value from lightcurve.cos.extract, but I couldn't import it for whatever reason... so I just copied and pasted
inputfile= sys.argv[1]
num_freq = 100000
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
times = (times- times[0])*second_per_mjd #rezeroing the time of observation.
#times = (times- times[0]) #rezeroing the time of observation.

#times= times * u.day
#times= times.to(u.s)
#print(times)
#period_range = np.array([10,20]) #hours
#period_range= [100, 1400] #seconds
#period_range= [70., 1400.]
#period_range= [2*second_per_mjd, 4000*second_per_mjd]
period_range=[2,2000]
#period_range= [ 5.00001,2000.]
#period_range = [70.,190.]
#period_range= period_range *3600. #seconds
frequency_range= np.linspace(1./period_range[1],1./period_range[0], num_freq)
frequency = frequency_range
print("max freq: " , np.nanmax(frequency))
print("freq min: ", np.nanmin(frequency))
print("freq step: ", frequency[1]-frequency[0])
#frequency, power = LombScargle(times, all_array['flux']).autopower()
#power = LombScargle(times, all_array['flux']).power(frequency_range)
power = LombScargle(times, all_array['flux']).power(frequency)

clean_indices= ~np.isnan(power)
clean_power= power[clean_indices]
clean_frequency= frequency[clean_indices]

best_index= np.argmax(clean_power)
best_power= clean_power[best_index]
best_freq= clean_frequency[best_index]

#print("best_freq: ", best_freq, best_freq/86400.)
print("best_freq: ", best_freq, "Hz")
#print("period: ", 1./best_freq*86400., "s")
print("Best period: ", 1./best_freq, "s")
print("best_power: ", best_power)

sorted_indices= np.argsort(clean_power)
sorted_power= clean_power[sorted_indices]
sorted_freq= clean_frequency[sorted_indices]

print("top 10 frequencies: ")
counter= 1
for top_freq in reversed(sorted_freq[-10:]):
    print(top_freq, "Hz\tPeriod: ", 1./top_freq, "s\t""\t=", 1./top_freq/second_per_mjd, "d\t", "Power: ", sorted_power[counter*-1])
    counter+=1


plt.figure(figsize= (20,9))
plt.plot(frequency, power)
plt.xlabel('frequency (Hz)')
#plt.xlabel('frequency')
#plt.ylim([0,0.02])
plt.yscale('log')
plt.title(inputfile + ' LombScargle')
plt.show()


plt.figure(figsize= (20,9))
plt.plot(frequency, power)
plt.xlabel('frequency (Hz)')
#plt.xlabel('frequency')
#plt.ylim([0,0.02])
plt.title(inputfile + ' LombScargle')
plt.show()


#plt.figure(figsize= (20,9))
#plt.plot(1./frequency, power)
#plt.xlabel("Period (seconds)")
#plt.title(inputfile + ' LombScargle')
#plt.show()


#plt.figure(figsize= (20,9))
#plt.plot(np.log10(1./frequency), power)
#plt.xlabel("log10(Period (seconds))")
#plt.title(inputfile + ' LombScargle')
#plt.show()


#plt.figure(figsize= (20,9))
#plt.plot((1./frequency)/3600., power)
#plt.xlabel("Period (hrs)")
#plt.title(inputfile + ' LombScargle')
#plt.show()
