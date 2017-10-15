import numpy as np
import os
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

from calcos import calcos
from costools import timefilter
import lightcurve as lc

target_dir= sys.argv[1]
#lcfile= sys.argv[2]
stepsize = int(sys.argv[2])
period= float(sys.argv[3])
unit_arg = sys.argv[4]

lcbase= target_dir+ '_grid_lightcurves/' + '*step' + str(stepsize) + '_*'
lcfile= glob(lcbase)[0]
dest_dir = 'dual_plots/'
if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
wave_min= 1130
lyman = [1208, 1225]
oxygen= [1295, 1312] #airglow wavelengths to be filtered according to lightcurve
wave_max= 1900

for dataset in glob(target_dir+ '/*x1d.fits'):
    hdu = fits.open(dataset)
    print dataset, hdu[0].header['OPT_ELEM'], hdu[0].header['CENWAVE'], hdu[1].header['EXPTIME']    

fig = plt.figure(figsize=(20,9))
ax1 = fig.add_subplot(2,1,1)

def remove_range(wave_array, other_array, bound_list):
    """
    Removes wavelengths and flux values from the array that fall in the range specified by bound_list
    """
    lower_bound = bound_list[0]
    upper_bound= bound_list[1]
    low_mask = np.where(wave_array < lower_bound)
    high_mask= np.where(wave_array > upper_bound)
    low_waves= wave_array[low_mask]
    high_waves= wave_array[high_mask]
    low_other = other_array[low_mask]
    high_other= other_array[high_mask]
    merge_waves= np.append(low_waves, high_waves)
    merge_other= np.append(low_other, high_other)
    return merge_waves, merge_other

counter= 0
for dataset in glob(target_dir+ '/*x1d.fits'):
    hdu = fits.open(dataset)
    #print hdu[1]
    #print hdu[1]['FUVA']
    #print hdu[1].data[0]
    #print hdu[1].data.shape
    #for thing in hdu[1].data:
        #print "thing: ", thing
        #print thing['wavelength']
    wavelengths= np.copy(hdu[1].data['wavelength'].ravel())
    print wavelengths.shape
    fluxes= np.copy(hdu[1].data['flux'].ravel())
    #print "fluxes",  fluxes
    arrshape= wavelengths.shape
    #unmasked= np.where((wavelengths < lyman[0]*np.ones(arrshape)) or ((wavelengths > lyman[1]*np.ones(arrshape)) and (wavelengths < oxygen[0])*np.ones(arrshape)) or (wavelengths > oxygen[1]*np.ones(arrshape)))
    lower_mask = np.where(wavelengths > wave_min)
    wavelengths= np.copy(wavelengths[lower_mask])
    fluxes = np.copy(fluxes[lower_mask])
    upper_mask= np.where(wavelengths < wave_max)
    wavelengths= np.copy(wavelengths[upper_mask])
    fluxes= np.copy(fluxes[upper_mask])
    wavelengths, fluxes = remove_range(wavelengths, fluxes, lyman)
    wavelengths, fluxes = remove_range(wavelengths, fluxes, oxygen)
    #unmasked= np.where(wavelengths < lyman[0])
    #wavelengths1= wavelengths[unmasked]
    #print wavelengths1.shape
    #print "fluxes: ",fluxes
    #fluxes1= fluxes[unmasked]
    #newmask= np.where(wavelengths > lyman[1])
    #wavelengths2= wavelengths[newmask]
    #print wavelengths2.shape
    #fluxes2= fluxes[newmask]
    #midmask= np.where(wavelengths2 < oxygen[0])
    #wavelengths2= wavelengths2[midmask]
    #print wavelengths2.shape
    #fluxes2= fluxes2[midmask]
    #lastmask= np.where(wavelengths > oxygen[1])
    #wavelengths3= wavelengths[lastmask]
    #print wavelengths3.shape
    #fluxes3= fluxes[lastmask]
    #lastermask= np.where(wavelengths3 < wave_max)
    #wavelengths3= wavelengths3[lastermask]
    #fluxes3= fluxes3[lastermask]
    ##now recombine all of those masked shenanigans
    #wavelengths= np.append(wavelengths1, wavelengths2)
    #wavelengths= np.append(wavelengths, wavelengths3)
    #fluxes= np.append(fluxes1, fluxes2)
    #fluxes = np.append(fluxes, fluxes3)
    print "fluxes.shape: ", fluxes.shape
    print "wavelengths.shape" , wavelengths.shape
    if counter == 0:
        flux_all= np.zeros(fluxes.shape)
        plot_waves = wavelengths
    try:
        flux_all = np.copy(fluxes+flux_all)
    except ValueError as error:
        print error
        flux_all= np.copy(flux_all +np.append(fluxes, 0.))
    counter += 1
    print "sum fluxes: ", np.sum(fluxes)
    #ax1.plot(wavelengths, fluxes/np.nanmean(fluxes), label=hdu[0].header['rootname'])
ax1.plot(plot_waves, flux_all/np.nanmean(flux_all), label= 'summed spectra')
    
    
    
ax1.legend(numpoints=1, fontsize=14, loc='best' )
    



ax1.set_ylabel('Flux (normed)')
ax1.set_xlabel('Wavelength $(\AA)$')
ax1.set_title(target_dir)
#ax1.set_yscale('log')

ax2= fig.add_subplot(2,1,2)

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
    

all_array = np.genfromtxt(lcfile, names=True)
#times= Time(all_array['mjd'], format='mjd')
times= np.copy(all_array['mjd'])
fluxes= np.copy(all_array['flux'])

times= (times - times[0]) *time_converter
fold_times = times%period


#plt.figure(figsize= (20,9))
#ax2.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
ax2.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
ax2.scatter(fold_times, fluxes)
ax2.set_xlabel("Time ("+ time_string + ")")
ax2.set_ylabel("Flux (normed and zeroed)")
ax2.set_xlim(0, period)
ax2.set_title(lcfile + ' Period fold '+ str(period) + ' ' + time_string)
#plt.show()








fig.tight_layout()
#plt.show()

fig.savefig(dest_dir+ target_dir + '_dual_plot_fold_'+ sys.argv[3] + sys.argv[4]+'_step' + str(stepsize)+'.pdf', bbox_inches = 'tight')
