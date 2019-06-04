"""
Created by Ben Kaiser (UNC - Chapel Hill) 2018-06-21

Plot two different COS x1dsum spectra that are each from different objects. I suppose the way this will work is to just scale the spectra using a specific region like the date way for me doing the model spectral comparison (or the ones that don't really have lines, so I guess the way I'm still going to have to do the models for the COS data that only features quasi-molecular satellites... but I digress).

"""


import numpy as np
import sys
import matplotlib.pyplot as plt
import config
import os
from glob import glob
from astropy.io import fits
from astropy import convolution as conv
import config
import spec_plot_tools as spt

target1 = sys.argv[1]
target2 = sys.argv[2]

scaling_range_dict = {"CIV":  [1540, 1560],
                      "Si":  [1380, 1420],
                      "N": [1230, 1260],
                      "CIII": [1170,  1180]}

scaling_range = scaling_range_dict[sys.argv[3]]
#scaling_range = [1700, 1750]
#scaling_range = [1230, 1260] #Nitrogen V
#scaling_range = [1380, 1420] #Si IV
#scaling_range = [1540, 1560] #C IV
#min_wave = 1180
min_wave = 1110
max_wave = 1800

def chi_square(spec1,spec2, err1, err2):
    """
    Should get a chi-square value by comparing one spectrum to the other.
    """
    interp2 = np.interp(spec1[0], spec2[0], spec2[1])
    interp1 = np.interp(spec2[0], spec1[0], spec1[1])
    return

def retrieve_target_spec(target_name):
    """
    Input:
    -target_name : the directory name of the target in HST_COS. It should exclude the slash
    """
    target_path = target_name+'/' 
    target_x1ds= glob(target_path+'*x1dsum.fits')
    target_file= target_x1ds[-1] #I only want the first one for now; when I eventually make this include all the visits for a given target I'll introduce some sort of function for combining them or whatever
    target_hdu= fits.open(target_file)
    cenwave=target_hdu[0].header['CENWAVE']
    grating = target_hdu[0].header['OPT_ELEM']
    seg_gap = config.seg_gap_dict[str(grating)][str(cenwave)]
    lyman_alpha = config.lyman_mask
    oxygen= config.oxygen_mask
    nitrogen = config.nitrogen_mask
    print "grating (OPT_ELEM):", grating, "Centwave: ", target_hdu[0].header['CENWAVE']
    #exposure_time = target_hdu[0].header['EXPTIME']
    target_waves = np.copy(target_hdu[1].data['wavelength'].ravel())
    target_flux= np.copy(target_hdu[1].data['flux'].ravel())
    target_error = np.copy(target_hdu[1].data['ERROR'].ravel())
    target_spec= np.vstack([target_waves, target_flux])
    target_err = np.vstack([target_waves, target_error])
    mask_list= [lyman_alpha]+[oxygen]+[seg_gap]+[nitrogen]
    target_spec = spt.clean_spectrum(target_spec, min_wave, max_wave, mask_list)
    target_err = spt.clean_spectrum(target_err, min_wave, max_wave, mask_list)
    return target_spec, target_err


spec1, error1 = retrieve_target_spec(target1)
spec2, error2 = retrieve_target_spec(target2)

plt.plot(spec2[0], spec2[1], label = target2, color = 'r')
plt.plot(spec1[0], spec1[1], label= target1, color = 'b')
spt.plot_emission_lines(spec1)

plt.title('COS Spectra; NOT RESCALED')
plt.legend()
plt.show()
#conv_kernel = conv.Box1DKernel(2, mode = 'oversample')
#spec2[1] = conv.convolve(spec2[1], conv_kernel)
#spec1[1] = conv.convolve(spec1[1], conv_kernel)

plt.plot(spec2[0], spec2[1]/error2[1], label = target2, color = 'r')
plt.plot(spec1[0], spec1[1]/error1[1], label= target1, color = 'b')
spt.plot_emission_lines(spec1)

plt.title('COS Spectra; NOT RESCALED, SNR')
plt.legend()
plt.ylabel('S/N')
plt.show()


plt.errorbar(spec2[0], spec2[1], yerr=error2[1], label = target2, color = 'r')
plt.errorbar(spec1[0], spec1[1], yerr=error1[1], label= target1, color = 'b')
spt.plot_emission_lines(spec1)

plt.title('COS Spectra; NOT RESCALED')
plt.legend()
plt.show()

#scale_factor= spt.get_scale_factor(spec1, spec2, scaling_range)
scale_factor = spt.get_scale_factor_max(spec1, spec2, scaling_range)

spec2[1] = spec2[1]*scale_factor



plt.plot(spec2[0], spec2[1], label = target2+str(' (rescaled)'), color = 'r')
plt.plot(spec1[0], spec1[1], label= target1, color = 'b')
spt.plot_emission_lines(spec1)

plt.title('COS Spectra; ' + target2 + ' is rescaled to flux values of ' + target1 + ' in ' + str(scaling_range[0]) + '-' + str(scaling_range[1]) + r' $\AA$')
plt.legend()
plt.show()


interp2 = np.interp(spec1[0], spec2[0], spec2[1])
interp1 = np.interp(spec2[0], spec1[0], spec1[1])

plt.plot(spec1[0], interp2)
plt.plot(spec1[0], spec1[1])
plt.show()

plt.plot(spec2[0], interp1)
plt.plot(spec2[0], spec2[1])
plt.show()


#plt.plot(spec1[0], spec1[1], label= target1, color = 'b')
plt.plot(spec1[0], spec1[1]/interp2, label=target1 + '/'  +target2)
spt.plot_emission_lines(spec1)

#plt.title('COS Spectra; ' + target1 + '/ ' + target2 + ' in ' + str(scaling_range[0]) + '-' + str(scaling_range[1]) + r' $\AA$')
#plt.legend()
#plt.show()


plt.plot(spec2[0], spec2[1]/interp1, label = target2 + '/' + target1)
spt.plot_emission_lines(spec2)
plt.axhline(y=1, color = 'k', alpha = 0.5)

#plt.title('COS Spectra; ' + target2 + '/ ' + target1 + ' in ' + str(scaling_range[0]) + '-' + str(scaling_range[1]) + r' $\AA$')
plt.legend()
plt.show()
