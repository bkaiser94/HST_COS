import numpy as np
import sys
import matplotlib.pyplot as plt
import config
import os
from glob import glob
from astropy.io import fits

"""
In the future I'll have to change the model fitting to accomodate the float or integer input values for logg and teff, but for tonight, I'm just going to input them as correctly padded strings to save myself a headache for the moment
"""

target_path = sys.argv[1]+'/' #the input should be the target name just like all of the other times this has been implemented in this repository
#path to the model files
teff= "013500" #hard coded for test target
#logg= 8.0 #model values to pull (should be changed to be read from somewhere in the future)
logg = "800" #hard coded for test target(it was 7.8 in the paper, but that is technically the same as this gridpoint to two sig figs)
lyman_alpha = [1214,1217] #probably should put this into the config file as the masking for the spectroscopic fitting (but really this goes back to the need for individual text files by target for the various maskings that are required)
oxygen= [1300, 1308]
seg_gap = [1268,1298] #segment gap mask used for photometry for the G130M at 1291 Angs
scale_wave_range= [1400, 1420] 
#scale factor to be just typed in for the moment but should be replaced by a function at some point that calculates it for you
target_x1ds= glob(target_path+'*x1dsum.fits')
target_file= target_x1ds[0] #I only want the first one for now; when I eventually make this include all the visits for a given target I'll introduce some sort of function for combining them or whatever
model_path= '/Users/BenKaiser/Desktop/DK_2010/'
#scaling_coefficient= 5.2e21
#wavelength range to consider for the scaling

target_hdu= fits.open(target_file)
print "Centwave: ", target_hdu[0].header['CENWAVE']
target_waves = np.copy(target_hdu[1].data['wavelength'].ravel())
target_flux= np.copy(target_hdu[1].data['flux'].ravel())
target_spec= np.vstack([target_waves, target_flux])
print target_spec.shape
print target_spec[0]


#function to trim the model plot to limit it to the range of the target's spectrum
#At some point I'll have to implement masking

def get_model_spec(fitted_teff, fitted_logg):
    """
    inputs: fitted_teff, fitted_logg
    
    output: 2-d numpy array containing the wavelengths and spectral values from the model grid for the input teff and logg
    """
    padded_teff = fitted_teff #needs to be changed to pad/remove decimals (and convert to string)
    padded_logg= fitted_logg #needs to be changed to pad/remove decimals (and convert to string)
    model_file= 'da' +padded_teff + '_' + padded_logg + '.dat'
    model_file= model_path + model_file
    model_file= glob(model_file)[0]
    model_spec= np.genfromtxt(model_file)
    model_spec= model_spec.T #transpose because numpy is weird on the axes it reads in.
    print "model_spec.shape:",  model_spec.shape
    return model_spec

def trim_spec(input_spec, min_wave, max_wave):
    lower_indices = np.where(input_spec[0]< max_wave)
    trimmed_waves= input_spec[0][lower_indices]
    trimmed_flux= input_spec[1][lower_indices]
    upper_indices= np.where(trimmed_waves > min_wave)
    trimmed_waves= trimmed_waves[upper_indices]
    trimmed_flux = trimmed_flux[upper_indices]
    trimmed_spec= np.vstack([trimmed_waves, trimmed_flux])
    print trimmed_spec.shape
    return trimmed_spec

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
    return np.vstack([merge_waves, merge_other])
#########

def get_scale_factor(target_spec, model_spec, wave_range):
    allowed_target = np.where((target_spec[0] > wave_range[0])& (target_spec[0] < wave_range[1]))
    allowed_model = np.where((model_spec[0] > wave_range[0]) & (model_spec[0] < wave_range[1]))
    target_mean= np.mean(target_spec[1][allowed_target])
    model_mean = np.mean(model_spec[1][allowed_model])
    scaling_factor =target_mean/model_mean
    return scaling_factor

###
def plot_overlays(spec1, spec2, scaling_coefficient):
    plt.plot(spec1[0], spec1[1], label = 'observed')
    plt.plot(spec2[0], spec2[1]*scaling_coefficient, label= 'model', color = 'r')
    plt.legend(numpoints=1, fontsize=14, loc='best' )
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (cgs units)')
    plt.title(target_path[:-1] +'_' +  teff + '_' + logg)
    plt.show()
    return ''

def sort_spectrum(input_spec):
    """
    Hopefully this fixes the apparent errant lines throughout the line plots of target spectra.
    """
    sort_indices = np.argsort(input_spec[0])
    sorted_waves= input_spec[0][sort_indices]
    sorted_flux= input_spec[1][sort_indices]
    sorted_spectrum= np.vstack([sorted_waves, sorted_flux])
    
    return sorted_spectrum
##########
print os.getcwd()

model_spec= get_model_spec(teff, logg)
model_spec= trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
target_spec= remove_range(target_spec[0], target_spec[1], lyman_alpha)
target_spec= remove_range(target_spec[0], target_spec[1], oxygen)
target_spec= remove_range(target_spec[0], target_spec[1], seg_gap)
target_spec = sort_spectrum(target_spec)
scaling_coefficient= get_scale_factor(target_spec, model_spec, scale_wave_range)
plot_overlays(target_spec, model_spec, scaling_coefficient)
