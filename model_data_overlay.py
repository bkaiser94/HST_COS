import numpy as np
import sys
import matplotlib.pyplot as plt
import config
import os
from glob import glob
from astropy.io import fits

import config
import spec_plot_tools as spt

"""
In the future I'll have to change the model fitting to accomodate the float or integer input values for logg and teff, but for tonight, I'm just going to input them as correctly padded strings to save myself a headache for the moment
"""

target_path = sys.argv[1]+'/' #the input should be the target name just like all of the other times this has been implemented in this repository
#path to the model files
teff= "013500" #hard coded for test target
#logg= 8.0 #model values to pull (should be changed to be read from somewhere in the future)
logg = "800" #hard coded for test target(it was 7.8 in the paper, but that is technically the same as this gridpoint to two sig figs)
#min_wave= 1142
#max_wave= 1422
min_wave= 1180
max_wave= 1800
#lyman_alpha = [1214,1217] #probably should put this into the config file as the masking for the spectroscopic fitting (but really this goes back to the need for individual text files by target for the various maskings that are required)
#oxygen= [1300, 1308]
#nitrogen=[1199,1201]
lyman_alpha = config.lyman_mask
oxygen= config.oxygen_mask
nitrogen = config.nitrogen_mask
#seg_gap = [1268,1298] #segment gap mask used for photometry for the G130M at 1291 Angs
#scale_wave_range= [1400, 1420] 
scale_wave_range= [1450,1550] 
#scale factor to be just typed in for the moment but should be replaced by a function at some point that calculates it for you
target_x1ds= glob(target_path+'*x1dsum.fits')
target_file= target_x1ds[0] #I only want the first one for now; when I eventually make this include all the visits for a given target I'll introduce some sort of function for combining them or whatever
model_path= '/Users/BenKaiser/Desktop/DK_2010/'
model_extension = '.dat'
#scaling_coefficient= 5.2e21
#wavelength range to consider for the scaling

target_hdu= fits.open(target_file)
cenwave=target_hdu[0].header['CENWAVE']
grating = target_hdu[0].header['OPT_ELEM']
seg_gap = config.seg_gap_dict[str(grating)][str(cenwave)]
print "grating (OPT_ELEM):", grating, "Centwave: ", target_hdu[0].header['CENWAVE']
target_waves = np.copy(target_hdu[1].data['wavelength'].ravel())
target_flux= np.copy(target_hdu[1].data['flux'].ravel())
target_spec= np.vstack([target_waves, target_flux])
#print target_spec.shape
#print target_spec[0]


#function to trim the model plot to limit it to the range of the target's spectrum
#At some point I'll have to implement masking


def get_model_spec(fitted_teff, fitted_logg):
    """
    inputs: fitted_teff, fitted_logg
    
    output: 2-d numpy array containing the wavelengths and spectral values from the model grid for the input teff and logg
    """
    padded_teff = fitted_teff #needs to be changed to pad/remove decimals (and convert to string)
    padded_logg= fitted_logg #needs to be changed to pad/remove decimals (and convert to string)
    model_file= 'da' +padded_teff + '_' + padded_logg + model_extension
    model_file= model_path + model_file
    model_file= glob(model_file)[0]
    model_spec= np.genfromtxt(model_file)
    model_spec= model_spec.T #transpose because numpy is weird on the axes it reads in.
    #print "model_spec.shape:",  model_spec.shape
    return model_spec

def get_model_fromfile(model_file):
    """
    inputs: model file name
    
    output: 2-d numpy array containing the wavelengths and spectral values from the model grid for the input teff and logg
    """
    model_spec= np.genfromtxt(model_file)
    model_spec= model_spec.T #transpose because numpy is weird on the axes it reads in.
    return model_spec

#########

def get_scale_factor(target_spec, model_spec, wave_range):
    allowed_target = np.where((target_spec[0] > wave_range[0])& (target_spec[0] < wave_range[1]))
    allowed_model = np.where((model_spec[0] > wave_range[0]) & (model_spec[0] < wave_range[1]))
    target_mean= np.mean(target_spec[1][allowed_target])
    model_mean = np.mean(model_spec[1][allowed_model])
    scaling_factor =target_mean/model_mean
    return scaling_factor
####

def calc_sq_dist(target_spec, model_spec):
    interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    #print "interp_model.shape", interp_model.shape
    norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(interp_model[1])
    nan_remove = np.isinf(norm_difs)
    norm_difs= norm_difs[~nan_remove]
    dif = np.sum(norm_difs)/norm_difs.shape[0]
    
    return dif
    
###
def plot_overlays(spec1, spec2):
    plt.plot(spec1[0], spec1[1], label = 'observed')
    plt.plot(spec2[0], spec2[1], label= 'model', color = 'r')
    plt.legend(numpoints=1, fontsize=14, loc='best' )
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (cgs units)')
    plt.title(target_path[:-1] )
    plt.show()
    return ''


        


##########
#model_spec= get_model_spec(teff, logg)
#target_spec= remove_range(target_spec[0], target_spec[1], lyman_alpha)
#target_spec= remove_range(target_spec[0], target_spec[1], oxygen)
#target_spec= remove_range(target_spec[0], target_spec[1], seg_gap)
#target_spec = sort_spectrum(target_spec)
print os.getcwd()
def run_model_grid(target_spec):
    mask_list= [lyman_alpha]+[oxygen]+[seg_gap]+[nitrogen]
    target_spec = spt.clean_spectrum(target_spec, min_wave, max_wave, mask_list)
    model_file_list = glob(model_path+'*'+model_extension)
    dist_list = []
    for model_file in model_file_list:
        model_spec = get_model_fromfile(model_file)
        model_spec= spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
        scaling_coefficient= get_scale_factor(target_spec, model_spec, scale_wave_range)
        model_spec[1]=model_spec[1]*scaling_coefficient
        new_dist = calc_sq_dist(target_spec, model_spec)
        dist_list.append(new_dist)
       
    dist_array = np.array(dist_list)
    for mod_file, dist_mod in zip(model_file_list, dist_list):
        print "model_file:", mod_file, "difference:", dist_mod
    min_index = np.argmin(dist_list)
    min_model = model_file_list[min_index]
    min_dist = dist_array[min_index]
    print "best fit model:", min_model
    model_spec= get_model_fromfile(min_model)
    model_spec = spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
    scaling_coefficient= get_scale_factor(target_spec, model_spec, scale_wave_range)
    model_spec[1]= model_spec[1]*scaling_coefficient
    
    plot_overlays(target_spec, model_spec)
    
run_model_grid(target_spec)
