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
#scale_wave_range= [1400, 1420] 
scale_wave_range= [1700,1750] 
lyman_alpha = config.lyman_mask
oxygen= config.oxygen_mask
nitrogen = config.nitrogen_mask
target_x1ds= glob(target_path+'*x1dsum.fits')
target_file= target_x1ds[-1] #I only want the first one for now; when I eventually make this include all the visits for a given target I'll introduce some sort of function for combining them or whatever
model_path= '/Users/BenKaiser/Desktop/DK_2010/'
model_extension = '.dat'
target_hdu= fits.open(target_file)
cenwave=target_hdu[0].header['CENWAVE']
grating = target_hdu[0].header['OPT_ELEM']
seg_gap = config.seg_gap_dict[str(grating)][str(cenwave)]
print "grating (OPT_ELEM):", grating, "Centwave: ", target_hdu[0].header['CENWAVE']
#exposure_time = target_hdu[0].header['EXPTIME']
target_waves = np.copy(target_hdu[1].data['wavelength'].ravel())
target_flux= np.copy(target_hdu[1].data['flux'].ravel())
target_error = np.copy(target_hdu[1].data['ERROR'].ravel())
target_spec= np.vstack([target_waves, target_flux])
target_err = np.vstack([target_waves, target_error])
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
    target_median= np.median(target_spec[1][allowed_target])
    model_median = np.median(model_spec[1][allowed_model])
    scaling_factor =target_median/model_median
    return scaling_factor
####

def calc_sq_dist(target_spec, model_spec, error_spec = np.array([])):
    interp_model_flux = np.interp(target_spec[0], model_spec[0], model_spec[1])
    interp_model= np.vstack([np.copy(target_spec[0]),interp_model_flux])
    #print "interp_model.shape", interp_model.shape
    if error_spec.shape[0] != 0:
        #norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(error_spec[1])
        norm_difs = (interp_model[1]-target_spec[1])**2/np.float_(error_spec[1])**2
        #norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(interp_model[1])
    else:
        print "no uncertainties provided"
        norm_difs = np.abs(interp_model[1]-target_spec[1])/np.float_(interp_model[1])
    #norm_difs = np.abs(interp_model[1]-target_spec[1])

    nan_remove = np.isinf(norm_difs)
    norm_difs= norm_difs[~nan_remove]
    dif = np.sum(norm_difs)/norm_difs.shape[0]
    
    return dif
    
###
def plot_overlays(spec1, spec2, errors):
    plt.plot(spec1[0], spec1[1], label = 'observed')
    #plt.errorbar(spec1[0],spec1[1], yerr = errors[1], label='observed')
    plt.plot(spec2[0], spec2[1], label= 'model', color = 'r')
    plt.legend(numpoints=1, fontsize=14, loc='best' )
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux (cgs units)')
    plt.title(target_path[:-1] )
    plt.show()
    return ''

def make_model_params_arrays(list_model_files):
    teff_list = []
    logg_list = []
    for file_path in list_model_files:
        file_name = file_path.split('/')[-1] #takes the characters after the last slash
        sub_string = file_name.split('da')[1] #the not-empty portion
        divided_string = sub_string.split('_')
        teff = int(divided_string[0] ) #the effective temperature for the model
        logg = int(divided_string[1].split('.')[0])
        teff_list.append(teff)
        logg_list.append(logg)
    teff_array = np.array(teff_list)
    logg_array = np.array(logg_list)
    return teff_array, logg_array

def chi_square_countours(teff_array, logg_array, dist_array):
    min_index = np.argmin(dist_array)
    print "Teff and logg min chi-squared values: ", teff_array[min_index],logg_array[min_index], "|chi-sq:", dist_array[min_index]
    contour_array = np.vstack([teff_array,logg_array, dist_array])
    #plt.imshow(contour_array, aspect= 100)
    #plt.contour(teff_array, logg_array, dist_array)
    marker_scale = 1/dist_array* dist_array.min() *40.
    #plt.scatter(teff_array, logg_array, s= 1./dist_array*30, c = 1./dist_array*20)
    plt.scatter(teff_array, logg_array, s=marker_scale, c = marker_scale)
    plt.plot(teff_array[min_index],logg_array[min_index], marker = '*', markersize = 14)
    plt.xlabel('T_eff')
    plt.ylabel('logg')
    plt.show()

def calc_rdist(scale_factor):
    dobs_over_dmod = 1./np.sqrt(scale_factor)
    assumed_model_dist = 1e-2 #kpc
    print "Target is ", dobs_over_dmod, "times farther away from us than the model"
    print "d_target ~", dobs_over_dmod*assumed_model_dist, "kpc (assuming d_obs = "+str(assumed_model_dist *1000)+")"

##########
#model_spec= get_model_spec(teff, logg)
#target_spec= remove_range(target_spec[0], target_spec[1], lyman_alpha)
#target_spec= remove_range(target_spec[0], target_spec[1], oxygen)
#target_spec= remove_range(target_spec[0], target_spec[1], seg_gap)
#target_spec = sort_spectrum(target_spec)
print os.getcwd()
def run_model_grid(target_spec,target_err=None):
    mask_list= [lyman_alpha]+[oxygen]+[seg_gap]+[nitrogen]
    target_spec = spt.clean_spectrum(target_spec, min_wave, max_wave, mask_list)
    target_err = spt.clean_spectrum(target_err, min_wave, max_wave, mask_list)
    model_file_list = glob(model_path+'*'+model_extension)
    dist_list = []
    for model_file in model_file_list:
        model_spec = get_model_fromfile(model_file)
        model_spec= spt.trim_spec(model_spec, np.min(target_spec[0]), np.max(target_spec[0]))
        scaling_coefficient= get_scale_factor(target_spec, model_spec, scale_wave_range)
        model_spec[1]=model_spec[1]*scaling_coefficient
        new_dist = calc_sq_dist(target_spec, model_spec, error_spec = target_err)
        dist_list.append(new_dist)
    teff_array, logg_array = make_model_params_arrays(model_file_list)
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
    calc_rdist(scaling_coefficient)
    plot_overlays(target_spec, model_spec, target_err)
    chi_square_countours(teff_array,logg_array, dist_array)
    
    
run_model_grid(target_spec, target_err)
