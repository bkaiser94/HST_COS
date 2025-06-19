from __future__ import print_function


import numpy as np
import os
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from astropy import convolution as conv
import sys
import config

#from calcos import calcos
#from costools import timefilter
import local_lightcurve as lc
#import spec_plot_tools as spt

sys.path.append('/Users/BenKaiser/Desktop/radial_velocity_calculations/')
import spec_plot_tools as spt

plt.rc('lines',linewidth=0.5)

fp_pos_colors = ['magenta', 'r', 'g', 'b', 'k']
def plot_all_x1d(target_dir, low_lim, high_lim, log_scale):
    wave_min= low_lim
    #lyman = [1208, 1225]
    #oxygen= [1295, 1312] #airglow wavelengths to be filtered according to lightcurve
    wave_max= high_lim
    #lcbase= target_dir+ '_grid_lightcurves/' + '*step' + str(stepsize) + '_*'+str(wave_min)+',' + str(wave_max)+'*'
    #lcfile= glob(lcbase)[0]
    #dest_dir = 'dual_plots/'
    #if not os.path.exists(dest_dir):
                #os.makedirs(dest_dir)
    

    for dataset in glob(target_dir+ '/*x1d.fits'):
        hdu = fits.open(dataset)
        #try:
            #print(dataset, hdu[0].header['OPT_ELEM'], hdu[0].header['CENWAVE'], hdu[1].header['EXPTIME'], "LP-POS: ", hdu[0].header['LIFE_ADJ'], "FP-POS: ", hdu[0].header['FPPOS'])
        #except KeyError as error:
            #print('KeyError:', error)
            #print('filename that produced error:', dataset)
            #print(dataset, hdu[0].header['OPT_ELEM'], hdu[0].header['CENWAVE'], hdu[1].header['EXPTIME'], "LP-POS: ", "FP-POS: ", hdu[0].header['FPPOS'])

    #fig = plt.figure(figsize=(200,9))
    fig = plt.figure()
    #silicon_lines= [1190, 1193, 1195, 1197, 1207, 1260, 1265, 1304, 1309]
    ax1 = fig.add_subplot(1,1,1)

    counter= 0
    mask_list= [config.lyman_mask]+[config.oxygen_mask]+[config.nitrogen_mask] #need to add the segment gap (not implemented directly in spec_plot_tools because it will eventually be called from a dict based on grating and central wavelength
    for mask in mask_list:
        plt.axvspan(mask[0],mask[1],alpha=0.1, color='k')
    for dataset in glob(target_dir+ '/*x1d.fits')[:6]:
        #fig = plt.figure()
        #ax1 = fig.add_subplot(1,1,1)
        hdu = fits.open(dataset)
        #print(hdu[1])
        #print(hdu[1]['FUVA'])
        #print(hdu[1].data[0])
        #print(hdu[1].data.shape)
        #for thing in hdu[1].data:
            #print("thing: ", thing)
            #print(thing['wavelength'])
        #fppos = hdu[0].header['FPPOS']
        #fppos_color = fp_pos_colors[fppos]
        wavelengths= np.copy(hdu[1].data['wavelength'].ravel())
        print(wavelengths.shape)
        fluxes= np.copy(hdu[1].data['flux'].ravel())
        #fluxes= np.copy(hdu[1].data['Gcounts'].ravel())

        #print("fluxes",  fluxes)
        arrshape= wavelengths.shape
 
        target_spec= np.vstack([wavelengths, fluxes])
        #mask_list= [config.lyman_mask]+[config.oxygen_mask]+[config.nitrogen_mask] #need to add the segment gap (not implemented directly in spec_plot_tools because it will eventually be called from a dict based on grating and central wavelength
        #cleaned_spec = spt.clean_spectrum(target_spec, wave_min, wave_max, mask_list)
        pix_kernel = conv.Box1DKernel(width = int(5), mode = 'oversample')
        pix_kernel.normalize()
        target_spec[1] = conv.convolve(target_spec[1], pix_kernel)
        if np.nanmin(target_spec[0])>1500:
            cleaned_spec = spt.clean_spectrum(target_spec, 1700, 3141, [])
        else:
            cleaned_spec = spt.clean_spectrum(target_spec, 1145, 3141, [])  
            #cleaned_spec=target_spec
        fluxes= cleaned_spec[1]
        wavelengths= cleaned_spec[0]
        #print("fluxes.shape: ", fluxes.shape)
        #print("wavelengths.shape" , wavelengths.shape)
        
        if counter == 0:
            flux_all= np.array([np.zeros(fluxes.shape)])
            print("flux_all.shape: ", flux_all.shape)
            plot_waves = wavelengths
        try:
            flux_all = np.append(flux_all, [fluxes], axis=0)
            print(flux_all.shape)
        except ValueError as error:
            print(error)
            shape_dif = flux_all.shape[1] - fluxes.shape[0]
            if shape_dif > 0 :
                padarray= np.zeros(shape_dif)
                fluxes_pad= np.append(fluxes,padarray)
                flux_all = np.append(flux_all, [fluxes_pad], axis=0)
            elif shape_dif < 0 :
                flux_all= np.append(flux_all,[fluxes[:shape_dif]], axis=0)
                print("truncating new spectrum to mach previous dimensions")
            
        
        counter += 1
        #print("sum fluxes: ", np.sum(fluxes))
        #ax1.plot(wavelengths, fluxes, label=hdu[1].header['EXPSTART'], color = fppos_color )
        #ax1.plot(wavelengths, fluxes, label=hdu[1].header['EXPSTART'],marker='.', linestyle='None',markersize=2)
        ax1.plot(wavelengths, fluxes, label=hdu[1].header['EXPSTART'])

        #ax1.set_ylabel('Flux (cgs)')
        #ax1.set_xlabel('Wavelength $(\AA)$')
        #ax1.set_title(target_dir)
        #plt.ylim(0,)
        #plt.show()
    #for this_line in config.silicon_lines:
        #if ((this_line >= wave_min) & (this_line <= wave_max)):
            #ax1.axvline(x= this_line ,linestyle = '-', color = 'g' , ymin = 0, ymax = 100000, linewidth = 1, alpha = 0.2)
    print(flux_all.shape)
    flux_all= flux_all[1:, :] #remove the first row of zeros
    flux_med= np.nanmedian(flux_all, axis = 0)
    #flux_med= np.sum(flux_all, axis=0)
    #ax1.plot(plot_waves, flux_med, label= 'median combined spectra', color =fp_pos_colors[0])
    #ax1.legend(numpoints=1, fontsize=14, loc='best' )
    ax1.set_ylabel('Flux (cgs)')
    ax1.set_xlabel('Wavelength $(\AA)$')
    ax1.set_title('GD 356 STIS data')
    #mwdd_filename='GD356_DAHe_IUE_fnu_MWDD_spec.csv'
    #mwdd_spec, mwdd_err_spec=spt.retrieve_mwdd_spec(mwdd_filename,convert_to_flambda=True,iue_spec=True)

    #pix_kernel = conv.Box1DKernel(width = int(5), mode = 'oversample')
    #pix_kernel.normalize()
    #mwdd_spec[1] = conv.convolve(mwdd_spec[1], pix_kernel)
    #plt.plot(mwdd_spec[0],mwdd_spec[1]*1e-16, label='IUE',linestyle=':')
    plt.legend()
    if log_scale:
        ax1.set_yscale('log')
        #bad_inds = np.where(
    #fig.tight_layout()
    #fig.savefig(target_dir+'x1dtruesum.pdf')
    plt.show()
    
    
if __name__ == '__main__':
    target_dir= sys.argv[1]
    #lcfile= sys.argv[2]
    #low_lim = int(sys.argv[2])
    #high_lim = int(sys.argv[3])
    log_scale = False
    try:
        scale_log = sys.argv[4]
        if scale_log.lower() == 'true':
            log_scale = True
    except IndexError:
        pass
    #period= float(sys.argv[3])
    #unit_arg = sys.argv[4]
    #nullstring = make_dual_plots(target_dir, stepsize, period, unit_arg)
    low_lim=0
    high_lim=5000
    nullstring = plot_all_x1d(target_dir, low_lim, high_lim, log_scale)
