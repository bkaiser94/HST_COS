import numpy as np
import os
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
import sys
import config

#from calcos import calcos
from costools import timefilter
import local_lightcurve as lc


def plot_all_x1d(target_dir, low_lim, high_lim):
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
        print dataset, hdu[0].header['OPT_ELEM'], hdu[0].header['CENWAVE'], hdu[1].header['EXPTIME']    

    fig = plt.figure(figsize=(200,9))
    #silicon_lines= [1190, 1193, 1195, 1197, 1207, 1260, 1265, 1304, 1309]
    ax1 = fig.add_subplot(1,1,1)

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
        #fluxes= np.copy(hdu[1].data['Gcounts'].ravel())

        #print "fluxes",  fluxes
        arrshape= wavelengths.shape
        #unmasked= np.where((wavelengths < lyman[0]*np.ones(arrshape)) or ((wavelengths > lyman[1]*np.ones(arrshape)) and (wavelengths < oxygen[0])*np.ones(arrshape)) or (wavelengths > oxygen[1]*np.ones(arrshape)))
        lower_mask = np.where(wavelengths > wave_min)
        wavelengths= np.copy(wavelengths[lower_mask])
        fluxes = np.copy(fluxes[lower_mask])
        upper_mask= np.where(wavelengths < wave_max)
        wavelengths= np.copy(wavelengths[upper_mask])
        fluxes= np.copy(fluxes[upper_mask])
        wavelengths, fluxes = remove_range(wavelengths, fluxes, config.lyman_mask)
        wavelengths, fluxes = remove_range(wavelengths, fluxes, config.oxygen_mask)
        wavelengths, fluxes = remove_range(wavelengths, fluxes, config.nitrogen_mask)
        
        print "fluxes.shape: ", fluxes.shape
        print "wavelengths.shape" , wavelengths.shape
        
        if counter == 0:
            flux_all= np.array([np.zeros(fluxes.shape)])
            print "flux_all.shape: ", flux_all.shape
            plot_waves = wavelengths
        try:
            flux_all = np.append(flux_all, [fluxes], axis=0)
            print flux_all.shape
        except ValueError as error:
            print error
            shape_dif = flux_all.shape[1] - fluxes.shape[0]
            if shape_dif > 0 :
                padarray= np.zeros(shape_dif)
                fluxes_pad= np.append(fluxes,padarray)
                flux_all = np.append(flux_all, [fluxes_pad], axis=0)
            elif shape_dif < 0 :
                flux_all= np.append(flux_all,[fluxes[:shape_dif]], axis=0)
                print "truncating new spectrum to mach previous dimensions"
            
        #try:
            ##flux_all = np.copy(fluxes+flux_all)
        #except ValueError as error:
            #print error
            #shape_dif = flux_all.shape[0] - fluxes.shape[0]
            #if shape_dif > 0 :
                #padarray= np.zeros(shape_dif)
                #flux_all= np.copy(flux_all + np.append(fluxes, padarray))
            #elif shape_dif < 0 :
                #flux_all= np.copy(flux_all + fluxes[:shape_dif])
        counter += 1
        #print "sum fluxes: ", np.sum(fluxes)
        ax1.plot(wavelengths, fluxes, label=hdu[1].header['EXPSTART'])
    for this_line in config.silicon_lines:
        ax1.axvline(x= this_line ,linestyle = '-', color = 'g' , ymin = 0, ymax = 100000, linewidth = 1, alpha = 0.2)
    print flux_all.shape
    flux_all= flux_all[1:, :] #remove the first row of zeros
    flux_med= np.nanmedian(flux_all, axis = 0)
    #flux_med= np.sum(flux_all, axis=0)
    ax1.plot(plot_waves, flux_med, label= 'median combined spectra')
    ax1.legend(numpoints=1, fontsize=14, loc='best' )
    ax1.set_ylabel('Flux (normed)')
    ax1.set_xlabel('Wavelength $(\AA)$')
    ax1.set_title(target_dir)
    #ax1.set_yscale('log')
    fig.tight_layout()
    fig.savefig(target_dir+'x1dtruesum.pdf')
    plt.show()
    
    
if __name__ == '__main__':
    target_dir= sys.argv[1]
    #lcfile= sys.argv[2]
    low_lim = int(sys.argv[2])
    high_lim = int(sys.argv[3])
    #period= float(sys.argv[3])
    #unit_arg = sys.argv[4]
    #nullstring = make_dual_plots(target_dir, stepsize, period, unit_arg)
    nullstring = plot_all_x1d(target_dir, low_lim, high_lim)
