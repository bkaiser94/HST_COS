import numpy as np
import os
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import config

#from calcos import calcos
from costools import timefilter
import local_lightcurve as lc
import spec_plot_tools as spt

target_dir= sys.argv[1]

#wave_limits= [1130,1800]
def run_spectral_exam(wave_limits):
    
    wave_min= wave_limits[0]
    #lyman = [1208, 1225]
    #oxygen= [1295, 1312] #airglow wavelengths to be filtered according to lightcurve
    wave_max= wave_limits[1]
    #lcbase= target_dir+ '_grid_lightcurves/' + '*step' + str(stepsize) + '_*'+str(wave_min)+',' + str(wave_max)+'*'
    #lcfile= glob(lcbase)[0]
    dest_dir = 'segmented_spectra/'+target_dir +'/'
    if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)


    for dataset in glob(target_dir+ '/*x1dsum.fits'):
        hdu = fits.open(dataset)
        print dataset, hdu[0].header['OPT_ELEM'], hdu[0].header['CENWAVE'], hdu[1].header['EXPTIME']    

    def plot_element_lines(wavelength_array, axis_variable):
        for this_line in config.silicon_lines:
            if ((this_line > np.nanmin(wavelength_array)) & (this_line< np.nanmax(wavelength_array))):
                axis_variable.axvline(x= this_line ,linestyle = '-', color = 'g' , ymin = 0, ymax = 100000, linewidth = 1, alpha = 0.9)


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

    def zero_mask(wave_array, other_array, bound_list):
        lower_bound = bound_list[0]
        upper_bound = bound_list[1]
        masking = np.where((wave_array > lower_bound) & (wave_array < upper_bound))
        other_array[masking]= 0.
        return other_array

    counter= 0
    vertical_dimension = 6
    max_points_per_window= 501

    for dataset in glob(target_dir+ '/*x1dsum.fits'):
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
        target_spec= np.vstack([wavelengths, fluxes])
        mask_list= [config.lyman_mask]+[config.oxygen_mask]+[config.nitrogen_mask] #need to add the segment gap in the future
        cleaned_spec = spt.clean_spectrum(target_spec, wave_min, wave_max, mask_list)
        wavelengths= cleaned_spec[0]
        fluxes= cleaned_spec[1]
        #lower_mask = np.where(wavelengths > wave_min)
        #wavelengths= np.copy(wavelengths[lower_mask])
        #fluxes = np.copy(fluxes[lower_mask])
        #upper_mask= np.where(wavelengths < wave_max)
        #wavelengths= np.copy(wavelengths[upper_mask])
        #fluxes= np.copy(fluxes[upper_mask])
        #print "fluxes.shape: ", fluxes.shape
        #wavelengths, fluxes = remove_range(wavelengths, fluxes, config.lyman_mask)
        #wavelengths, fluxes = remove_range(wavelengths, fluxes, config.oxygen_mask)
        #wavelengths, fluxes = remove_range(wavelengths, fluxes, config.nitrogen_mask)
        #fluxes= zero_mask(wavelengths, fluxes, config.lyman_mask)
        #fluxes= zero_mask(wavelengths, fluxes, config.oxygen_mask)
        
        print "fluxes.shape: ", fluxes.shape
        print "wavelengths.shape" , wavelengths.shape
        
        if counter == 0:
            flux_all= np.array([np.zeros(fluxes.shape)])
            print "flux_all.shape: ", flux_all.shape
            plot_waves = wavelengths
        print "wavelengths differences", plot_waves[:5] - wavelengths[:5]
        print plot_waves[:5]
        print wavelengths[:5]
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
                print "truncating new spectrum to match previous dimensions"
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
        #ax1.plot(wavelengths, fluxes/np.nanmean(fluxes), label=hdu[0].header['rootname'])
    num_segs= plot_waves.shape[0] / max_points_per_window +1

    #silicon_lines= [1190, 1193, 1195, 1197, 1207, 1260, 1265, 1304, 1309]
    flux_all= flux_all[1:, :] #remove the first row of zeros
    flux_med= np.nanmedian(flux_all, axis = 0)
    #flux_normed = flux_med/np.nanmean(flux_med) #changed this so it's actually not normalized
    flux_normed= flux_med
    needed_zeros= max_points_per_window- plot_waves.shape[0] % max_points_per_window
    flux_normed = np.append(flux_normed, np.zeros(needed_zeros))
    dlambda= plot_waves[-1]-plot_waves[-2]    #find the wavelength bin size
    plus_lambda = dlambda* np.indices([needed_zeros,])[0] #find vaues that need to be added to get out to the appropriate range
    plus_lambda = plus_lambda + plot_waves[-1]
    print "plus_lambda.shape:", plus_lambda.shape
    plot_waves= np.append(plot_waves, plus_lambda) #tack on the necessary wavelength values on the far end of the spectrum to get it farther without having to deal with all of the associated noise.
    #now reshape both plot_waves and flux_normed to be multidimensional arrays that we can iterate through.
    #output a long version of the spectrum as well
    longfig= plt.figure(figsize= (200, 12))
    longax= longfig.add_subplot(1,1,1)
    longax.plot(plot_waves, flux_normed, color= 'k')
    plot_element_lines(plot_waves, longax)
    plt.grid()
    #longax.set_ylabel('Flux (normed)')
    longax.set_ylabel('Flux (cgs units)')
    longax.set_xlabel('Wavelength $(\AA)$')
    longfig.savefig(dest_dir+target_dir+'long_spectrum_'+ str(wave_limits[0])+','+str(wave_limits[1]) + '.pdf')

    plot_waves= plot_waves.reshape(num_segs, max_points_per_window)
    flux_normed= flux_normed.reshape(num_segs, max_points_per_window)

    num_mults_42= num_segs/42 #42 is the maximum number of subplots allowed

    figsaves= 0
    if num_segs > 42:
        for gorounds in range(0,num_mults_42):
            offset = 42* gorounds
            new_max= 42
            fig = plt.figure(figsize=(20,new_max* vertical_dimension))
            for segment in range(0,new_max):
                ax = fig.add_subplot(new_max,1,segment+1)
                waves_for_plot= plot_waves[segment+offset]
                flux_for_plot = flux_normed[segment+offset]
                plot_element_lines(waves_for_plot, ax)
                #for this_line in config.silicon_lines:
                    #if ((this_line > np.nanmin(waves_for_plot)) & (this_line< np.nanmax(waves_for_plot))):
                        #ax.axvline(x= this_line ,linestyle = '-', color = 'g' , ymin = 0, ymax = 100000, linewidth = 1, alpha = 0.2)
                #print flux_all.shape
                ax.plot(waves_for_plot, flux_for_plot, color = 'k')
                #ax.set_ylabel('Flux (normed)')
                ax.set_ylabel('Flux (cgs units)')
                ax.set_xlabel('Wavelength $(\AA)$')
                ax.set_xlim([np.nanmin(waves_for_plot), np.nanmax(waves_for_plot)])


            fig.tight_layout()
            plt.title(target_dir)
            fig.savefig(dest_dir+target_dir+ 'segmented_spectra_wlim'+ str(wave_limits[0])+','+str(wave_limits[1]) +'_fig'+str(figsaves)+ '.pdf')
            figsaves +=1
            #plt.show()
        #and now the remaining segments:
        offset= 42*num_mults_42
        fig = plt.figure(figsize=(20,num_segs* vertical_dimension))
        for segment in range(0,num_segs%42):
            ax = fig.add_subplot(num_segs,1,segment+1)
            waves_for_plot= plot_waves[segment+offset]
            flux_for_plot = flux_normed[segment+offset]
            for this_line in config.silicon_lines:
                if ((this_line > np.nanmin(waves_for_plot)) & (this_line< np.nanmax(waves_for_plot))):
                    ax.axvline(x= this_line ,linestyle = '-', color = 'g' , ymin = 0, ymax = 100000, linewidth = 1, alpha = 0.2)
            #print flux_all.shape
            ax.plot(waves_for_plot, flux_for_plot, color ='k')
            ax.set_ylabel('Flux (normed)')
            ax.set_xlabel('Wavelength $(\AA)$')
            ax.set_xlim([np.nanmin(waves_for_plot), np.nanmax(waves_for_plot)])
            

        fig.tight_layout()
        plt.title(target_dir)
        fig.savefig(dest_dir+target_dir+ 'segmented_spectra_wlim'+ str(wave_limits[0])+','+str(wave_limits[1]) +'_fig'+str(figsaves)+ '.pdf')
        figsaves +=1
    else:
        print "Currently not well-configured to handle fewer than 42 segments"
        print "That corresponds to " + str(max_points_per_window * 42) + ' data points in the set.'
        fig = plt.figure(figsize=(20,num_segs* vertical_dimension))
        for segment in range(0,num_segs):
            ax = fig.add_subplot(num_segs,1,segment+1)
            waves_for_plot= plot_waves[segment]
            flux_for_plot = flux_normed[segment]
            for this_line in config.silicon_lines:
                if ((this_line > np.nanmin(waves_for_plot)) & (this_line< np.nanmax(waves_for_plot))):
                    ax.axvline(x= this_line ,linestyle = '-', color = 'g' , ymin = 0, ymax = 100000, linewidth = 1, alpha = 0.2)
        #print flux_all.shape
            ax.plot(waves_for_plot, flux_for_plot, color ='k')
            ax.set_ylabel('Flux (normed)')
            ax.set_xlabel('Wavelength $(\AA)$')
            ax.set_xlim([np.nanmin(waves_for_plot), np.nanmax(waves_for_plot)])


        fig.tight_layout()
        plt.title(target_dir)
        fig.savefig(dest_dir+target_dir+ 'segmented_spectra_wlim'+ str(wave_limits[0])+','+str(wave_limits[1]) +'_fig'+str(figsaves)+ '.pdf')
        figsaves +=1
        
for wave_limits in config.wave_limit_list:
    try:
        run_spectral_exam(wave_limits)
    except IndexError as error:
        print "Issue producing a spectrum for " + str(wave_limits)
        print error




