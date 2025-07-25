"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-06-19

This is going to make lightcurves from STIS data. I duplicated try_lightcurve.py, which I made 8 years ago to run on COS data, so it is a bit of a mess. I'm planning to do some cleaning in here to handle the STIS data, but I didn't want to ruin the working COS light curve generation functions.

Seemingly has to run in Python2 because lightcurve needs calcos.costools, but I guess newer versions of calcos don't have that module, so it doesn't work.
"""


from __future__ import print_function


import numpy as np
import os
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.time import Time
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Table, Column
from calcos import calcos
from costools import timefilter
#import lightcurve as lc
import local_lightcurve as lc #I have created a second directory within here that I could edit the 

import config
#target_dir = sys.argv[1] + '/'
#stepsize = int(sys.argv[2]) 
second_per_mjd= 1./1.15741e-5     #SECOND_PER_MJD value from lightcurve.cos.extract, but I couldn't import it for whatever reason... so I just copied and pasted
start_dir = os.getcwd()
print("start_dir: ", start_dir)
#def make_lightcurve(target_dir, stepsize, wlim, plotall=True, detector='FUV'):
def make_lightcurve(target_dir, stepsize, wlim, plotall=True, detector='FUV'):
    """
    """
    print("")
    print("===============")
    print("stepsize: ", stepsize)
    print("plotall=", plotall)
    print("target_dir: ", target_dir)
    #wlim= [1130, 1900] #3200 is the default max from what I could gather as is 915. The geocoronal emission lines are automatically excluded from the lightcurve
    if detector=='FUV':
        band= '_a' #should be '_a', '_b', '' for the combined band images, or  '*' to get all bands available (the '*' option actually causes cos.py to double-count bands)
    elif detector=='NUV':
        band=''
    else:
        print("Apparently there's not a detector header...or it didn't make it into make_lightcurve()...")

    mjd_array = np.array([])

    gross_array = np.array([])
    flux_table = np.array([])
    bkg_array=np.array([])
    count_array=np.array([])
    error_array=np.array([])
    ferror_array=np.array([])
    negative_flux_files= []
    negative_means= []
    if plotall:
        fig = plt.figure(figsize= (20,9))
        ax = fig.add_subplot(1,1,1)


    print(os.getcwd())
    assembled_files= glob(target_dir + '*corrtag' + band + '.fits')
    header=fits.getheader(assembled_files[0])
    for item in assembled_files:
        astrotable, astrometa = lc.cos.extract(item, step=stepsize, wlim= wlim)
        #print('astrometa: ', astrometa)
        
        mjd_array= np.append(mjd_array, astrotable['mjd'])
        gross_array= np.append(gross_array, astrotable['gross'])
        flux_table= np.append(flux_table, astrotable["flux"]) #for complete dataset mean normalization
        bkg_array=np.append(bkg_array, astrotable['background'])
        count_array=np.append(count_array, astrotable['gross']-astrotable['background'])
        error_array=np.append(error_array, np.sqrt(astrotable['gross']+astrotable['background']))
        ferror_array= np.append(ferror_array, np.sqrt(astrotable['gross'] + astrotable['background']) / (astrotable['gross'] - astrotable['background']) * astrotable['flux'])
        #flux_table= np.append(flux_table, (astrotable["flux"]/np.nanmedian(astrotable['flux'])-1)) #for individual fits file median normalization
        #flux_table= np.append(flux_table, (astrotable["flux"]/np.nanmean(astrotable['flux'])-1)) #for individual fits file mean normalization
        print("first mjd: " , astrotable['mjd'][0])
        #print("astrotable['gross']: " , astrotable['gross'])
        #print("astrotable['background']: ", astrotable['background'])
        #count_diff= astrotable['gross']-astrotable['background']
        #print("difference: ", count_diff)
        #print("normed counts: ", np.float_(count_diff)/np.nanmean(count_diff))
        print("file: ", item)
        print("np.nanmean(flux):", np.nanmean(astrotable['flux']))
        print('-------------')
        if np.nanmean(astrotable['flux'] < 0.):
            negative_flux_files.append(item)
            negative_means.append(np.nanmean(astrotable['flux']))
        
    print('===========')
    print("Negative flux files and means detected: ")
    for filename, nanmean in zip(negative_flux_files, negative_means):
        print("File: ", filename, "| nanmean: ", nanmean)
    print('===========')


    pre_flux = np.copy(flux_table)
    ferror_array= np.copy(ferror_array/flux_table) #normalizing the error to the flux values to whatever units flux happens to be in at a given time
    norm_ferror_array= np.copy(ferror_array*(flux_table/np.nanmean(flux_table))) #making the flux error array be in scaled units to the normalized flux that's about to be calculated
    #it had to be done without subtracting 1 because then you'd have negative values for the errors inexplicably
    #this is also why it couldn't be left in arbitrary units because it would later be multiplied by the negative flux values
    flux_array = np.copy( flux_table/np.nanmean(flux_table)-1.) #normalize flux array around zero #This was the mean all norming method
    #flux_array = np.copy(flux_table) #this should be uncommented for the individual fits normalization
    #flux_array= np.copy(flux_table/np.nanmedian(flux_table)-1.) #New version as of 2017-11-30
    #print("np.nanmean(flux_array)" , np.nanmean(flux_array))
    #print("np.nanmean(flux_table: " ,np.nanmean(flux_table))
    #print("np.sum(flux_array-pre_flux): ", np.sum(flux_array- pre_flux))
    #print("np.sum(flux_array + pre_flux): ", np.sum(flux_array + pre_flux))
    #print("")
    #print("np.nanmax(flux_array): ", np.nanmax(flux_array))
    #print("mjd value for that: ", mjd_array[np.argmax(flux_array)])
    #print("top 5 flux values: ", np.sort(flux_array)[-5:])

    #print('mjd_array.shape: ', mjd_array.shape)
    #print('gross_array.shape: ', gross_array.shape)
    #print('flux_array.shape: ', flux_array.shape)

    #print("mjd_diff: ", mjd_array - np.roll(mjd_array,1))
    if plotall:
        fig = plt.figure(figsize= (20,9))
        ax = fig.add_subplot(1,1,1)
        ax.scatter(mjd_array, gross_array)
        ax.set_xlabel('mjd')
        ax.set_ylabel('gross')
        ax.set_title(target_dir[:-1] + ' step= ' + str(stepsize))
        fig.savefig(target_dir[:-1]+ '_gross_step' + str(stepsize)+'_wlim'+ str(wlim[0])+',' + str(wlim[1])+ '.pdf')
        plt.show()
        #plt.clf()

        #fig2 = plt.figure(figsize= (20,9))
        #ax = fig2.add_subplot(1,1,1)

        #ax.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
        #ax.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
        #ax.scatter(mjd_array, flux_array)
        #ax.scatter(mjd_array, pre_flux, color = 'r', alpha = 0.5)
        #ax.set_xlabel('mjd')
        #ax.set_ylabel('flux')
        ##plt.xlim(56921.256, 56921.258)
        ##ax.set_ylim(0,0.005)
        #plt.legend(['y=1','y = 0','normalized flux', 'not normed'])
        #ax.set_title(target_dir[:-1] + ' step= ' + str(stepsize))
        #fig2.savefig(target_dir[:-1]+ '_flux_step' + str(stepsize)+ '_wlim'+ str(wlim[0])+',' + str(wlim[1])+ '.pdf')
        #plt.show()

        plt.title('normed')
        plt.hist(flux_array, bins = 200)
        plt.show()

        plt.title('not normed')
        plt.hist(gross_array, bins= 20)
        plt.show()
    else:
        pass
    #see if this is being run exernally to see if we should be gathering the lightcurves somewhere
    
    #time_sec= (mjd_array- mjd_array[0])*second_per_mjd #This needs to be changed for the times in seconds to be in terms of bmjd_tdb
    sub_dir,sub_file = assembled_files[0].split('/')
    
    os.chdir(sub_dir+'/') #change into the directory that contains the files so we can look at the header
    hdu= fits.open(sub_file) #look at the header
    os.chdir('../') #change back to the original directory we were in
    target_ra= hdu[0].header['RA_TARG'] #header in degrees
    target_dec= hdu[0].header['DEC_TARG'] #headers in degrees
    print("starting time corrections")
    target_coord = coord.SkyCoord(target_ra, target_dec, unit= (u.deg, u.deg), frame= 'icrs')
    time_mjd= Time(mjd_array, format= 'mjd', scale= 'utc', location= (0.*u.deg, 0.* u.deg)) #assumes location at lon = 0 , lat =0, elev= sea level
    ltt_bary= time_mjd.light_travel_time(target_coord) #light travel time for the target for the barycentric correction
    print("light travel time calculated")
    print("max time correction: ", np.max(np.abs(ltt_bary)))
    bmjd_array= (time_mjd.tdb + ltt_bary).tdb.mjd #barycentric correction
    time_sec= (bmjd_array- bmjd_array[0])*second_per_mjd #This should correctly output the seconds times in the BMJD_tdb version
    #textheader= 'mjd\tgross\tflux\ttime(s)(bmjd_tdb)\tbmjd_tdb'
    textheader= 'mjd\tgross\tflux\ttime(s)(bmjd_tdb)\tbmjd_tdb\tbackground\tcounts\terror\tflux_error'
    comment_string= 'mask_deg='+str(config.mask_deg) #20190610
    textheader=comment_string+'\n' + textheader #20190610
    #textcomment= 'step = ' + str(stepsize)
    if __name__ != "__main__":
        dest_dir = target_dir[:-1] + "_grid_lightcurves/"
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        os.chdir(dest_dir)
    print(os.getcwd())
    out_array= np.append([mjd_array], [gross_array], axis = 0)
    out_array = np.append(out_array, [flux_array], axis = 0)
    out_array = np.append(out_array, [time_sec], axis=0)
    out_array= np.append(out_array, [bmjd_array], axis=0)
    out_array= np.append(out_array, [bkg_array], axis=0)
    out_array= np.append(out_array, [count_array], axis=0)
    out_array= np.append(out_array, [error_array], axis=0)
    out_array= np.append(out_array, [norm_ferror_array], axis=0)

    print(out_array.shape)
    out_array= out_array.T
    print(out_array.shape)
    print("writing textfile")
    #np.savetxt(target_dir[:-1]+'_lightcurve_step' + str(stepsize)+'_wlim'+ str(wlim[0])+',' + str(wlim[1])+'.txt', out_array, delimiter = '\t', header = textheader)
    np.savetxt(target_dir[:-1]+'_lightcurve_step' + str(stepsize)+'_wlim'+ str(wlim[0])+',' + str(wlim[1])+'.txt', out_array, delimiter = '\t', header = textheader) #20190610
    if __name__ != '__main__':
        os.chdir('../')
    print(os.getcwd())
    return "done"


def do_lightcurve(thisfile, stepsize):
    astrotable, astrometa=lc.stis.extract(thisfile,step=stepsize)
    mjd_array= astrotable['mjd']
    gross_array= astrotable['gross']
    flux_table=  astrotable["flux"] #for complete dataset mean normalization
    bkg_array=astrotable['background']
    count_array=astrotable['gross']-astrotable['background']
    error_array=np.sqrt(astrotable['gross']+astrotable['background'])
    ferror_array= np.sqrt(astrotable['gross'] + astrotable['background']) / (astrotable['gross'] - astrotable['background']) * astrotable['flux']
    
    pre_flux = np.copy(flux_table)
    ferror_array= np.copy(ferror_array/flux_table) #normalizing the error to the flux values to whatever units flux happens to be in at a given time
    #norm_ferror_array= np.copy(ferror_array*(flux_table/np.nanmean(flux_table))) #making the flux error array be in scaled units to the normalized flux that's about to be calculated
    #it had to be done without subtracting 1 because then you'd have negative values for the errors inexplicably
    #this is also why it couldn't be left in arbitrary units because it would later be multiplied by the negative flux values
    #flux_array = np.copy( flux_table/np.nanmean(flux_table)-1.) #normalize flux array

    
    return mjd_array, gross_array,bkg_array,error_array,ferror_array, count_array, flux_table




def make_easy_lightcurve(filename='*_corrtag.fits', stepsize=5,target_dir=os.getcwd(),plotall=False, name_append=''):
    mjd_array = np.array([])
    gross_array = np.array([])
    flux_table = np.array([])
    bkg_array=np.array([])
    count_array=np.array([])
    error_array=np.array([])
    ferror_array=np.array([])
    
    sorted_filename=sorted(glob(filename))
    

    if len(sorted_filename)==0:
        tempfilename=filename.replace('corrtag','tag')
        sorted_filename=sorted(glob(tempfilename))
        print('Creating corrtag stis files.')
        for thisfile in sorted_filename:
            lc.stis.stis_corrtag(thisfile,clean=False)#clean=False means it won't delete the *tag.fits files when it's done with them, which it wants to do because memory saving I suppose, but then I have to re-retrieve it.
        sorted_filename=sorted(glob(filename))
    else:
        pass
    header=fits.getheader(sorted_filename[0])
    for thisfile in sorted_filename:
        #sorted(glob(filename))
        mjd_array_new, gross_array_new,bkg_array_new,error_array_new,ferror_array_new, count_array_new, flux_table_new = do_lightcurve(thisfile,stepsize)
        #astrotable, astrometa=lc.stis.extract(thisfile,step=stepsize)
        mjd_array= np.append(mjd_array, mjd_array_new)
        gross_array= np.append(gross_array, gross_array_new)
        flux_table= np.append(flux_table, flux_table_new) #for complete dataset mean normalization
        bkg_array=np.append(bkg_array, bkg_array_new)
        count_array=np.append(count_array, count_array_new)
        error_array=np.append(error_array, error_array_new)
        ferror_array= np.append(ferror_array, ferror_array_new)
        
    pre_flux = np.copy(flux_table)
    ferror_array= np.copy(ferror_array/flux_table) #normalizing the error to the flux values to whatever units flux happens to be in at a given time
    norm_ferror_array= np.copy(ferror_array*(flux_table/np.nanmean(flux_table))) #making the flux error array be in scaled units to the normalized flux that's about to be calculated
    #it had to be done without subtracting 1 because then you'd have negative values for the errors inexplicably
    #this is also why it couldn't be left in arbitrary units because it would later be multiplied by the negative flux values
    flux_array = np.copy( flux_table/np.nanmean(flux_table)-1.) #normalize flux array
    
    
    print("starting time corrections")
    target_coord = coord.SkyCoord(header['RA_TARG'], header['DEC_TARG'], unit= (u.deg, u.deg), frame= 'icrs')
    time_mjd= Time(mjd_array, format= 'mjd', scale= 'utc', location= (0.*u.deg, 0.* u.deg)) #assumes location at lon = 0 , lat =0, elev= sea level
    ltt_bary= time_mjd.light_travel_time(target_coord) #light travel time for the target for the barycentric correction
    print("light travel time calculated")
    print("max time correction: ", np.max(np.abs(ltt_bary)))
    bmjd_array= (time_mjd.tdb + ltt_bary).tdb.mjd #barycentric correction
    time_sec= (bmjd_array- bmjd_array[0])*second_per_mjd #This should correctly output the seconds times in the BMJD_tdb version
    #textheader= 'mjd\tgross\tflux\ttime(s)(bmjd_tdb)\tbmjd_tdb'
    
    out_array= np.append([mjd_array], [gross_array], axis = 0)
    out_array = np.append(out_array, [flux_array], axis = 0)
    out_array = np.append(out_array, [time_sec], axis=0)
    out_array= np.append(out_array, [bmjd_array], axis=0)
    out_array= np.append(out_array, [bkg_array], axis=0)
    out_array= np.append(out_array, [count_array], axis=0)
    out_array= np.append(out_array, [error_array], axis=0)
    out_array= np.append(out_array, [norm_ferror_array], axis=0)

    print(out_array.shape)
    out_array= out_array.T
    print(out_array.shape)
    print("writing textfile")
    
    
    
    
    textheader= 'mjd\tgross\tflux\ttime(s)(bmjd_tdb)\tbmjd_tdb\tbackground\tcounts\terror\tflux_error'
    comment_string= 'mask_deg='+str(config.mask_deg) #20190610
    textheader=comment_string+'\n' + textheader #20190610
    np.savetxt(target_dir+'_lightcurve_step' + str(stepsize)+name_append+'.txt', out_array, delimiter = '\t', header = textheader) #20190610
    
    if plotall:
        fig = plt.figure(figsize= (20,9))
        ax = fig.add_subplot(1,1,1)
        ax.scatter(mjd_array, gross_array)
        ax.set_xlabel('mjd')
        ax.set_ylabel('gross')
        ax.set_title(target_dir[:] + ' step= ' + str(stepsize))
        fig.savefig(target_dir+ '_gross_step' + str(stepsize)+ '.pdf')
        plt.show()
        
        fig = plt.figure(figsize= (20,9))
        ax = fig.add_subplot(1,1,1)
        ax.scatter(mjd_array, count_array)
        ax.set_xlabel('mjd')
        ax.set_ylabel('Net')
        ax.set_title(target_dir[:] + ' step= ' + str(stepsize))
        fig.savefig(target_dir[:]+ '_net_step' + str(stepsize)+ '.pdf')
        plt.show()
        
        fig = plt.figure(figsize= (20,9))
        ax = fig.add_subplot(1,1,1)
        ax.scatter(mjd_array, flux_array)
        ax.set_xlabel('mjd')
        ax.set_ylabel('Flux')
        ax.set_title(target_dir[:] + ' step= ' + str(stepsize))
        fig.savefig(target_dir[:]+ '_flux_step' + str(stepsize)+ '.pdf')
        plt.show()
        #plt.clf()

        #fig2 = plt.figure(figsize= (20,9))
        #ax = fig2.add_subplot(1,1,1)

        #ax.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
        #ax.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
        #ax.scatter(mjd_array, flux_array)
        #ax.scatter(mjd_array, pre_flux, color = 'r', alpha = 0.5)
        #ax.set_xlabel('mjd')
        #ax.set_ylabel('flux')
        ##plt.xlim(56921.256, 56921.258)
        ##ax.set_ylim(0,0.005)
        #plt.legend(['y=1','y = 0','normalized flux', 'not normed'])
        #ax.set_title(target_dir[:-1] + ' step= ' + str(stepsize))
        #fig2.savefig(target_dir[:-1]+ '_flux_step' + str(stepsize)+ '_wlim'+ str(wlim[0])+',' + str(wlim[1])+ '.pdf')
        #plt.show()

        #plt.title('normed')
        #plt.hist(flux_array, bins = 200)
        #plt.show()

        #plt.title('not normed')
        #plt.hist(gross_array, bins= 20)
        #plt.show()
    else:
        pass
    
    return

def make_lightcurves_by_wavelengths():
    sorted_filenames=sorted(glob('*corrtag.fits'))
    cenwave_list=[]
    for filename in sorted_filenames:
        header=fits.getheader(filename)
        cenwave_list.append(header['cenwave'])
        
    print('cenwave_list',cenwave_list)
    unique_cenwaves=np.unique(cenwave_list)
    print('unique_cenwaves',unique_cenwaves)
    filename_array=np.array(sorted_filenames)
    
    
    
    return

#Check if this is being executed alone
if __name__ == '__main__':
    
    if len(sys.argv)>1:
        print('User input provided')
    else:
        print('No user input provided')
        
    make_easy_lightcurve(plotall=True)
    #target_dir = sys.argv[1] + '/'
    #stepsize = float(sys.argv[2]) #still not sure what that actually means. 
    ##wlim= [1130, 1900] 
    #lowerlim= int(sys.argv[3])
    #upperlim= int(sys.argv[4])
    #wlim= [lowerlim, upperlim]
    #outstring = make_lightcurve(target_dir, stepsize, wlim)

    
    #plt.plot(astrotable['mjd'], astrotable['gross'])
    #plt.title(item)
    #plt.show()

    #plt.plot(astrotable['mjd'], astrotable['flux'])
    #plt.title(item)
    #plt.show()
        
        
#for item in glob(target_dir + '*corrtag*.fits'):
    #astrotable, astrometa = lc.cos.extract(item, step = 2)
    
    
    #plt.plot(astrotable['mjd'], astrotable['gross'])
    #plt.title(item)
    #plt.show()

    #plt.plot(astrotable['mjd'], astrotable['flux'])
    #plt.title(item)
    #plt.show()
#print(astrotable)
##print(astrotable.shape)


#for j in astrotable:
    #print(j)
    #print(astrotable[j])
    
