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


#target_dir = sys.argv[1] + '/'
#stepsize = int(sys.argv[2]) #still not sure what that actually means. Is it in seconds?
def make_lightcurve(target_dir, stepsize, plotall=True):
    """
    """
    print ""
    print "==============="
    print "stepsize: ", stepsize
    print "plotall=", plotall
    print "target_dir: ", target_dir
    wlim= [900, 3200] #3200 is the default max from what I could gather. The geocoronal emission lines are automatically excluded from the lightcurve
    band= '*' #should be '_a', '_b', '' for the combined band images, or  '*' to get all bands available

    mjd_array = np.array([])

    gross_array = np.array([])
    flux_table = np.array([])


    negative_flux_files= []
    negative_means= []
    if plotall:
        fig = plt.figure(figsize= (20,9))
        ax = fig.add_subplot(1,1,1)


    print os.getcwd()
    for item in glob(target_dir + '*corrtag' + band + '.fits'):
        astrotable, astrometa = lc.cos.extract(item, step=stepsize, wlim= wlim)
        #print 'astrometa: ', astrometa
        
        mjd_array= np.append(mjd_array, astrotable['mjd'])
        gross_array= np.append(gross_array, astrotable['gross'])
        flux_table= np.append(flux_table, astrotable['flux'])
        print "file: ", item
        print "np.nanmean(flux):", np.nanmean(astrotable['flux'])
        print '-------------'
        if np.nanmean(astrotable['flux'] < 0.):
            negative_flux_files.append(item)
            negative_means.append(np.nanmean(astrotable['flux']))
        
    print '==========='
    print "Negative flux files and means detected: "
    for filename, nanmean in zip(negative_flux_files, negative_means):
        print "File: ", filename, "| nanmean: ", nanmean
    print '==========='


    pre_flux = np.copy(flux_table)
    flux_array = np.copy( flux_table/np.nanmean(flux_table)) #normalize flux array

    print "np.nanmean(flux_array)" , np.nanmean(flux_array)
    print "np.nanmean(flux_table: " ,np.nanmean(flux_table)
    print "np.sum(flux_array-pre_flux): ", np.sum(flux_array- pre_flux)
    print "np.sum(flux_array + pre_flux): ", np.sum(flux_array + pre_flux)
    print ""
    print "np.nanmax(flux_array): ", np.nanmax(flux_array)
    print "mjd value for that: ", mjd_array[np.argmax(flux_array)]
    print "top 5 flux values: ", np.sort(flux_array)[-5:]

    #print 'mjd_array.shape: ', mjd_array.shape
    #print 'gross_array.shape: ', gross_array.shape
    #print 'flux_array.shape: ', flux_array.shape

    #print "mjd_diff: ", mjd_array - np.roll(mjd_array,1)
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

        fig2 = plt.figure(figsize= (20,9))
        ax = fig2.add_subplot(1,1,1)

        ax.axhline(y= 1, linestyle= '-', color = 'm', xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
        ax.axhline(y= 0,linestyle = '-', color = 'g' , xmin = 0, xmax = 100000, linewidth = 1, alpha = 0.2)
        ax.scatter(mjd_array, flux_array)
        ax.scatter(mjd_array, pre_flux, color = 'r', alpha = 0.5)
        ax.set_xlabel('mjd')
        ax.set_ylabel('flux')
        #plt.xlim(56921.256, 56921.258)
        #ax.set_ylim(0,0.005)
        plt.legend(['y=1','y = 0','normalized flux', 'not normed'])
        ax.set_title(target_dir[:-1] + ' step= ' + str(stepsize))
        fig2.savefig(target_dir[:-1]+ '_flux_step' + str(stepsize)+ '_wlim'+ str(wlim[0])+',' + str(wlim[1])+ '.pdf')
        plt.show()

        plt.title('normed')
        plt.hist(flux_array, bins = 70)
        plt.show()

        plt.title('not normed')
        plt.hist(pre_flux, bins= 70)
        plt.show()
    else:
        pass
    #see if this is being run exernally to see if we should be gathering the lightcurves somewhere
    if __name__ != "__main__":
        dest_dir = target_dir[:-1] + "_grid_lightcurves/"
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
        os.chdir(dest_dir)
    print os.getcwd()
    textheader= 'mjd\tgross\tflux'
    #textcomment= 'step = ' + str(stepsize)
    out_array= np.append([mjd_array], [gross_array], axis = 0)
    out_array = np.append(out_array, [flux_array], axis = 0)
    print out_array.shape
    out_array= out_array.T
    print out_array.shape
    print "writing textfile"
    np.savetxt(target_dir[:-1]+'_lightcurve_step' + str(stepsize)+'_wlim'+ str(wlim[0])+',' + str(wlim[1])+'.txt', out_array, delimiter = '\t', header = textheader)
    if __name__ != '__main__':
        os.chdir('../')
    print os.getcwd()
    return "done"

#Check if this is being executed alone
if __name__ == '__main__':
    target_dir = sys.argv[1] + '/'
    stepsize = int(sys.argv[2]) #still not sure what that actually means. Is it in seconds?
    outstring = make_lightcurve(target_dir, stepsize)

    
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
#print astrotable
##print astrotable.shape


#for j in astrotable:
    #print j
    #print astrotable[j]
    
