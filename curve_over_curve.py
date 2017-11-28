"""
written by Ben Kaiser

Takes the text file inputs of two light curves produced by try_lightcurve.py (or run_step_grid.py, which calls the previous program)
and divides the first by the second after first adding one to the zeroed flux values.

The input lightcurve files must have the same timestep and time bins at the moment. No interpolation is allowed. The timestamps will
be taken from the first input file.

"""

__author__ ='Ben Kaiser'


import numpy as np
import os
from glob import glob
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from astropy.io import fits
import sys
#from astropy.time import Time


dest_dir = 'divided_lightcurves/'
if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

def make_output_name(path1, path2):
    """
    takes the input names that include a pointing to a sub-directory and extracts the relevant information to produce the output filename that we'll use
    """
    def extract_wlims(inpath):
       inname= inpath.split('/')[1]
       prefix= inname.split('.')[0]
       wlim_with_word= prefix.split('_')[-1]
       wlims= wlim_with_word.split('m')[-1]
       return wlims
    
    name1= path1.split('/')[-1]
    out_base= name1.split('_')[:-1]
    out_base= '_'.join(out_base)
    wlim1= extract_wlims(path1)
    wlim2= extract_wlims(path2)
    #print out_base
    out_name = out_base + '_wlim'+wlim1+ '_over' + wlim2 + '.txt'
    #print "out_name: ", out_name
    return out_name
    

def lightcurve_divide(path1, path2):
    all1= np.genfromtxt(path1, names= True)
    print np.genfromtxt(path1,skip_header = 0)[0]
    all2= np.genfromtxt(path2, names = True)
    #print all1
    bmjd_array = all1['bmjd_tdb']
    #time_s = all1['time(s)(bmjd_t']
    time_s= np.genfromtxt(path1).T[3] #I don't know why I can't successfully extract the header but I can't
    print ("WARNING!")
    print ("For whatever reason, np.genfromtxt can't recognize 'time(s)(bmjd_tdb)' as the header,")
    print ("so we're just going off indices. If you've changed any of the text file layouts, this could")
    print ("be bad. If the numbers printed below vaguely correspond to your timestep, you should")
    print ("be good.")
    lc1= all1['flux']+1
    lc2= all2['flux']+1
    div_curve = lc1/lc2-1
    #print all1.shape
    #print all1[0]
    #print all1['mjd']
    #print all1[:, -2]
    print time_s[0:5]
    #print bmjd_array
    out_name = make_output_name(path1,path2)
    out_name= dest_dir+out_name
    out_array= np.vstack(([time_s], [div_curve], [bmjd_array]))
    out_array = out_array.T
    print "saving ", out_name
    np.savetxt(out_name, out_array, delimiter= '\t', header= 'time(s)(bmjd_tdb)\tflux\tbmjd_tdb')
    #subtracted version
    new_out= np.vstack(([time_s], [lc1-lc2], [bmjd_array]))
    new_out = new_out.T
    out_name= out_name.replace('_over', '_sub')
    print "saving the subtracted version", out_name
    np.savetxt(out_name,new_out, delimiter= '\t', header= 'time(s)(bmjd_tdb)\tflux\tbmjd_tdb')










###################################################
if __name__ == '__main__':
    path1 = sys.argv[1]
    path2= sys.argv[2]
    #make_output_name(path1, path2)
    lightcurve_divide(path1,path2)
