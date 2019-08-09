"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-08-08


Take the extensive lightcurves and pare them down to bare-bones output versions per the request for sake of 
ease


This should be called inside the directory of the lightcurves you'd like to simplify, so you'll probably need to run

$ python ../simplify_lightcurves.py

because you'll most likely be one layer down assuming this is you, Ben.
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


import config


filename_string= '*lightcurve*.txt'


filenames = glob(filename_string)
wanted_headers= [
    'time_seconds',
    'flux',
    'flux_error'
    ]
delimiter= '\t'

output_header=delimiter.join(wanted_headers)
print('output_header', output_header)

def simplify_lightcurve(input_filename):
    print('simplifying', input_filename)
    all_array = np.genfromtxt(input_filename, names=True, skip_header=1)
    #all_array = np.genfromtxt(input_filename, skip_header=1)
    print(all_array.dtype.names)
    bmjd_tdb= all_array['bmjd_tdb']
    start_bmjd= bmjd_tdb[0]
    time_seconds=all_array['timesbmjd_tdb']
    flux= all_array['flux']
    flux_error=all_array['flux_error']
    output_array=np.vstack([time_seconds, flux, flux_error])
    output_array= output_array.T
    full_header= 'initial bmjd_tdb:' +str(start_bmjd)+'\n'+output_header
    #np.savetxt(input_filename, output_array, delimiter=delimiter, header= output_header)
    np.savetxt(input_filename, output_array, delimiter=delimiter, header= full_header)
    print(input_filename, 'simplified')
    return

    
    
for filename in filenames:
    simplify_lightcurve(filename)
    
    
