import numpy as np
import sys
import matplotlib.pyplot as plt
import config
import os
from glob import glob
from astropy.io import fits


target_path = sys.argv[1]+'/' #the input should be the target name just like all of the other times this has been implemented in this repository
#path to the model files
teff= 12000
logg= 8.0 #model values to pull (should be changed to be read from somewhere in the future)
lyman_alpha = [1214,1217] #probably should put this into the config file as the masking for the spectroscopic fitting (but really this goes back to the need for individual text files by target for the various maskings that are required)
#scale factor to be just typed in for the moment but should be replaced by a function at some point that calculates it for you
target_x1ds= glob(target_path+'*x1dsum.fits')
target_file= target_x1ds[0] #I only want the first one for now; when I eventually make this include all the visits for a given target I'll introduce some sort of function for combining them or whatever
model_spec= ''
#wavelength range to consider for the scaling

target_spec

#function to trim the model plot to limit it to the range of the target's spectrum

#At some point I'll have to implement masking
