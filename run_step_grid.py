import numpy as np
from glob import glob
import os
import sys
import config
from astropy.io import fits
import try_lightcurve as tlc


"""
calls try_lightcurve.py and then runs it over a grid of step sizes, and potentially will iterate over the different target directories.

"""


#step_grid= np.arange(1,121, 1) #by going one above the desired value, it is included in the range. For example to include 120 make it 121.
#step_grid= np.arange(15,121, 15) #by going one above the desired value, it is included in the range. For example to include 120 make it 121.
#limit_list= [[1130,1430],[1130,1800],[1130,1850],[1130,1900]]
step_grid= np.arange(5,31,5)
app_array= np.arange(45,61,15)
low_steps= np.arange(1,5,1)
step_grid= np.append(low_steps, step_grid)
step_grid= np.append(step_grid, app_array)
step_grid= np.append(step_grid, [90,120,180,240,300])
def get_wave_limits(target_dir):
    #open a fits file and get the central wavelength and grating headers
    #look through the different if statements to find the appropriate grating
    #then look in the if-statements under that category for the correct wave_limit_list
    #This could probably be rewritten to branch through a dictionary of dictionaries of lists... but whatever
    checked_file = glob(target_dir+ '*corrtag*.fits')[0]
    hdu = fits.open(checked_file)
    grating= hdu[0].header['OPT_ELEM']
    center_wave= hdu[0].header['CENWAVE']
    if grating == "G130M":
        if center_wave == 1291:
            list_limits= config.G130M_1291A_wave_list
        elif center_wave == 1300:
            list_limits = config.G130M_1300A_wave_list
        elif center_wave == 1096:
            list_limits= config.G130M_1096A_wave_list
    elif grating == "G140L":
        if center_wave == 1105:
            list_limits = config.G140L_1105A_wave_list
    elif grating == "G160M":
        if center_wave == 1600:
            list_limits = config.G160M_1600A_wave_list
    return list_limits


def run_step_grid(target_dir):
    list_limits = get_wave_limits(target_dir)
    for wave_limit in list_limits:
        #try:
            #outstring = tlc.make_lightcurve(target_dir, 1, wave_limit, plotall = False) #one second binning before running the other steps in the grid
        #except IndexError as error:
            #print "\n******************************"
            #print error
            #print "******************************\n"
        for stepsize in step_grid:
            try:
                outstring = tlc.make_lightcurve(target_dir, stepsize, wave_limit, plotall = False)
            except IndexError as error:
                print "\n******************************"
                print error
                print "******************************\n"
    return ''


if __name__ == '__main__':
    target_dir = sys.argv[1] + '/'
    run_step_grid(target_dir)

