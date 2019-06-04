"""
Test the different methods of producing the 'correct' bin edges for binning the wavelengths together

"""

from __future__ import print_function


import numpy as np
import matplotlib.pyplot as plt
import config
import sys
from glob import glob
from astropy.io import fits
import run_step_grid as rsg

target_dir = sys.argv[1] + '/'
file_target = '*corrtag_a.fits'
get_thing = glob(target_dir+ file_target)[0]
hdu = fits.open(get_thing)
original_timestamps= hdu['events'].data['time']
start= 0
end = np.max(original_timestamps)


def get_step(step):
    instep = step
    step = step - (step%config.cos_refresh_rate) #method of calculating the stepsize that conforms to the native time bins
    #print("input step: ", instep)
    #print("calculated step: ", step)
    #print("difference: ", (instep-step))
    return step

def calc_arange1(step):
    return np.arange(start, end, step)

def calc_arange2(step):
    return np.float32(np.arange(start,end, step))

def calc_arange3(step):
    return np.arange(start, end, step, dtype = 'float32')

def calc_arange4(step):
    return np.arange(np.float32(start), np.float32(end), np.float32(step), dtype= 'float32')

def loop_range(step):
    compiled_list = []
    value= start
    while value < end:
        compiled_list.append(value)
        value+= step
    return np.array(compiled_list, dtype= 'float32')
        
def loop_rangef32(step):
    """
    this one should represent the actual values that I get out... maybe
    """
    compiled_list= []
    value= np.float32(start)
    step= np.float32(step)
    while value < np.float32(end):
        compiled_list.append(value)
        value+= step
    return np.array(compiled_list, dtype = 'float32')


def run_tests(instep):
    step= get_step(instep)
    arange1= calc_arange1(step)
    arange2= calc_arange2(step)
    arange3= calc_arange3(step)
    arange4 = calc_arange4(step)
    range_loop = loop_range(instep)
    range_loopf32= loop_rangef32(instep)
    plt.plot(arange1-arange2)
    plt.title('arange1-arange2')
    plt.show()
    plt.plot(arange1-arange3)
    plt.title('arange1-arange3')
    plt.show()
    plt.plot(arange2-arange3)
    plt.title('arange2-arange3')
    plt.show()
    plt.plot(arange4- arange2)
    plt.title('arange4-arange2')
    plt.show()
    plt.plot(range_loop[:range_loopf32.shape[0]]- range_loopf32)
    plt.title('range_loop-range_loopf32')
    plt.show()
    plt.plot(arange4[:range_loopf32.shape[0]]-range_loopf32)
    plt.title('arange4-range_loopf32')
    plt.show()
    
def hist_tests(instep):
    step = get_step(instep)
    arange1= calc_arange1(step)
    arange2= calc_arange2(step)
    arange3= calc_arange3(step)
    arange4 = calc_arange4(step)
    range_loop = loop_range(instep)
    range_loopf32= loop_rangef32(instep)
    hist_arange1= np.histogram(original_timestamps, arange1)[0]
    hist_arange2= np.histogram(original_timestamps, arange2)[0]
    hist_arange3= np.histogram(original_timestamps, arange3)[0]
    hist_arange4= np.histogram(original_timestamps, arange4)[0]
    hist_range_loop = np.histogram(original_timestamps, range_loop)[0]
    hist_range_loopf32= np.histogram(original_timestamps, range_loopf32)[0]
    #plt.plot(hist_arange1-hist_arange2)
    #plt.title('arange1-arange2')
    #plt.show()
    #plt.plot(hist_arange1-hist_arange3)
    #plt.title('arange1-arange3')
    #plt.show()
    #plt.plot(hist_arange2-hist_arange3)
    #plt.title('arange2-arange3')
    #plt.show()
    #plt.plot(hist_arange4- hist_arange2)
    #plt.title('arange4-arange2')
    #plt.show()
    #plt.plot(hist_range_loop[:hist_range_loopf32.shape[0]]- hist_range_loopf32)
    #plt.title('range_loop-range_loopf32')
    #plt.show()
    #plt.plot(hist_arange4[:hist_range_loopf32.shape[0]]-hist_range_loopf32)
    #plt.title('arange4-range_loopf32')
    #plt.show()
    ####3
    plt.scatter( np.arange(0,hist_arange1.shape[0],1),hist_arange1)
    plt.title('arange1')
    plt.show()
    plt.scatter(np.arange(0,hist_arange3.shape[0],1),hist_arange3)
    plt.title('arange3')
    plt.show()
    plt.scatter(np.arange(0,hist_arange2.shape[0],1),hist_arange2)
    print("poisson noise expected: ", np.sqrt(np.mean(hist_arange2)))
    print("standard deviation calculated: ", np.std(hist_arange2))
    plt.title('arange2')
    plt.show()
    plt.scatter(np.arange(0,hist_arange4.shape[0],1),hist_arange4)
    plt.title('arange4')
    plt.show()
    plt.scatter(np.arange(0,hist_range_loop.shape[0],1),hist_range_loop)
    plt.title('range_loop')
    plt.show()
    plt.scatter(np.arange(0,hist_range_loopf32.shape[0],1),hist_range_loopf32)
    plt.title('range_loopf32')
    plt.show()

def test_noise(instep):
    step = get_step(instep)
    arange2= calc_arange2(step)
    hist_arange2= np.histogram(original_timestamps, arange2)[0]
    sigma= np.std(hist_arange2)/np.mean(hist_arange2)
    poisson= np.sqrt(np.mean(hist_arange2))/np.mean(hist_arange2)
    return sigma, poisson

def test_noise1(instep):
    #step = get_step(instep)
    step= instep
    arange2= calc_arange1(step)
    hist_arange2= np.histogram(original_timestamps, arange2)[0]
    sigma= np.std(hist_arange2)/np.mean(hist_arange2)
    poisson= np.sqrt(np.mean(hist_arange2))/np.mean(hist_arange2)
    return sigma, poisson

def relative_error_calc(times, test_num= 2):
    all_sigma= np.array([])
    all_poisson= np.array([])
    all_times= np.array([])
    for time in times:
    #for multiple in range(1,20):
        if test_num == 2:
            sigma, poisson = test_noise(time)
        elif test_num == 1:
            sigma, poisson = test_noise1(time)
        all_times= np.append(all_times, time)
        all_poisson= np.append(all_poisson, poisson)
        all_sigma= np.append(all_sigma, sigma)
    combined_array = np.vstack([all_times, all_poisson,all_sigma])
    combined_array = combined_array.T
    textheader= 'time\tpoisson\tstd_dev'
    np.savetxt('arange'+str(test_num)+'_relative_errors.txt', combined_array, delimiter= '\t', header = textheader)
    for row in combined_array:
        print(row)

    plt.scatter(all_times, all_poisson, label= 'poisson')
    plt.scatter(all_times, all_sigma, label = 'std deviation', color = 'r')
    plt.legend(numpoints=1, fontsize=14, loc='best')
    plt.title('arange'+ str(test_num)+'_'+get_thing + '_noise_estimates')
    plt.show()
    
times= rsg.step_grid
relative_error_calc(times, test_num=1)
relative_error_calc(times, test_num =2)

    

#run_tests(config.cos_refresh_rate)
#run_tests(config.cos_refresh_rate)
#hist_tests(config.cos_refresh_rate)
#hist_tests(0.1)
#print(1)
#hist_tests(1)
#print(1.1)
#hist_tests(1.1)
#print(0.64)
#hist_tests(0.64)

#for value in [config.cos_refresh_rate, np.float32(config.cos_refresh_rate), 1]:
    #print(np.min(np.abs(original_timestamps - get_step(value))))
