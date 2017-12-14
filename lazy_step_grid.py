import numpy as np
import sys
import run_step_grid as rsg 

list_targs = np.genfromtxt('targets_for_lc.txt', 'str')

for target in list_targs:
    rsg.run_step_grid(target) #this one takes the target directroy input with the slash tacked on... yeah I should have standardized this across the different things, but what do you want from me?
