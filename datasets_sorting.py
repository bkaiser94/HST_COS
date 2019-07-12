"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-07-11


Should take the myriad datasets from '1retrieval/' and then put them into new directories that should follow the 
standard naming convention. The new directories will have to be generated from the .txt file that contains all of 
the datasets and affiliated information about the observations.



"""

from __future__ import print_function

import numpy as np
import os
from glob import glob
from astropy.table import Table, Column
import sys
import csv



import config

home_dir= os.getenv('HOME')

info_file= 'SNAP_done_smaller.csv'
massive_dir= '2retrieval/' #directory that contains the datasets that need to be sorted into the other directories
targeting_file = 'targets_for_lc.txt' #filename for the collected directories to be evaluated by lazy_step_grid.py
#all_array=np.genfromtxt(info_file, delimiter='/t', names=True).T
#all_table= Table.read(info_file, format='ascii.csv')
#all_array= np.genfromtxt(info_file, delimiter=',', names=True, dtype='S32')
#print(all_array['Grating'])



#############################
replacement_list=[
    ['+', 'p'],
    ['-','m']
    ]


#################################
def make_dir_name(input_row):
    wd_name = str(row['Target Name'])
    grating= str(row['Grating'])
    print(wd_name)
    print('replacing characters')
    for rep in replacement_list:
        wd_name=wd_name.replace(rep[0], rep[1])
        print('cleaned name:', wd_name)
    dir_name= wd_name+'_'+grating+'/'
    print('dir_name:', dir_name)
    return dir_name

def check_directory(input_row):
    dest_dir= make_dir_name(input_row)
    if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
    else:
        pass
    return dest_dir

def get_files_associated(input_row):
    file_association= input_row['File']
    core_assoc= file_association[:6]
    print('core_assoc', core_assoc)
    search_string= home_dir+massive_dir+core_assoc+'*.fits'
    print('search_string', search_string)
    coll_files= glob(search_string)
    print('coll_files', coll_files)
    return coll_files

def move_fileset(input_row):
    coll_files= get_files_associated(input_row)
    dest_dir= check_directory(row)
    for origin in coll_files:
        parts = origin.split('/')
        dest_filename= parts[-1]
        print('dest_filename', dest_filename)
        full_dest= home_dir+dest_dir+dest_filename
        print('moving', origin, ' to ', full_dest)
    return

###############################

counter=0
collected_lists= []
dirs_for_eval=[]
with open(info_file, 'rb') as csvfile:
    reader= csv.reader(csvfile, delimiter=',')
    for row in reader:
        if counter==0:
            names_list= row
        else:
            collected_lists.append(row)
        counter+=1
        print(row)
        
print(collected_lists)
print(len(collected_lists))

info_array= np.array(collected_lists)
print('shape', info_array.shape)
print('len(names_list)', len(names_list))
names_tuple= tuple(names_list)
print('names_tuple', names_tuple)
#info_array=info_array.T

info_table= Table(info_array, names=names_tuple)

info_table.pprint()

for row in info_table:
    #make_dir_name(row)
    dest_dir= check_directory(row)
    print('output dest_dir', dest_dir)
    dirs_for_eval.append([dest_dir])
    #get_files_associated(input_row)
    move_fileset(row)
    
    
print('saving', targeting_file)
np.savetxt(targeting_file, dirs_for_eval)
print(targeting_file, 'saved.')


