import numpy as np
from nexusformat.nexus import *
from matplotlib import pyplot as plt

def load_nxs_data_snapshot(file_path, snapshot=0):
    '''
    Reads .nxs file from a given file path. Returns a 2d numpy array
    of a single snapshot. Uses the nexusformat interface:
    https://nexpy.github.io/nexpy/pythonshell.html
    '''
    a=nxload(file_path)
    data = np.array(a.entry.instrument.detector.data)[snapshot,:,:]
    return data

def load_nxs_data_series(file_path):
    '''
    Reads .nxs file from a given file path. Returns a 3d numpy array
    of a 2x2 data time series. Uses the nexusformat interface:
    https://nexpy.github.io/nexpy/pythonshell.html
    '''
    b = nxload(file_path) 
    data = b.NXentry[0].instrument.detector.data 
    return data

def load_nxs_metadata(file_path):
    '''
    Reads .nxs metada from a given file path. Returns a 2d numpy array
    of a data series. Uses the nexusformat interface:
    https://nexpy.github.io/nexpy/pythonshell.html
    '''
    a=nxload(file_path) 
    print a.tree 
    #count_time = a.entry.instrument.detector.count_time
    #shutter_time = a.entry.instrument.detector.collection.shutter_time
    #px = a.entry.instrument.detector.count_time
    #py = a.entry.instrument.detector.count_time
    #add here more metadata

    #metadata = np.array((count_time,shutter_time,px,py))
    return #metadata

def format_filepath(base_dir,data_dir,dataname,detector_name,start_number,end_number):
    '''
    Assistive function to contsruct the file pathe in the right format
    '''
    filesep = '/'
    filenum = '%.5d-%.5d'%(int(start_number),int(end_number))
    dataext = '.nxs'
    filename = dataname+'_'+filenum+dataext
    filepath = base_dir+data_dir+dataname+filesep+ detector_name+filesep+filename

    return filepath

