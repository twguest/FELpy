#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.1.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import requests
import time
import json
import zipfile

import pandas as pd


def load_options():
    return pd.read_excel("../../data/FASTXPD_SA1_Options.xlsx")    
    
def get_FASTXPD_pulses(input_folder = 'XFEL_S1_04.96keV_12.0GeV_0020pC_SASE_U_BLI_2014-05-01_FAST',
                       number_xy = 50,
                       skip_slices = 10,
                       z_output_point = 20, 
                       time_begin = 0,
                       time_end = 35, 
                       from_run_number = 1,
                       to_run_number = 1, 
                       user_prefix = "xfel_pulse_test",
                       my_path = '/tmp',
                       my_email_address = 'twguest@students.latrobe.edu.au'):
    
     # where to download and extract result files
    
    url = 'https://in.xfel.eu/fastxpd/rest/requests'
    request_id = ''
    
    request_data = {
      'input_folder': input_folder ,
      'time_begin': time_begin,
      'time_end': time_end,
      'number_xy': number_xy,
      'skip_slices': skip_slices,
      'z_output_point': z_output_point,
      'from_run_number': from_run_number,
      'to_run_number': to_run_number,
      'email': my_email_address,
      'user_prefix': user_prefix
    }
    
    response = requests.post(url, data=request_data)
    
    if not response.ok:
      # If response code is not ok (200), print(the resulting http error code with description
      response.raise_for_status() 
    
    answer = json.loads(response.content)
    
    if answer['status'] > 0:
      # got an answer
      request_id = answer['request_id']
      status = 'Queued'
      print("Request id: %s" % request_id) 
    else:
      # got an error
      print(answer['messages'])  
    
    # poll request status
    if request_id:
      
      url_request = "%s/%s" % (url, request_id)
      
      while (status == 'Queued' or status == 'Processing...'):
        time.sleep(10) # poll status every 10 seconds
        response = requests.get(url_request)
        answer = json.loads(response.content)
        status = answer['current_status']
        print("Current status: %s" % status)
    
      if status == 'Completed':
        print("Request has been completed")
    
        # show available files
        files = answer['files']
        print("Available result files:")
        print(files)
        
        # download and extract HDF5 files if present
        if 'zipfile_hdf5' in files:
          zip_file_url = files['zipfile_hdf5']['link']
          download_path = "%s/%s" % (my_path, zip_file_url.split('/')[-1])
          r = requests.get(zip_file_url, stream=True, verify=False)
          print("Downloading %s in %s" % (zip_file_url, download_path))
          
          # downloading zip file
          with open(download_path, 'wb') as f:
            for chunk in r.iter_content(1024):
              f.write(chunk)
    
          print("File downloaded, extracting")    
    
          # extract files
          z = zipfile.ZipFile(download_path, allowZip64=True)
          z.extractall(my_path)
          
if __name__ == '__main__':
    #get_FASTXPD_pulses()
    
    options = load_options()

    for idx in range(options.shape[0]):
        input_folder = options['Pulse Name'][idx]
        time_end = options['End Time'][idx]
        z_output_point = options['Start Point'][idx]
        to_run_number  = options['Runs'][idx]
        
        name_str = input_folder.split("XFEL_")[1].split("_SASE")[0] + "_{}m".format(z_output_point)
        
        print(input_folder)
        print(z_output_point)
        
        get_FASTXPD_pulses(input_folder = input_folder,
                           number_xy = 50,
                           skip_slices = 10,
                           z_output_point = z_output_point, 
                           time_begin = 0,
                           time_end = time_end, 
                           from_run_number = 1,
                           to_run_number = 1, ###to_run_number, ###hardcoded for test
                           user_prefix = name_str,
                           my_path = "/gpfs/exfel/data/user/guestt/dCache/",
                           my_email_address = 'twguest@students.latrobe.edu.au')
        