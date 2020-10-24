import sys
sys.path.insert(0, './functions')

import ngl

import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from pattern_cor import *


nc_f = './testout.nc'  # Your filename
nc_fid = Dataset(nc_f, 'r')  # Dataset is the class behavior to open the file
                             # and create an instance of the ncCDF4 class

# Extract data from NetCDF file
lats = nc_fid.variables['lat'][:]  # extract/copy the data
lons = nc_fid.variables['lon'][:]
trackdens = nc_fid.variables['fulldens'][:]  # shape is time, lat, lon as shown above

trackdensdims = trackdens.shape
nmods = trackdensdims[0]

pi = 3.141592653589793
deg2rad = pi / 180.
denslatwgt = np.cos(deg2rad*lats)

for ii in range(nmods):
  tmp = pattern_cor(trackdens[0,:,:],trackdens[ii,:,:],denslatwgt,0)
  print(tmp)
  
for ii in range(nmods):
  #print("avg "+str(np.nanmean(trackdens[0,:,:]))+" "+str(np.nanmean(trackdens[ii,:,:])))
  tmp = wgt_arearmse2(trackdens[0,:,:],trackdens[ii,:,:],denslatwgt,0)
  print(tmp)
  
for ii in range(nmods):
  #print("avg "+str(np.nanmean(trackdens[0,:,:]))+" "+str(np.nanmean(trackdens[ii,:,:])))
  ratio = taylor_stats(trackdens[ii,:,:],trackdens[0,:,:],denslatwgt,0)
  print(ratio)
