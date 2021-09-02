#!/usr/bin/env python
# Script to create a wrf profile netCDF file

from netCDF4 import Dataset
import numpy as     np
import sys
from scipy import interpolate

wrffilename = "wrfforcing.nc"

rootgrp = Dataset(wrffilename, "w", format="NETCDF4")
print(rootgrp.data_model)

# Set some variables
NaluDensity = 1.00
NaluTflux   = 0.01220961466456118644
t0          = 0.0
t1          = 100000

# --- Load the profile data ---
# download data from https://github.com/lawrenceccheung/IEATask29_precursor/blob/main/REF1/NaluWind_newBC.refT
velprof  = np.loadtxt('NaluWind_velocity.dat')
tempprof = np.loadtxt('NaluWind_temperature.dat')

# doublecheck that the z-profiles match
if np.linalg.norm(velprof[:,0]-tempprof[:,0])>1.0E-8:
    print("Z profiles do not match!")
    sys.exit(1)
else:
    zprof_nalu = velprof[:,0]
    uprof_nalu = velprof[:,1]
    vprof_nalu = velprof[:,2]
    Tprof_nalu = tempprof[:,1]

# Get the AMR-Wind grid
zamr = np.linspace(6,1914,160)
# Map from nalu grid to AMR-Wind grid
interpu = interpolate.interp1d(zprof_nalu, uprof_nalu)
interpv = interpolate.interp1d(zprof_nalu, vprof_nalu)
interpT = interpolate.interp1d(zprof_nalu, Tprof_nalu)

uprof   = interpu(zamr)
vprof   = interpv(zamr)
Tprof   = interpT(zamr)

print(len(uprof), len(vprof), len(Tprof))

# Create the heights
heights   = zamr
nheight   = rootgrp.createDimension("nheight", len(heights))
ncheights = rootgrp.createVariable("heights", "f8", ("nheight",))
ncheights[:] = heights

# Create the times
times     = np.array([t0, t1])
ntime     = rootgrp.createDimension("ntime", len(times))
nctimes   = rootgrp.createVariable("times", "f8", ("ntime",))
nctimes[:] = times

# Add momentum u
nc_momu     = rootgrp.createVariable("wrf_momentum_u", "f8", 
                                     ("ntime", "nheight",))
for i in range(len(times)):
    nc_momu[i,:] = NaluDensity*uprof

# Add momentum v
nc_momv     = rootgrp.createVariable("wrf_momentum_v", "f8", 
                                     ("ntime", "nheight",))
for i in range(len(times)):
    nc_momv[i,:] = NaluDensity*vprof

# Add the temperature
nc_temp     = rootgrp.createVariable("wrf_temperature", "f8", 
                                     ("ntime", "nheight",))
for i in range(len(times)):
    nc_temp[i,:] = Tprof

# Add the temperature fluxes
tflux       = NaluTflux*np.ones(len(times))
nc_tflux    = rootgrp.createVariable("wrf_tflux", "f8", ("ntime",))
nc_tflux[:] = tflux

rootgrp.close()
