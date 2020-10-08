#!/usr/bin/env python
#
#


import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

def timeaverage(t, dat, t1, t2):
    tfiltered   = t[(t>=t1)&(t<=t2)]
    datfiltered = dat[(t>=t1)&(t<=t2),:]
    Nvars  = len(dat[0,:])
    tstart = tfiltered[0]
    tend   = tfiltered[-1]
    avgdat = np.zeros(Nvars)
    #print('Nvars = %i tstart = %f tend = %f'%(Nvars, tstart, tend))
    #print(len(tfiltered))
    #print(len(datfiltered))
    for i in range(len(tfiltered)-1):
        dt     = tfiltered[i+1] - tfiltered[i]
        avgdat = avgdat + 0.5*dt*(datfiltered[i+1,:] + datfiltered[i,:])
    return avgdat/(tend-tstart)

stdvars = ['u',         'v',      'w',        'theta', 
           u"u'u'_r",  u"u'v'_r", u"u'w'_r", 
           u"v'v'_r",  u"v'w'_r", u"w'w'_r"]

def loadData(filename, varslist=stdvars, group='mean_profiles', avgt=[]):
    alldat={}
    with Dataset(filename) as d:
        #print(d['mean_profiles'].variables)
        t = d.variables['time'][:]
        alldat['t'] = t
        alldat['z'] = d['mean_profiles'].variables['h'][:]
        for var in varslist:
            print('Loading '+var)
            x = d[group].variables[var][:,:]
            if len(avgt)>=2:
                t1 = avgt[0]
                t2 = avgt[1]
                alldat[var] = timeaverage(t, x, t1, t2)
            else:
                alldat[var] = x
    return alldat

def saveAsText(alldat, varslist, filename, header):
    # Make the list of vars
    writedat = alldat['z']
    for var in varslist:
        writedat = np.vstack((writedat, alldat[var]))
    np.savetxt(filename, writedat.transpose(), header=header)
    return

# def pullData(filename):
#     alldat={}
#     with Dataset(filename) as d:
#         t = d.variables['time'][:]
#         alldat['t'] = t
#         alldat['z'] = d['mean_profiles'].variables['h'][:]
#         U    = d['mean_profiles'].variables['u'][:,:]
#         V    = d['mean_profiles'].variables['v'][:,:]
#         Temp = d['mean_profiles'].variables['theta'][:,:]
#         avgU = timeaverage(t, U, 15000.0, 20000.0)
#         avgV = timeaverage(t, V, 15000.0, 20000.0)
#         avgT = timeaverage(t, Temp, 15000.0, 20000.0)
#         alldat['U'] = avgU
#         alldat['V'] = avgV
#         alldat['T'] = avgT
#     return alldat

# avgdat=pullData('abl_statistics00000.nc')

# DanAeroZ  = [17,        28.5,      41,        57,        77,        90]
# DanAeroWS = [5.884,     5.973,     5.931,     6.128,     6.028,     6.088]
# plt.plot(DanAeroWS, DanAeroZ, 'bs')
# plt.plot(avgdat['U'], avgdat['z'])
# plt.plot(avgdat['V'], avgdat['z'])
# plt.plot(np.sqrt(avgdat['U']**2 + avgdat['V']**2), avgdat['z'])
# plt.xlim([0,10])
# plt.ylim([0,160])
# plt.show()

# plt.plot(avgdat['T'], avgdat['z'])
# plt.show()


# with Dataset('abl_statistics00000.nc') as d:
#     t    = d.variables['time'][:]
#     z    = d['mean_profiles'].variables['h'][:]
#     U    = d['mean_profiles'].variables['u'][:,:]
#     V    = d['mean_profiles'].variables['v'][:,:]
#     #W    = d['mean_profiles'].variables['w'][:,:]
#     #Temp = d['mean_profiles'].variables['theta'][:,:]
#     avgU = timeaverage(t, U, 15000.0, 20000.0)
#     avgV = timeaverage(t, V, 15000.0, 20000.0)
#     plt.plot(avgU, z)
#     plt.plot(avgV, z)
#     plt.plot(np.sqrt(avgU**2 + avgV**2), z)
#     plt.show()
#     #print(t)
#     #print(z)
