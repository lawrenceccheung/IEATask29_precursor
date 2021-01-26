import sys
import numpy as np
from netCDF4 import Dataset
import os.path
import argparse

"""You can use this script via the command line or in your own python
code via:

import scanncaverages
from netCDF4 import Dataset
netcdffile='/ascldap/users/lcheung/GPFS1/2020/scitech/AIAAScitech2021.github/AMRWindRuns/stable/05ms/sim_dx2.5/post_processing_dx2.5/sampling60000.nc'
group          = 'p_h'
outputfile     = 'test.out'
ipts           = [0, 10]
jpts           = [0]
kpt            = 0
kpt_alpha      = 1
dx             = 2
dy             = 2
avgtimes       = [[15000.0, 15100.0], [15100, 15200]]
ncdata   = Dataset(netcdffile, 'r')
scanncaverages.scanTimeSeries(ncdata, ipts, jpts, kpt, dx, dy, avgtimes,outputfile, kpt_alpha=kpt_alpha, group=group, verbose=True)

"""

def findMatchingPt(ptlist, p, eps):
    for ipt, xpt in enumerate(ptlist):
        if (np.linalg.norm(np.array(xpt)-np.array(p)))<eps: return ipt
    # Error out
    raise Exception("error in findMatchingPt") 
        
def writeTimeSeries(ncdat, xvec, yvec, zvec, savestring, 
                    group='p_h', useindices=False, verbose=True):
    """
    Average the spectra over multiple x, y, z locations
    """
    #Nt     = ncdat.dimensions['num_time_steps'].size
    Npts   = ncdat[group].dimensions['num_points'].size
    allpts = ncdat[group].variables['coordinates']

    t      = ncdat['time'][:]
    allvx  = ncdat[group].variables['velocityx']
    allvy  = ncdat[group].variables['velocityy']
    allvz  = ncdat[group].variables['velocityz']
    Navg   = 0
    zerocol = np.zeros(len(t))
    all_ulongavgs = []
    for x in xvec:
        for y in yvec:
            for z in zvec:
                if useindices:
                    Nijk   = ncdat[group].ijk_dims
                    Ni     = Nijk[0]
                    Nj     = Nijk[1]
                    ipt    = int(x) + int(y)*Ni + int(z)*Ni*Nj
                    xyzpt = np.array([allpts[ipt, 0], 
                                      allpts[ipt, 1], 
                                      allpts[ipt, 2]])
                    icol   = x*np.ones(len(t))
                    jcol   = y*np.ones(len(t))
                    kcol   = z*np.ones(len(t))
                else:
                    xyzpt = np.array([x,y,z])
                    ipt   = findMatchingPt(allpts, xyzpt, 1.0E-6)
                    icol  = zerocol
                    jcol  = zerocol
                    kcol  = zerocol

                #print(allpts[ipt, :])
                u = allvx[:,ipt]
                v = allvy[:,ipt]
                w = allvz[:,ipt]
                xcol = xyzpt[0]*np.ones(len(t))
                ycol = xyzpt[1]*np.ones(len(t))
                zcol = xyzpt[2]*np.ones(len(t))
                savedat = np.vstack((t, 
                                     icol, jcol, kcol, 
                                     xcol, ycol, zcol, 
                                     u, v, w))
                fname=savestring%(x,y,z)
                np.savetxt(fname, savedat.transpose())
                if verbose: print("saved "+fname)

def scanTimeSeries(ncdat, ivec, jvec, kpt, dx, dy, avgtimes,
                   savefile, kpt_alpha=None, group='p_h', verbose=True):
    """
    Average the spectra over multiple x, y, z locations
    """
    Nt     = ncdat.dimensions['num_time_steps'].size
    Npts   = ncdat[group].dimensions['num_points'].size
    allpts = ncdat[group].variables['coordinates']

    # Construct time windows
    #twindows=[]
    #[twindows.append([times[i], times[i+1]]) for i in range(len(times)-1)]
    t      = ncdat['time'][:]
    allvx  = ncdat[group].variables['velocityx']
    allvy  = ncdat[group].variables['velocityy']
    allvz  = ncdat[group].variables['velocityz']

    if len(avgtimes)==0: twindows = [[t[0], t[-1]]]
    else:                twindows = avgtimes
    #print(twindows)

    Nijk   = ncdat[group].ijk_dims
    Ni     = Nijk[0]
    Nj     = Nijk[1]
    #print('Nijk = ',Nijk)
    #print('Nt   = ',Nt, len(t))
    header="# T1 T2    I J K   X Y Z   AVGU AVGV AVGW  STDU STDV STDW"
    if kpt_alpha is not None: header += " ALPHA" 
    if verbose: print(header)
    if len(savefile)>0:
        with open(savefile, 'w') as f:
            f.write(header+'\n')
            f.close()
    # Loop over time segments
    for tavg in twindows:
        #print(tavg)
        # Loop over all points
        for ix in ivec:
            for iy in jvec:
                # Get the base point
                basept = (ix) + (iy)*Ni + kpt*Ni*Nj
                xyzpt = np.array([allpts[basept, 0], 
                                  allpts[basept, 1], 
                                  allpts[basept, 2]])
                # Get all other points around it
                ipts=[]
                for deltai in range(dx): 
                    for deltaj in range(dy):
                        ipt    = (ix+deltai) + (iy+deltaj)*Ni + kpt*Ni*Nj
                        ipts.append(ipt)
                # Get average velocity over points in window
                tfilter = (tavg[0]<=t)&(t <= tavg[1])
                meanu_vec, meanv_vec, meanw_vec =[], [], []
                for ipt in ipts:
                    meanu_vec.append(np.mean(allvx[tfilter, ipt]))
                    meanv_vec.append(np.mean(allvy[tfilter, ipt]))
                    meanw_vec.append(np.mean(allvz[tfilter, ipt]))
                avgu = np.mean(meanu_vec)
                avgv = np.mean(meanv_vec)
                avgw = np.mean(meanw_vec)

                # Get stddev of velocity over points in window
                stdu_vec, stdv_vec, stdw_vec =[], [], []
                for ipt in ipts:
                    stdu_vec.append(np.std(allvx[tfilter, ipt])**2)
                    stdv_vec.append(np.std(allvy[tfilter, ipt])**2)
                    stdw_vec.append(np.std(allvz[tfilter, ipt])**2)
                stdu = np.sqrt(np.mean(stdu_vec))
                stdv = np.sqrt(np.mean(stdv_vec))
                stdw = np.sqrt(np.mean(stdw_vec))

                # Get the alpha exponent (if needed)
                if kpt_alpha is not None:
                    # Get the basic point
                    basept_alpha = (ix) + (iy)*Ni + kpt_alpha*Ni*Nj
                    xyzpt_alpha = np.array([allpts[basept_alpha, 0], 
                                            allpts[basept_alpha, 1], 
                                            allpts[basept_alpha, 2]])
                    
                    # Get all points around it
                    ipts_alpha=[]
                    for deltai in range(dx): 
                        for deltaj in range(dy):
                            ipt    = (ix+deltai) + (iy+deltaj)*Ni + kpt_alpha*Ni*Nj
                            ipts_alpha.append(ipt)
                    alphau_vec, alphav_vec, alphaw_vec =[], [], []
                    for ipt in ipts_alpha:
                        alphau_vec.append(np.mean(allvx[tfilter, ipt]))
                        alphav_vec.append(np.mean(allvy[tfilter, ipt]))
                        alphaw_vec.append(np.mean(allvz[tfilter, ipt]))
                    avgu_alpha = np.mean(alphau_vec)
                    avgv_alpha = np.mean(alphav_vec)
                    avgw_alpha = np.mean(alphaw_vec)
                    Uh_0       = np.sqrt(avgu**2 + avgv**2)
                    Uh_alpha   = np.sqrt(avgu_alpha**2 + avgv_alpha**2)
                    #print("z1 = %f z2 = %f"%(xyzpt[2], xyzpt_alpha[2]))
                    # ------------
                    # Following the formula (V2/V1) = (z2/z1)^alpha
                    #                       alpha = log(V2/V1)/log(z2/z1)
                    if kpt_alpha > kpt:
                        alpha_shear = np.log(Uh_alpha/Uh_0)/np.log(xyzpt_alpha[2]/xyzpt[2])
                    else:
                        alpha_shear = np.log(Uh_0/Uh_alpha)/np.log(xyzpt[2]/xyzpt_alpha[2])

                outstring=" ".join([repr(x) for x in 
                                    [tavg[0], tavg[1], ix, iy, kpt,
                                     xyzpt[0], xyzpt[1], xyzpt[2], 
                                     avgu, avgv, avgw, stdu, stdv, stdw]])
                if kpt_alpha is not None: outstring += " "+repr(alpha_shear)
                if verbose:
                    print(outstring)
                if len(savefile)>0:
                    with open(savefile, 'a') as f:
                        f.write(outstring+'\n')
                        f.close()
    return        

# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":

    helpstring="""
"""
    # Parse arguments
    parser = argparse.ArgumentParser(description=helpstring, 
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "ncfile",
        help="netcdf file",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--outputfile",
        help="output file format",
        dest='outputfile',
        default='',
        type=str,
    )
    parser.add_argument(
        "-g",
        "--group",
        help="NetCDF group, default=p_h",
        default='p_h',
        type=str,
    )
    parser.add_argument(
        "-i",
        nargs='+',
        help="i points to extract",
        dest='ipts',
        required=True,
    )

    parser.add_argument(
        "-j",
        nargs='+',
        help="j points to extract",
        dest='jpts',
        required=True,
    )

    parser.add_argument(
        "-k",
        #nargs='+',
        help="k point to extract",
        dest='kpt',
        required=True,
    )
    parser.add_argument(
        "--dx",
        help="Number of points in x",
        dest='dx',
        default=1,
        #type=str,
    )
    parser.add_argument(
        "--dy",
        help="Number of points in y",
        dest='dy',
        default=1,
        #type=str,
    )
    parser.add_argument(
        "--avgtimes",
        help="A comma separated string with pairs of averaging time windows",
        dest='avgtimes',
        default='',
        type=str,
    )

    # Load the options
    args      = parser.parse_args()
    filename  = args.ncfile
    outputfile= args.outputfile
    group     = args.group
    ipts      = [int(x) for x in args.ipts]
    jpts      = [int(y) for y in args.jpts]
    kpt       = int(args.kpt)
    dx        = int(args.dx)
    dy        = int(args.dy)
    if len(args.avgtimes)==0: 
        avgtimes = []
    else: 
        avgtimes = [[float(y) for y in x.split()] for x in args.avgtimes.split(',')]

    print("netcdf file    = "+filename)
    print("group          = \'"+group+"\'")
    print("outputfile     = \'"+outputfile+"\'")
    #print("using indices  = "+repr(useindices))
    print("ipts           = "+repr(ipts))
    print("jpts           = "+repr(jpts))    
    print("kpt            = "+repr(kpt))
    print("dx             = "+repr(dx))
    print("dy             = "+repr(dy))
    print("avgtimes       = "+repr(avgtimes))

    ncdata   = Dataset(filename, 'r')
    scanTimeSeries(ncdata, ipts, jpts, kpt, dx, dy, avgtimes,
                   outputfile, group=group, verbose=True)
