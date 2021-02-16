#!/usr/bin/env python
#

import numpy as np
import sys, os
import argparse
from netCDF4 import Dataset
#from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Process
from multiprocessing import Semaphore

utildir='../../utilities'
sys.path.insert(1, utildir)
import scanncaverages

# Set the defaults here
# ==================================
default_netcdffile ='./sampling30000.nc' 
default_group      = 'p_hub'
default_outputfile = 'output_j'
default_ipts       = np.arange(0,255,4) #[0, 10] 
default_jpts       = [0, 10] 
default_kpt        = 3
default_kpt_alpha  = 4
default_dx         = 1
default_dy         = 1
default_tstart     = 15000
default_tend       = 20000
default_deltat     = 240   # 4 min
# ==================================


# See https://stackoverflow.com/questions/38987/how-do-i-merge-two-dictionaries-in-a-single-expression-in-python-taking-union-o
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def scanpoint(arg, sema, verbose):
    # For semaphore usage: 
    # https://stackoverflow.com/questions/20886565/using-multiprocessing-process-with-a-maximum-number-of-simultaneous-processes
    print("Working on j="+repr(arg['jpts']))
    ncdata         = Dataset(arg['netcdffile'], 'r')
    scanncaverages.scanTimeSeries(ncdata, 
                                  arg['ipts'], arg['jpts'], arg['kpt'], 
                                  arg['dx'], arg['dy'], arg['avgtimes'],
                                  arg['outputfile'], 
                                  kpt_alpha=arg['kpt_alpha'], 
                                  group=arg['group'], verbose=verbose)
    ncdata.close()
    sema.release()
    return

def main():
    # Handle arguments
    parser = argparse.ArgumentParser(description='Find the right conditions')
    parser.add_argument('--ipts',        
                        help="i indices to probe",
                        type=int, nargs='+')
    parser.add_argument('--jpts',        
                        help="j indices to probe",
                        type=int, nargs='+')
    parser.add_argument('--outfileprefix',        
                        help="Output filename prefix",
                        default=default_outputfile, nargs='+')
    parser.add_argument('--np',        
                        help="Number of process threads to use",
                        type=int, default=2)    
    parser.add_argument('-v','--verbose',        
                        help="Turn on verbose mode",
                        action='store_true', default=False)    
    parser.add_argument(
        "-o",
        "--outputfile",
        help="output file format",
        dest='outputfile',
        default='',
        type=str,
    )
    # Get the default and user arguments
    args=parser.parse_args()

    netcdffile     = default_netcdffile   # '/ascldap/users/lcheung/GPFS1/2020/scitech/AIAAScitech2021.github/AMRWindRuns/stable/05ms/sim_dx2.5/post_processing_dx2.5/sampling60000.nc'
    group          = default_group # 'p_h'
    outputfile     = args.outfileprefix
    ipts           = default_ipts if args.ipts == None else args.ipts
    jpts           = default_jpts if args.jpts == None else args.jpts
    kpt            = default_kpt
    kpt_alpha      = default_kpt_alpha
    dx             = default_dx
    dy             = default_dy

    # Options for how to run
    verbose        = args.verbose
    nprocs         = args.np

    # Build the average times list
    tstart         = default_tstart
    tend           = default_tend
    deltat         = default_deltat
    avgtimes       = []
    for t in np.arange(tstart, tend, deltat):
        if t+deltat > tend: continue
        avgtimes.append([t, t+deltat])
    #avgtimes       = [[15000.0, 15100.0], [15100, 15200]]

    print("outputfile="+outputfile)
    print("ipts="+repr(ipts))
    print("jpts="+repr(jpts))
    print("nprocs="+repr(nprocs))
    print("avgtimes="+repr(avgtimes))

    # Set up the dictionary
    baseargs = {'netcdffile':netcdffile,
                'group':group, 
                'outputfile': outputfile,
                'ipts':ipts,
                'kpt':kpt,
                'kpt_alpha':kpt_alpha,
                'dx':dx,
                'dy':dy,
                'avgtimes':avgtimes,
            }

    # Divide up work into a list
    worklist=[]
    for j in jpts:
        joutfile=outputfile+repr(j)+'.dat'
        jdict = merge_two_dicts(baseargs, {'outputfile':joutfile, 'jpts':[j]})
        worklist.append(jdict)

    # Print the work list
    sema = Semaphore(nprocs)
    process=[]
    for d in worklist:
        sema.acquire()
        p = Process(target=scanpoint, args=(d, sema, verbose))
        process.append(p)
        p.start()

    # join all the processes when end
    for p in process:
        p.join()

    return

if __name__ == "__main__":
    # Call the main function
    main()
