# Important header information
naluhelperdir = '../../utilities'
# Import libraries
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(1, naluhelperdir)
import postproamrwind
# %matplotlib inline


# Load the data
filename='post_processing/abl_statistics00000.nc'
avgtimes=[15000,20000]
ablstats=postproamrwind.loadData(filename,avgt=avgtimes)

#print(postproamrwind.stdvars)
# Save the data into text files
postproamrwind.saveAsText(ablstats, ['u','v','w'], 'amrwind_velocity.dat', header='z, u, v, w')
postproamrwind.saveAsText(ablstats, ['theta'],     'amrwind_temperature.dat', header='z, T')
postproamrwind.saveAsText(ablstats, [u"u'u'_r", u"u'v'_r", u"u'w'_r", u"v'v'_r", u"v'w'_r", u"w'w'_r"],     
                        'amrwind_reynoldsstresses.dat', header='z, uu, uv, uw, vv, vw, ww')
