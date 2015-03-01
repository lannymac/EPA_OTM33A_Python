import numpy as np
from scipy.optimize import curve_fit
import defs as d
import pylab as pl; pl.close('all')
import sys
'''
This software uses EPA's OTA33A Method. From concentration and wind data
the emission rate from a point source can be estimated.

Any functions that begin with the prefix "d.", information can be found
within the defs.py file.

This program was coded without the ability to automatically read in all
different types of files. This is left up to the user. This program will
only provide a function to perform the data analysis.
'''
################
#    INPUTS    #
################

# WOULD YOU LIKE TO SEE A PLOT?
make_plot      = True

# DATA FILES AND SOURCE-RECEPTOR DISTANCE
ch4,ws3,wd3,ws2,wd2,temp,pres,ws3z,ws3y,ws3x,time = d.load_excel('STR_3061611_01.xls')

distance = 42. # meters

# TRACER INFORMATION
chemical_name = 'CH4' # for plots
mw_chemical   = 16.04 # [g mol-1]

# BIN LIMITS
theta_start   = 5.
theta_end     = 365.
delta_theta   = 10.

# DATA CUTOFFS
wslimit       = 0.0 # set wind speed cut limit [m s-1]
wdlimit       = 60.0 # set wind angle cut limit [+/- deg]. 180 indicates no filter
cutoff        = 2 # bin density cutoff limit

############################
#    BEGIN CALCULATION     #
############################

volumeRate, massRate = d.OTA33A(ch4,temp,pres,ws3z,ws3x,ws3y,time,wslimit,wdlimit,cutoff,distance,theta_start,theta_end,delta_theta,mw_chemical,chemical_name,make_plot=make_plot)

# let us  print our results to the screen
print('%s is predicted to have an emission rate of:\n' % (chemical_name))
print('\t%.3f g/s' % (massRate))
print('\t%.3f SLPM' % (volumeRate))

