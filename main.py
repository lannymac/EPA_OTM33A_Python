import numpy as np
from scipy.optimize import curve_fit
from defs import *
import sys
import pylab as pl; pl.close('all')

'''
This software is based on EPA's OTA33A Method. From concentration and wind data
the emission rate from a point source can be estimated. This program will assume
a point source gaussian plume to estimate the emission rate.
'''

make_plot = True
print_2_screen = True

files = ['CHRISTMAN_2014_03_14_1010_TO_1034.npz']
d = [42.]

chemical_1 = 'CH4'
chemical_2 = 'C2H2'

mw_chemical_1 = 16.04 # g mol-1
mw_chemical_2 = 26.04 # g mol-1

theta_start = 0.
theta_end = 360.
delta_theta = 5.

cutoff = 0.


for k1 in range(len(files)):


    tracer_1,tracer_2,lat,lon,ws3,wd3,ws2,wd2,temp,\
        pres,ws3z,ws3x,ws3y,ws3t,time = custom_load_files(files[k1],chemical_1,chemical_2)

    # convert tracers from ppm to g m-3
    tracer_1 = ppm2gm3(tracer_1,mw_chemical_1,np.mean(temp),np.mean(pres))
    tracer_2 = ppm2gm3(tracer_2,mw_chemical_2,np.mean(temp),np.mean(pres))

    # calculate mean wind speed and turbulent intensity 
    # for use in finding the PG stability class
    mean_wind_speed = np.mean(ws3)
    turbulent_intensity = np.std(ws3z)/mean_wind_speed

    #subtract tracer background
    tracer_1 = tracer_1 - np.percentile(tracer_1,05)
    tracer_1[np.where(tracer_1 <0.)] = 0.

    tracer_2 = tracer_2 - np.percentile(tracer_2,05)
    tracer_2[np.where(tracer_2 <0.)] = 0.

    # create bins
    bins = np.arange(theta_start,theta_end+1,delta_theta)

    # get array for which bin every measured wind speed is in
    indices= np.digitize(wd2,bins)-1
    
    # make array of the middle of each bin; len(mid_bins) = len(bins) - 1
    mid_bins = (bins[:-1] + bins[1:])/2.

    # create empty arrays to store "histograms"
    hist_tracer_1 = np.zeros_like(mid_bins)
    hist_tracer_2 = np.zeros_like(mid_bins)    

    # for each wind direction bin, find the corresponding mean tracer concentration
    for i in range(indices.min(),indices.max()+1):
        temp_vals = tracer_1[np.where(indices == i)]
        if len(temp_vals) > cutoff:
            hist_tracer_1[i] = temp_vals.mean()

        temp_vals = tracer_2[np.where(indices == i)]
        if len(temp_vals) > cutoff:
            hist_tracer_2[i] = temp_vals.mean()
        
    # ensure that the peak concentration is around 180 degrees for fitting
    roll_amount_1 = int(len(hist_tracer_1)/2. -1) - np.where(hist_tracer_1 == hist_tracer_1.max())[0][0]
    roll_amount_2 = int(len(hist_tracer_2)/2. -1) - np.where(hist_tracer_2 == hist_tracer_2.max())[0][0]

    hist_tracer_1 = np.roll(hist_tracer_1,roll_amount_1)
    hist_tracer_2 = np.roll(hist_tracer_2,roll_amount_2)

    # fitting procedure
    const_0 = [tracer_1.max(),24,180]
    fit_tracer_1,cov_tracer_1 = curve_fit(func,mid_bins,hist_tracer_1,p0 = const_0)
    if make_plot: fit_plot(mid_bins,hist_tracer_1,fit_tracer_1,chemical_1)

    const_0 = [tracer_2.max(),24,180]
    fit_tracer_2,cov_tracer_2 = curve_fit(func,mid_bins,hist_tracer_2,p0 = const_0)
    if make_plot: fit_plot(mid_bins,hist_tracer_2,fit_tracer_2,chemical_2)

    # calcualte the vertical and horizontal dispersion of the gaussian plume
    # using PGT stability classes
    sy, sz = sigma(d[k1],std_wind = np.std(wd3),turb = turbulent_intensity)

    # calculate the density of the tracers using the ideal gas law
    R = 0.082 # L atm K-1 mol-1
    rho_tracer_1_stp = (1.*mw_chemical_1)/(R*298.)
    rho_tracer_2_stp = (1.*mw_chemical_2)/(R*298.)

    rho_tracer_1 = (np.mean(pres)*mw_chemical_1)/(R*np.mean(temp))
    rho_tracer_2 = (np.mean(pres)*mw_chemical_2)/(R*np.mean(temp))
    
    # putting it all together to calculate the emission rate in L/min
    emission_tracer_1 = 2*np.pi*fit_tracer_1[0]*mean_wind_speed*sy*sz/rho_tracer_1*60
    emission_tracer_2 = 2*np.pi*fit_tracer_2[0]*mean_wind_speed*sy*sz/rho_tracer_2*60

    # now L/min at STP
    emission_tracer_1_stp = 2*np.pi*fit_tracer_1[0]*mean_wind_speed*sy*sz/rho_tracer_1_stp*60
    emission_tracer_2_stp = 2*np.pi*fit_tracer_2[0]*mean_wind_speed*sy*sz/rho_tracer_2_stp*60

    if print_2_screen: print_filename(files[k1]); \
       print_screen(chemical_1,emission_tracer_1,emission_tracer_1_stp); \
       print_screen(chemical_2,emission_tracer_2,emission_tracer_2_stp)
