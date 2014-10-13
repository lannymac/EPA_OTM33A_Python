import numpy as np
from scipy.optimize import curve_fit
import defs as d
import pylab as pl; pl.close('all')
'''
This software uses EPA's OTA33A Method. From concentration and wind data
the emission rate from a point source can be estimated.

Any functions that begin with the prefix "d.", information can be found
within the defs.py file.

'''
################
#    INPUTS    #
################

make_plot      = True

# DATA FILES AND SOURCE-RECEPTOR DISTANCE
files = ['Christman_20140919/CHRISTMAN_2014_09_19_measurement_%d.npz'% (k) for k in range(1,10)]


distance     = [45,75,75,75,175,175,175,150,75]
actual_flow_rate = [5,7.5,7.5,10,5,7.5,7.5,7.5,2.5]


# TRACER INFORMATION
chemical_name = 'CH4'
mw_chemical   = 16.04 # [g mol-1]
#mw_chemical   = 26.04
# BIN LIMITS
theta_start   = 0.
theta_end     = 360.
delta_theta   = 10.

# DATA CUTOFFS
wslimit       = 0.0 # set wind speed cut limit [m s-1]
wdlimit       = 60.0 # set wind angle cut limit [+/- deg]. 180 indicates no filter
cutoff        = 2 # bin density cutoff limit

############################
#    BEGIN CALCULATION     #
############################
pd = np.zeros((len(files)))
cc = np.zeros((len(files)))
for k1 in range(len(files)):
    print('+'*(len(files[k1])+10))
    print('+    %s    +' % (files[k1]))
    print('+'*(len(files[k1])+10)+'\n')

    # load in the data
    #tracer,ws3,wd3,ws2,wd2,temp,pres,ws3z,ws3x,ws3y,time = d.load_excel(files[k1])
    tracer,lat,lon,ws3,wd3,ws2,wd2,temp,pres,ws3z,ws3x,ws3y,ws3t,time = d.custom_load_files(files[k1],chemical_name)
    # correct wind from sonic anemometer
    ws3x, ws3y, ws3z, wd3, ws3 = d.sonic_correction(ws3x,ws3y,ws3z)

    # engage wind speed limit cutoff
    if wslimit > 0.0:
        # create a mask array where False indices are where wind speed < wslimit
        ws_mask = np.ma.masked_less_equal(ws3,wslimit).mask
        
        # apply the created mask to the relevant variables
        time.mask   = ws_mask
        tracer.mask = ws_mask
        ws3.mask    = ws_mask
        wd3.mask    = ws_mask
        ws2.mask    = ws_mask
        wd2.mask    = ws_mask
        temp.mask   = ws_mask
        pres.mask   = ws_mask
        ws3z.mask   = ws_mask

    else:
        # if we don't want to mask any wind speeds, we have to make sure that we 
        # don't lose any data during the np.percentile step below
        ws_mask = np.ones_like(tracer.mask)*False
        
    # subtract tracer background. I had to manually apply the mask within the percentile function
    # because currently np.percentile does not play well with masked arrays.
    tracer -= np.mean(tracer[np.where(tracer < np.percentile(tracer[~ws_mask],5))]) # [g m-3]
    # there are bound to be some negative values in the tracer array now
    # let's set them to zero now

    tracer[np.where(tracer <0.)] = 0.

    # create wind direction bins from input
    bins = np.arange(theta_start,theta_end+1,delta_theta)

    # get indices for which bin every measured wind speed is in
    indices= np.digitize(wd3[~ws_mask],bins)-1
    
    # make array of the middle of each bin; len(mid_bins) = len(bins) - 1
    mid_bins = (bins[:-1] + bins[1:])/2.

    # create empty arrays to store mean concentration as a function of wind direction
    tracer_avg = np.zeros_like(mid_bins)

    # for each wind direction bin, find the corresponding mean tracer concentration
    for i in range(indices.min(),indices.max()+1):
        temp_vals = tracer[~ws_mask][np.where(indices == i)]

        # ensure that the amount of data points exceeds the cutoff value
        # we don't want bins with a couple of values to dominate the "gaussian" curve
        if len(temp_vals) > int(cutoff*len(tracer)/100.):
            tracer_avg[i] = temp_vals.mean()

    # ensure that the peak concentration is around 180 degrees for fitting
    roll_amount = int(len(tracer_avg)/2. -1) - np.where(tracer_avg == tracer_avg.max())[0][0]
    tracer_avg = np.roll(tracer_avg,roll_amount)
    
    # get the bin with peak average concentration
    max_bin = mid_bins[np.where(tracer_avg == tracer_avg.max())[0][0]-1]
    
    # calculate wind direction cut off based on input value
    bin_cut_lo = max_bin - wdlimit
    bin_cut_hi = max_bin + wdlimit

    # add the roll amount to wind direction to correspond correctly to the
    # average concentration array
    wd3 = wd3 + (roll_amount-1)*delta_theta

    # ensure that wind direction is within 0 to 360
    wd3 = d.wrap(wd3)

    # get the indices of where wd3 is outside of wind direction limit
    wd_mask = np.where((wd3 > bin_cut_hi) | (wd3 < bin_cut_lo))

    # mask arrays where wd3 is outside of wind direction cutoff
    wd3[wd_mask] = np.ma.masked
    wd2[wd_mask] = np.ma.masked
    ws3[wd_mask] = np.ma.masked
    ws3z[wd_mask] = np.ma.masked
    temp[wd_mask] = np.ma.masked
    pres[wd_mask] = np.ma.masked

    # fitting procedure
    # here are some initial guesses that produce good results
    const_0 = [tracer_avg.max(),24,180]
    
    # the curve fit function uses the Levenberg-Marquardt algorithm which
    # in my opionion does a great job
    fit_tracer,cov_tracer = curve_fit(d.gaussian_func,mid_bins,tracer_avg,p0 = const_0) # fit coefficients
    if make_plot: d.fit_plot(mid_bins,tracer_avg,fit_tracer,chemical_name) # make plot if you want

    # calculate the standard deviation of wind direction and turbulent intensity 
    # for use in finding the PG stability class
    turbulent_intensity = np.std(ws3z)/np.mean(ws3) # turbulent intensity
    std_wind_dir = d.yamartino_method(wd2) # st. dev. of wind direction [deg]

    # calcualte the vertical and horizontal dispersion of the gaussian plume
    # using PGT stability classes
    sy, sz = d.sigma(distance[k1],std_wind = std_wind_dir,turb = turbulent_intensity,tables=True,stab=None) # [m]

    # convert tracers from ppm to g m-3
    fit_amplitude = d.ppm2gm3(fit_tracer[0],mw_chemical,np.mean(temp),np.mean(pres))

    # calculate the density of the tracers using the ideal gas law for conversion back to L/min if desired
    rho_tracer_stp = (d.P_stp*mw_chemical)/(d.R*d.T_stp) # density of gas at STP [g L-1]
    rho_tracer = (np.mean(pres)*mw_chemical)/(d.R*np.mean(temp)) # density of gas at ambient conditions [g L-1]
   
    # putting it all together to calculate the emission rate in g/s
    emission_tracer_mass_per_time = 2*np.pi*fit_amplitude*np.mean(ws3)*sy*sz # [g s-1]

    # now calculate the emission rate in L/min
    emission_tracer_volume_per_time = 2*np.pi*fit_amplitude*np.mean(ws3)*sy*sz/rho_tracer*60 # [L min-1]

    # now L/min at STP
    emission_tracer_volume_per_time_stp = 2*np.pi*fit_amplitude*np.mean(ws3)*sy*sz/rho_tracer_stp*60 # [L min-1]

    # let us finally print our results to the screen
    print('%s is predicted to have an emission rate of:\n' % (chemical_name))
    print('\t%.3f g/s' % (emission_tracer_mass_per_time))
    print('\t%.3f LPM' % (emission_tracer_volume_per_time))
    print('\t%.3f SLPM\n' % (emission_tracer_volume_per_time_stp))
    pd[k1] = (emission_tracer_volume_per_time_stp - actual_flow_rate[k1])/actual_flow_rate[k1]*100.
    print('\t%.3f percent difference' % (pd[k1]))
    

    cc[k1] = np.corrcoef(tracer_avg,d.gaussian_func(mid_bins,fit_tracer[0],fit_tracer[1],fit_tracer[2]))[0][1]**2
    print('\tR^2 = %.3f' % (cc[k1]))


fig=pl.figure()
ax=fig.add_subplot(111)
ax.scatter(cc,pd,s=40)

rangex = np.diff((ax.get_xlim()))
rangey = np.diff((ax.get_ylim()))

for i in range(len(files)):
    ax.text(cc[i]+.01*rangex,pd[i]+.01*rangey,'%d' % (i+1))
pl.axis('tight')
ax.plot(ax.get_xlim(),[0,0],'--',c='k',lw=2)
ax.set_xlabel('R$^{2}$',fontsize=18)
ax.set_ylabel('Percent Difference',fontsize=18)
pl.axis('tight')
pl.show()
