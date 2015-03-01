import numpy as np
import pylab as pl
import xlrd
from scipy.optimize import curve_fit

###################
#    CONSTANTS    #
###################

R = 0.08205746 # universal gas constant [L atm K-1 mol-1]
P_stp = 1. # standard atmospheric pressure [atm]
T_stp = 298. # standard temperature [K]

###################
#    FUNCTIONS    #
###################

def gaussian_func(x,a,sigma,x_0):
    '''
    This function will return a gaussian with the
    proper coefficients
    '''
    return a*np.exp(-((x - x_0)**2)/(2*sigma**2))

def stability_class(std_wind,turb,v=False):
    '''
    This function will return an average stability class
    based on the standard deviation of the wind speed 
    and the turbulent intensity.

    INPUTS
    std_wind :: standard deviation on the horizontal wind direction [degrees]
    turb     :: turbulent intensity
    '''

    if std_wind > 27.5:
        pg_wind = 1
    elif std_wind <=27.5 and std_wind > 23.5:
        pg_wind = 2
    elif std_wind <=23.5 and std_wind >19.5:
        pg_wind = 3
    elif std_wind <= 19.5 and std_wind > 15.5:
        pg_wind = 4
    elif std_wind <= 15.5 and std_wind > 11.5:
        pg_wind = 5
    elif std_wind <= 11.5 and std_wind > 7.5:
        pg_wind = 6
    elif std_wind <= 7.50:
        pg_wind = 7
    if v:
        print("Horizontal Stability class = %d "%(pg_wind))
    if turb > 0.205:
        pg_turb = 1
    elif turb <= .205 and turb > 0.180:
        pg_turb = 2
    elif turb <= 0.180 and turb > 0.155:
        pg_turb = 3
    elif turb <=0.155 and turb > 0.130:
        pg_turb = 4
    elif turb <= 0.130 and turb > 0.105:
        pg_turb = 5
    elif turb <= 0.105 and turb > 0.080:
        pg_turb = 6
    elif turb <= 0.080:
        pg_turb = 7
    if v:
        print("Vertical Stability class = %d "%(pg_wind))
    final_pg = int(round(np.mean((pg_wind,pg_turb))))
    return final_pg


def sigma_func(I,J,K,dist):
    '''
    This function is used to estimate sigma_y and sigma_z based on 
    stability parameters.
    '''
    return np.exp(I + J*np.log(dist) + K*(np.log(dist)**2))
    
def sigma(dist,std_wind=None,turb=None,stab= None,tables=False):
    
    '''
    This definition will calculate the horizontal and vertical standard deviation of the emission
    distribution. The seven custom classes are based on the Pasquill-Gifford model in the following way:

    1 - PG class A
    2 - Linear interpolation between PG classes A and B
    3 - PG class B
    4 - Linear interpolation between PG classes B and C
    5 - PG class C
    6 - Linear interpolation between PG classes C and D
    7 - PG class D

    If you manually enter a stability class (i.e. 1-7) it will over-ride the conversion of 
    std_wind and turb into a stability class.

    INPUTS
    dist     :: distance between receptor and source
    std_wind :: standard deviation of the horizontal wind direction [degrees]
    turb     :: turbulent intensity
    stab     :: stability class

    '''
    if stab == None:
        stab = stability_class(std_wind,turb)
    file_y = np.load('pgtabley.npz')
    pgtabley = file_y['data']
    sigma_y = pgtabley[round(dist)-1,stab-1]

    file_y = np.load('pgtablez.npz')
    pgtablez = file_y['data']
    sigma_z = pgtablez[round(dist)-1,stab-1]

    return sigma_y,sigma_z

def ppm2gm3(conc,mw,T,P):
    '''
    This definition will convert a concentration given in ppm to
    a concentration in g/m3.
    
    INPUTS
    conc :: concentration [ppb]
    mw   :: molecular weight [g/mol]
    T    :: temperature [K]
    P    :: pressure [atm]
    '''
    R = 0.082 # L atm K-1 mol-1
    mass = conc*mw
    V = (1e6*R*T)/(P*1000.)
    return mass/V

def fit_plot(x,y_data,fit,name):
    '''
    This definition will plot the actual data versus the Gaussian fit performed

    INPUTS

    x      :: wind direction bins
    y_data :: calculated average concentration
    fit    :: the fit coefficients for a gaussian function
    name   :: this string will be made the title
    '''
    y_fit = gaussian_func(x,fit[0],fit[1],fit[2])
    
    fs=18
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y_data,c='k',label='Data',lw=2)
    ax.plot(x,y_fit,c='r',label='Fit',lw=2)    
    ax.set_xlabel('Wind Direction [degrees]',fontsize=fs)
    ax.set_ylabel('Average concentration [ppm]',fontsize=fs)
    ax.set_title(name,fontsize=fs+2)
    ax.set_xlim(0,360)
    leg=ax.legend(fontsize=fs)
    leg.get_frame().set_alpha(0.)
    pl.show()

def load_excel(filename):
    '''
    This function was created to read in the example data from an excel file.

    Ideally when using this program, the user will create their own small script
    to properly read in their data.
    '''
    n1=53 # Remove sampling time deplay of concetration and Sonic

    workbook = xlrd.open_workbook(filename) # load Excel file
    sheet = workbook.sheet_by_index(0) # use the first sheet in the file

    # create time array
    tmp = np.array([sheet.cell_value(row,0) for row in range(1+n1,sheet.nrows)])
    time = []
    for j in range(len(tmp)):
        time.append(xlrd.xldate_as_tuple(tmp[j],workbook.datemode))
    time = np.ma.array(time,mask=False)

    #grab all the other variables
    tracer = np.ma.array([sheet.cell_value(row,1) for row in range(n1+1,sheet.nrows)],mask=False)
    ch4 = np.ma.array([sheet.cell_value(row,2) for row in range(n1+1,sheet.nrows)],mask=False)

    lat = np.ma.array([sheet.cell_value(row,3) for row in range(1,sheet.nrows-n1)],mask=False)
    lon = np.ma.array([sheet.cell_value(row,4) for row in range(1,sheet.nrows-n1)],mask=False)
    ws3 = np.ma.array([sheet.cell_value(row,5) for row in range(1,sheet.nrows-n1)],mask=False)
    wd3 = np.ma.array([sheet.cell_value(row,6) for row in range(1,sheet.nrows-n1)],mask=False)
    ws2 = np.ma.array([sheet.cell_value(row,7) for row in range(1,sheet.nrows-n1)],mask=False)
    wd2 = np.ma.array([sheet.cell_value(row,8) for row in range(1,sheet.nrows-n1)],mask=False)
    temp =np.ma.array([sheet.cell_value(row,9) for row in range(1,sheet.nrows-n1)],mask=False)+273.15
    pres =np.ma.array([sheet.cell_value(row,11) for row in range(1,sheet.nrows-n1)],mask=False)/1000.
    ws3z =np.ma.array([sheet.cell_value(row,12) for row in range(1,sheet.nrows-n1)],mask=False)
    ws3y =np.ma.array([sheet.cell_value(row,14) for row in range(1,sheet.nrows-n1)],mask=False)
    ws3x =np.ma.array([sheet.cell_value(row,13) for row in range(1,sheet.nrows-n1)],mask=False)
    ws3t =np.ma.array([sheet.cell_value(row,15) for row in range(1,sheet.nrows-n1)],mask=False)

    return ch4,ws3,wd3,ws2,wd2,temp,pres,ws3z,ws3y,ws3x,time

def yamartino_method(wd):
    '''
    This function will calculate the standard deviation of the wind direction using the Yamartino
    Method. More information on the method can be found on Wikipedia (http://en.wikipedia.org/wiki/Yamartino_method)
    or from the original paper itself:

    Yamartino, R.J. (1984). "A Comparison of Several "Single-Pass" Estimators of the Standard Deviation of
    Wind Direction". Journal of Climate and Applied Meteorology 23 (9): 1362-1366.
    '''
    n = float(len(wd))
    sa = np.mean(np.sin(wd*np.pi/180))
    ca = np.mean(np.cos(wd*np.pi/180))
    
    epsilon = np.sqrt(1 - (sa**2 + ca**2))
    std_wd = np.sqrt(n/(n-1))*np.arcsin(epsilon)*(1 + (2./np.sqrt(3) - 1)*epsilon**3)

    return std_wd*180/np.pi


def wrap(wd):
    '''
    This function will ensure that any wind direction array fed into it will only
    have values between 0 and 360.

    INPUTS

    wd :: array of wind directions [degrees]
    '''
    wd[np.where(wd >360.)] -= 360.
    wd[np.where(wd <0.)] += 360.
    
    return wd

def sonic_correction(u,v,w):
    '''
    Rotate 3D sonic coordinate system to streamlined set (-180 deg)
    Field data should be near -180 deg by manual sonic rotation set
    If mean direction is in N hemispehre, rotation sets to 0 deg (warning)
    Output <x>, <z>, = 0 or check file
    '''
    mdir = 180 + (np.arctan2(-1*(np.mean(v)),-1*(np.mean(u)))*180/np.pi)

    RA = np.arctan(np.mean(v)/np.mean(u)) # First rotation (set <x> = 0)

    CRA = np.cos(RA)
    SRA = np.sin(RA)
    
    u_rot = u*CRA + v*SRA
    v_rot = -u*SRA + v*CRA
    
    mndir = 180 + (np.arctan2(-1*(np.mean(v_rot)),-1*(np.mean(u_rot)))*180/np.pi)
    
    
    RB = np.arctan(np.mean(w)/np.mean(u_rot)) # Second rotation (set <z> = 0)
    
    CRB = np.cos(RB)
    SRB = np.sin(RB)
    
    u_rot2 = u_rot*CRB + w*SRB
    w_rot = (-1)*u_rot*SRB + w*CRB
    
    mnnws3x = np.mean(u_rot2)
    mnws3z = np.mean(w_rot)
    mmndir = 180 + (np.arctan2(-1*(mnws3z),-1*(mnnws3x))*180/np.pi)

    # calculate and asign new 3D sonic 2D wind direction
    wd3 = 180 + (np.arctan2(-1*(v_rot),-1*(u_rot2))*180/np.pi)

    # calculate and asign new 3D sonic 2D wind speed
    ws3 = np.sqrt(u_rot2**2 + v_rot**2)

    return u_rot2,v_rot,w_rot, wd3, ws3


def custom_load_files(filename,chemical):
    '''
    This is a custom file loader that I made since
    I store my files in npz format. 
    '''
    a = np.load(filename)
    tracer = np.ma.array(a[chemical],mask=False)
    lat = np.ma.array(a['lat'],mask=False)
    lon = np.ma.array(a['lon'],mask=False)
    ws3 = np.ma.array(a['ws3'],mask=False)
    wd3 = np.ma.array(a['wd3'],mask=False)
    ws2 = np.ma.array(a['ws2'],mask=False)
    wd2 = np.ma.array(a['wd2'],mask=False)
    temp = np.ma.array(a['temp'],mask=False)+273.
    pres = np.ma.array(a['pres'],mask=False)/1000.
    ws3z = np.ma.array(a['ws3z'],mask=False)
    ws3x = np.ma.array(a['ws3x'],mask=False)
    ws3y = np.ma.array(a['ws3y'],mask=False)
    ws3t = np.ma.array(a['ws3t'],mask=False)
    time = np.ma.array(a['time'],mask=False)
    
    return tracer,lat,lon,ws3,wd3,ws2,wd2,temp,pres,ws3z,ws3x,ws3y,ws3t,time


def OTA33A(gasConc,temp,pres,ws3z,ws3x,ws3y,time,wslimit,wdlimit,cutoff,distance,theta_start,theta_end,delta_theta,mw_chemical,chemical_name,make_plot=False):
    '''
    This definition is the meat of the program. It will accept the following arguments:
    ################
    #    INPUTS    #
    ################

    gasConc     :: concentration [ppb] of gas as a function of time
    temp        :: temperature [K] for conversion to mass units
    pres        :: pressure [Pa] for conversion
    ws3z        :: vertical wind component [m/s] from sonic anemometer
    ws3x        :: x wind component [m/s] from sonic anemometer
    ws3y        :: y wind component [m/s] from sonic anemometer
    time        :: time array (not actually used, nice to have for debugging though)
    wslimit     :: lower limit on wind speeds [m/s]
    wdlimit     :: range of wind directions to accept (180 means no filter)
    cutoff      :: cutoff percentage for bins (i.e. if bin has less than cutoff % of all values, bin = 0)
    distance    :: distance from source to receptor [m]
    theta_start :: starting wind direction for bins
    theta_end   :: ending wind direction for bins
    delta_theta :: wind direction bin size
    mw_chemical :: molecular weight of gas [g/mol] for conversion to mass units
    make_plot   :: True if you want to see a plot of the distribution and the gaussian fit
    

    and the program will give you:

    #################
    #    OUTPUTS    #
    #################

    emission_gasConc_volume_per_time :: emission rate in SLPM
    emission_gasConc_mass_per_time   :: emission rate in [g/s]
    '''
    # correct 3D sonic data so that mean(wsz) ~ 0.
    ws3y, ws3x, ws3z, wd3, ws3 = sonic_correction(ws3x,ws3y,ws3z)

    # engage wind speed limit cutoff
    if wslimit > 0.0:
        # create a mask array where False indices are where wind speed < wslimit
        ws_mask = np.ma.masked_less_equal(ws3,wslimit).mask
        
        # apply the created mask to the relevant variables
        time.mask   = ws_mask
        gasConc.mask = ws_mask
        ws3.mask    = ws_mask
        wd3.mask    = ws_mask
        temp.mask   = ws_mask
        pres.mask   = ws_mask
        ws3z.mask   = ws_mask

    else:
        # if we don't want to mask any wind speeds, we have to make sure that we 
        # don't lose any data during the np.percentile step below
        ws_mask = np.ones_like(gasConc.mask)*False
        
    # subtract tracer background. I had to manually apply the mask within the percentile function
    # because currently np.percentile does not play well with masked arrays.
    gasConc -= np.mean(gasConc[np.where(gasConc < np.percentile(gasConc[~ws_mask],5))]) # [g m-3]
    # there are bound to be some negative values in the tracer array now
    # let's set them to zero now

    gasConc[np.where(gasConc <0.)] = 0.

    # create wind direction bins from input
    bins = np.arange(theta_start,theta_end+delta_theta,delta_theta)

    # get indices for which bin every measured wind speed is in
    indices= np.digitize(wd3[~ws_mask],bins)-1
    
    # make array of the middle of each bin; len(mid_bins) = len(bins) - 1
    mid_bins = (bins[:-1] + bins[1:])/2.

    # create empty arrays to store mean concentration as a function of wind direction
    gasConc_avg = np.zeros_like(mid_bins)

    # for each wind direction bin, find the corresponding mean tracer concentration
    for i in range(indices.min(),indices.max()+1):
        temp_vals = gasConc[~ws_mask][np.where(indices == i)]

        # ensure that the amount of data points exceeds the cutoff value
        # we don't want bins with a couple of values to dominate the fitting
        
        if len(temp_vals) > int(cutoff*len(gasConc)/100.):
            gasConc_avg[i] = temp_vals.mean()

    # ensure that the peak concentration is around 180 degrees for fitting
    roll_amount = int(len(gasConc_avg)/2. -3) - np.argmin(abs(gasConc_avg - np.average(wd3[~ws_mask],weights=gasConc[~ws_mask])))#np.where(tracer_avg == tracer_avg.max())[0][0]
    gasConc_avg = np.roll(gasConc_avg,roll_amount)
    
    # get the bin with peak average concentration
    max_bin = mid_bins[np.where(gasConc_avg == gasConc_avg.max())[0][0]]
    
    # calculate wind direction cut off based on input value
    bin_cut_lo = max_bin - wdlimit
    bin_cut_hi = max_bin + wdlimit

    # add the roll amount to wind direction to correspond correctly to the
    # average concentration array
    wd3 = wd3 + (roll_amount-1)*delta_theta

    # ensure that wind direction is within 0 to 360
    wd3 = wrap(wd3)

    # get the indices of where wd3 is outside of wind direction limit
    wd_mask = np.where((wd3 > bin_cut_hi) | (wd3 < bin_cut_lo))
    
    # mask arrays where wd3 is outside of wind direction cutoff
    wd3[wd_mask] = np.ma.masked
    ws3[wd_mask] = np.ma.masked
    ws3z[wd_mask] = np.ma.masked
    temp[wd_mask] = np.ma.masked
    pres[wd_mask] = np.ma.masked

    # fitting procedure
    # here are some initial guesses that produce good results
    const_0 = [gasConc_avg.max(),24,180]
    
    # the curve fit function uses the Levenberg-Marquardt algorithm which
    # in my opionion does a great job
    fit_gasConc,cov_gasConc = curve_fit(gaussian_func,mid_bins,gasConc_avg,p0 = const_0) # fit coefficients
    if make_plot: fit_plot(mid_bins,gasConc_avg,fit_gasConc,chemical_name) # make plot if you want

    # calculate the standard deviation of wind direction and turbulent intensity 
    # for use in finding the PG stability class
    turbulent_intensity = np.std(ws3z)/np.mean(ws3) # turbulent intensity
    std_wind_dir = yamartino_method(wd3) # st. dev. of wind direction [deg]

    # calcualte the vertical and horizontal dispersion of the gaussian plume
    # using PGT stability classes
    sy, sz = sigma(distance,std_wind = std_wind_dir,turb = turbulent_intensity,tables=True,stab=None) # [m]

    # convert tracers from ppm to g m-3
    fit_amplitude = ppm2gm3(fit_gasConc[0],mw_chemical,np.mean(temp),np.mean(pres))

    # calculate the density of the tracers using the ideal gas law for conversion back to L/min if desired
    rho_gasConc_stp = (P_stp*mw_chemical)/(R*T_stp) # density of gas at STP [g L-1]
    rho_gasConc = (np.mean(pres)*mw_chemical)/(R*np.mean(temp)) # density of gas at ambient conditions [g L-1]
   
    # putting it all together to calculate the emission rate in g/s
    emission_gasConc_mass_per_time = 2*np.pi*fit_amplitude*np.mean(ws3)*sy*sz # [g s-1]

    # now calculate the emission rate in L/min
    emission_gasConc_volume_per_time = 2*np.pi*fit_amplitude*np.mean(ws3)*sy*sz/rho_gasConc*60 # [L min-1]

    # now L/min at STP
    emission_gasConc_volume_per_time_stp = 2*np.pi*fit_amplitude*np.mean(ws3)*sy*sz/rho_gasConc_stp*60 # [L min-1]
    return emission_gasConc_volume_per_time, emission_gasConc_mass_per_time
