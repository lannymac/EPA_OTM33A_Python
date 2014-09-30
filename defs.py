import numpy as np
from defs import *
import pylab as pl

def func(x,a,sigma,x_0):
    '''
    This function will return a gaussian with the
    proper coefficients
    '''
    return a*np.exp(-((x - x_0)**2)/(2*sigma**2))

def stability_class(std_wind,turb):
    '''
    This function will return an average stability class
    based on the standard deviation of the wind speed 
    and the turbulent intensity.

    INPUTS
    std_wind :: standard deviation on the horizontal wind direction [degrees]
    turb     :: turbulent intensity
    '''

    if std_wind > 27.5:
        pg_wind = 1.
    elif std_wind <=27.5 and std_wind > 23.5:
        pg_wind = 2.
    elif std_wind <=23.5 and std_wind >19.5:
        pg_wind = 3.
    elif std_wind <= 19.5 and std_wind > 15.5:
        pg_wind = 4.
    elif std_wind <= 15.5 and std_wind > 11.5:
        pg_wind = 5.
    elif std_wind <= 11.5 and std_wind > 7.5:
        pg_wind = 6.
    elif std_wind <= 7.50:
        pg_wind = 7.

    if turb > 0.205:
        pg_turb = 1.
    elif turb <= .205 and turb > 0.180:
        pg_turb = 2.
    elif turb <= 0.180 and turb > 0.155:
        pg_turb = 3.
    elif turb <=0.155 and turb > 0.130:
        pg_turb = 4.
    elif turb <= 0.130 and turb > 0.105:
        pg_turb = 5.
    elif turb <= 0.105 and turb > 0.080:
        pg_turb = 6.
    elif turb <= 0.080:
        pg_turb = 7.

    # average together the two estimates of stability class and round
    final_pg = int(np.round(np.mean((pg_wind,pg_turb))))
    return final_pg


def sigma_func(I,J,K,dist):
    '''
    This function is used to estimate sigma_y and sigma_z based on 
    stability parameters.
    '''
    return np.exp(I + J*np.log(dist) + K*(np.log(dist)**2))
    
def sigma(dist,std_wind=None,turb=None,stab= None):
    
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

    if stab == 1:
        sigma_y = sigma_func(-1.104,.9879,-.0076,dist)
        sigma_z = sigma_func(4.679,-1.7172,0.2770,dist)
    elif stab == 2:
        sigma_y = np.mean([sigma_func(-1.104,.9879,-.0076,dist),sigma_func(-1.634,1.0350,-0.0096,dist)],axis=0)
        sigma_z = np.mean([sigma_func(-1.999,0.8752,0.0136,dist),sigma_func(4.679,-1.7172,0.2770,dist)],axis=0)
    elif stab == 3:
        sigma_y = sigma_func(-1.634,1.0350,-0.0096,dist)
        sigma_z = sigma_func(-1.999,0.8752,0.0136,dist)
    elif stab == 4:
        sigma_y = np.mean([sigma_func(-1.634,1.0350,-0.0096,dist),sigma_func(-2.054,1.0231,-0.0076,dist)],axis=0)
        sigma_z = np.mean([sigma_func(-1.999,0.8752,0.0136,dist),sigma_func(-2.341,0.9477,-0.0020,dist)],axis=0)
    elif stab == 5:
        sigma_y = sigma_func(-2.054,1.0231,-0.0076,dist)
        sigma_z = sigma_func(-2.341,0.9477,-0.0020,dist)
    elif stab == 6:
        sigma_y = np.mean([sigma_func(-2.054,1.0231,-0.0076,dist),sigma_func(-2.555,1.0423,-0.0087,dist)],axis=0)
        sigma_z = np.mean([sigma_func(-2.341,0.9477,-0.0020,dist),sigma_func(-3.186,1.1737,-0.0316,dist)],axis=0)
    elif stab == 7:
        sigma_y = sigma_func(-2.555,1.0423,-0.0087,dist)
        sigma_z = sigma_func(-3.186,1.1737,-0.0316,dist)
    else:
        sigma_y = np.nan
        sigma_z = np.nan

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
    '''
    y_fit = func(x,fit[0],fit[1],fit[2])
    
    fs=18
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y_data,c='k',label='Data',lw=2)
    ax.plot(x,y_fit,c='r',label='Fit',lw=2)    
    ax.set_xlabel('Wind Direction [degrees]',fontsize=fs)
    ax.set_ylabel('Average concentration [g/m3]',fontsize=fs)
    ax.set_title(name,fontsize=fs+2)
    leg=ax.legend(fontsize=fs)
    leg.get_frame().set_alpha(0.)
    pl.show()


def print_filename(filename):
    '''
    This definition will print a header to the screen
    '''
    print('+'*(len(filename)+10))
    print('+    %s    +' % (filename))
    print('+'*(len(filename)+10)+'\n')



def print_screen(name,emis,emis_stp):
    '''
    This definition will print results from the program
    '''
    print('%s is predicted to have an emission rate of %.3f LPM' % (name,emis))
    print('%s is predicted to have an emission rate of %.3f SLPM\n' % (name,emis_stp))
    


def custom_load_files(filename,chemical_1,chemical_2):
    '''
    This is a custom file loader that I made since
    I store my files in npz format. 
    '''
    fuck = np.load(filename)
    tracer_1 = np.ma.array(fuck[chemical_1],mask=False)
    tracer_2 = np.ma.array(fuck[chemical_2],mask=False)/1000.
    lat = np.ma.array(fuck['lat'],mask=False)
    lon = np.ma.array(fuck['lon'],mask=False)
    ws3 = np.ma.array(fuck['ws3'],mask=False)
    wd3 = np.ma.array(fuck['wd3'],mask=False)
    ws2 = np.ma.array(fuck['ws2'],mask=False)
    wd2 = np.ma.array(fuck['wd2'],mask=False)
    temp = np.ma.array(fuck['temp'],mask=False)+273.
    pres = np.ma.array(fuck['pres'],mask=False)/1000.
    ws3z = np.ma.array(fuck['ws3z'],mask=False)
    ws3x = np.ma.array(fuck['ws3x'],mask=False)
    ws3y = np.ma.array(fuck['ws3y'],mask=False)
    ws3t = np.ma.array(fuck['ws3t'],mask=False)
    time = np.ma.array(fuck['time'],mask=False)
    
    return tracer_1,tracer_2,lat,lon,ws3,wd3,ws2,wd2,temp,pres,ws3z,ws3x,ws3y,ws3t,time
