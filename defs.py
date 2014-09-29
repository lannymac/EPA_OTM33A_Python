import numpy as np
from defs import *
def func(x,a,sigma,x_0):
    '''
    This function will return a gaussian with the
    proper coefficients
    '''
    return a*np.exp(-((x - x_0)**2)/(2*sigma**2))

def stability_class(std_wind,turb):
    
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
    elif std_wind <= 11.5:
        pg_wind = 6.

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
    elif turb <= 0.105:
        pg_turb = 6.
    final_pg = np.mean((pg_wind,pg_turb))
    
    final_pg_letter = num2let(int(np.round(final_pg)))
    return final_pg_letter


def num2let(number):
    if number == 1:
        letter = 'A'
    elif number == 2:
        letter = 'B'
    elif number == 3:
        letter = 'C'
    elif number == 4:
        letter = 'D'
    elif number == 5:
        letter = 'E'
    elif number == 6:
        letter = 'F'
    return letter

def sigma_func(I,J,K,dist):
    return np.exp(I + J*np.log(dist) + K*(np.log(dist)**2))
    
def sigma_y(dist,std_wind,turb):

    stab = stability_class(std_wind,turb)

    if stab == 'A':
        sigma = sigma_func(-1.104,.9879,-.0076,dist)
    elif stab == 'B':
        sigma = sigma_func(-1.634,1.0350,-0.0096,dist)
    elif stab == 'C':
        sigma = sigma_func(-2.054,1.0231,-0.0076,dist)
    elif stab == 'D':
        sigma = sigma_func(-2.555,1.0423,-0.0087,dist)
    elif stab == 'E':
        sigma = sigma_func(-2.754,1.0106,-0.0064,dist)
    elif stab == 'F':
        sigma = sigma_func(-3.143,1.0148,-0.0070,dist)
    else:
        sigma = np.nan
    return sigma

def sigma_z(dist,std_wind,turb):

    stab = stability_class(std_wind,turb)
    
    if stab == 'A':
        sigma = sigma_func(4.679,-1.7172,0.2770,dist)
    elif stab == 'B':
        sigma = sigma_func(-1.999,0.8752,0.0136,dist)
    elif stab == 'C':
        sigma = sigma_func(-2.341,0.9477,-0.0020,dist)
    elif stab == 'D':
        sigma = sigma_func(-3.186,1.1737,-0.0316,dist)
    elif stab == 'E':
        sigma = sigma_func(-3.783,1.3010,-0.0450,dist)
    elif stab == 'F':
        sigma = sigma_func(-4.490,1.4024,-0.0540,dist)
    else:
        sigma = np.nan
    return sigma

