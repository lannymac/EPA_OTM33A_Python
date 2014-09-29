import numpy as np
from main import *
from scipy.optimize import curve_fit
from defs import *
import sys

files = ['CHRISTMAN_2014_03_14_1010_TO_1034.npz']

theta_start = 0.
theta_end = 360.
delta_theta = 10.

d = [40.]

for k1 in range(len(files)):

    fuck = np.load(files[k1])
    tracer = np.ma.array(fuck['C2H2'],mask=False)
    ch4 = np.ma.array(fuck['CH4'],mask=False)
    lat = np.ma.array(fuck['lat'],mask=False)
    lon = np.ma.array(fuck['lon'],mask=False)
    ws3 = np.ma.array(fuck['ws3'],mask=False)
    wd3 = np.ma.array(fuck['wd3'],mask=False)
    ws2 = np.ma.array(fuck['ws2'],mask=False)
    wd2 = np.ma.array(fuck['wd2'],mask=False)
    temp = np.ma.array(fuck['temp'],mask=False)
    pres = np.ma.array(fuck['pres'],mask=False)
    ws3z = np.ma.array(fuck['ws3z'],mask=False)
    ws3x = np.ma.array(fuck['ws3x'],mask=False)
    ws3y = np.ma.array(fuck['ws3y'],mask=False)
    ws3t = np.ma.array(fuck['ws3t'],mask=False)
    time = np.ma.array(fuck['time'],mask=False)


    mean_wind_speed = np.mean(ws3)
    turbulent_intensity = np.std(ws3z)/mean_wind_speed

    #subtract ch4 background
    ch4 = ch4 - np.percentile(ch4,05)

    # get bins
    bins = np.arange(theta_start,theta_end+1,delta_theta)
    indices= np.digitize(wd3,bins)-1
    mid_bins = (bins[:-1] + bins[1:])/2.
    hist = np.zeros_like(mid_bins)
    
    for i in range(indices.min(),indices.max()+1):
        hist[i] = ch4[np.where(indices == i)].mean()
        
    roll_amount = int(len(hist)/2. -1) - np.where(hist == hist.max())[0][0]
    hist = np.roll(hist,roll_amount)

    # fitting procedure
    const_0 = [ch4.max(),24,180]
    fit,cov = curve_fit(func,mid_bins,hist,p0 = const_0)


    sy = sigma_y(d[k1],np.std(ws3*180/np.pi),turbulent_intensity)
    sz = sigma_z(d[k1],np.std(ws3*180/np.pi),turbulent_intensity)
    
    mw = 16.0
    gasconst = 8.314510   # J mol-1 K-1
    rt0rp0 = (gasconst*298)/101.325
    emission_ch4 = 2*np.pi*fit[0]*mean_wind_speed*sy*(1e-06 * mw*1000)/rt0rp0
