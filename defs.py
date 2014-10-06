import numpy as np
import pylab as pl
import xlrd
from scipy import optimize
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


