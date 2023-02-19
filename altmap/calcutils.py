#standard libs
import pkgutil
import pathlib
import sys

#external libs
import numpy as np
from scipy import interpolate
import ndjson

def getdatapath():
    package = pkgutil.get_loader("altmap")
    datapath = str(pathlib.Path(package.get_filename()).parent)+"/data/"

    return datapath

def getmagdec(lat_c,lon_c,searchthres=0.33):

    path = getdatapath()+'gsigeomag2020_dec_ver1.0_20220111.csv'
    temp = np.loadtxt(path,delimiter=",")

    def diffang(lat,lon):
        return ((lat-lat_c)**2 + (lon-lon_c)**2 )**0.5

    res = [[poi[2]+poi[3]/60.,diffang(poi[0],poi[1])] for poi in temp if diffang(poi[0],poi[1]) < searchthres]
    
    return res[np.array(res)[:,1].argmin()][0] #[deg]

def getPrcslist():

    """
    ref: https://www.gsi.go.jp/LAW/heimencho.html
    """

    path = getdatapath()+"prcs_h14.ndjson"

    with open(path,"r") as f:
        input_data = ndjson.load(f)
    
    return input_data

def getPrcsorigin(system_num):

    temp = [ prcsdata for prcsdata in getPrcslist() if prcsdata['system_num']==system_num ]

    if len(temp)==1:
        lat = float(temp[0]["lat_deg"]) + float(temp[0]["lat_min"])/60. + float(temp[0]["lat_sec"])/3600.
        lon = float(temp[0]["lon_deg"]) + float(temp[0]["lon_min"])/60. + float(temp[0]["lon_sec"])/3600.
        
        return lat, lon #[deg]

    else:
        print('Error: "system_num" may be out of range.', file=sys.stderr)
        sys.exit(1)

def latlon2xy(phi_deg, lambda_deg, phi0_deg, lambda0_deg):

    """
    ref: https://qiita.com/sw1227/items/e7a590994ad7dcd0e8ab
    """
    
    phi_rad = np.deg2rad(phi_deg)
    lambda_rad = np.deg2rad(lambda_deg)
    phi0_rad = np.deg2rad(phi0_deg)
    lambda0_rad = np.deg2rad(lambda0_deg)

    def A_array(n):
        A0 = 1 + (n**2)/4. + (n**4)/64.
        A1 = -     (3./2)*( n - (n**3)/8. - (n**5)/64. ) 
        A2 =     (15./16)*( n**2 - (n**4)/4. )
        A3 = -   (35./48)*( n**3 - (5./16)*(n**5) )
        A4 =   (315./512)*( n**4 )
        A5 = -(693./1280)*( n**5 )
        return np.array([A0, A1, A2, A3, A4, A5])

    def alpha_array(n):
        a0 = np.nan # dummy
        a1 = (1./2)*n - (2./3)*(n**2) + (5./16)*(n**3) + (41./180)*(n**4) - (127./288)*(n**5)
        a2 = (13./48)*(n**2) - (3./5)*(n**3) + (557./1440)*(n**4) + (281./630)*(n**5)
        a3 = (61./240)*(n**3) - (103./140)*(n**4) + (15061./26880)*(n**5)
        a4 = (49561./161280)*(n**4) - (179./168)*(n**5)
        a5 = (34729./80640)*(n**5)
        return np.array([a0, a1, a2, a3, a4, a5])

    m0 = 0.9999 
    a = 6378137.
    F = 298.257222101

    n = 1. / (2*F - 1)
    A_array = A_array(n)
    alpha_array = alpha_array(n)

    A_ = ( (m0*a)/(1.+n) )*A_array[0] # [m]
    S_ = ( (m0*a)/(1.+n) )*( A_array[0]*phi0_rad + np.dot(A_array[1:], np.sin(2*phi0_rad*np.arange(1,6))) ) # [m]

    lambda_c = np.cos(lambda_rad - lambda0_rad)
    lambda_s = np.sin(lambda_rad - lambda0_rad)

    t = np.sinh( np.arctanh(np.sin(phi_rad)) - ((2*np.sqrt(n)) / (1+n))*np.arctanh(((2*np.sqrt(n)) / (1+n)) * np.sin(phi_rad)) )
    t_ = np.sqrt(1 + t*t)

    xi2  = np.arctan(t / lambda_c) # [rad]
    eta2 = np.arctanh(lambda_s / t_)

    x = A_ * (xi2 + np.sum(np.multiply(alpha_array[1:],
                                       np.multiply(np.sin(2*xi2*np.arange(1,6)),
                                                   np.cosh(2*eta2*np.arange(1,6)))))) - S_ # [m]
    y = A_ * (eta2 + np.sum(np.multiply(alpha_array[1:],
                                        np.multiply(np.cos(2*xi2*np.arange(1,6)),
                                                    np.sinh(2*eta2*np.arange(1,6)))))) # [m]
    
    return x, y # [m]


def xy2latlon(x, y, phi0_deg, lambda0_deg):

    """
    ref: https://qiita.com/sw1227/items/e7a590994ad7dcd0e8ab
    """
    
    phi0_rad = np.deg2rad(phi0_deg)
    lambda0_rad = np.deg2rad(lambda0_deg)

    def A_array(n):
        A0 = 1 + (n**2)/4. + (n**4)/64.
        A1 = -     (3./2)*( n - (n**3)/8. - (n**5)/64. ) 
        A2 =     (15./16)*( n**2 - (n**4)/4. )
        A3 = -   (35./48)*( n**3 - (5./16)*(n**5) )
        A4 =   (315./512)*( n**4 )
        A5 = -(693./1280)*( n**5 )
        return np.array([A0, A1, A2, A3, A4, A5])

    def beta_array(n):
        b0 = np.nan # dummy
        b1 = (1./2)*n - (2./3)*(n**2) + (37./96)*(n**3) - (1./360)*(n**4) - (81./512)*(n**5)
        b2 = (1./48)*(n**2) + (1./15)*(n**3) - (437./1440)*(n**4) + (46./105)*(n**5)
        b3 = (17./480)*(n**3) - (37./840)*(n**4) - (209./4480)*(n**5)
        b4 = (4397./161280)*(n**4) - (11./504)*(n**5)
        b5 = (4583./161280)*(n**5)
        return np.array([b0, b1, b2, b3, b4, b5])

    def delta_array(n):
        d0 = np.nan # dummy
        d1 = 2.*n - (2./3)*(n**2) - 2.*(n**3) + (116./45)*(n**4) + (26./45)*(n**5) - (2854./675)*(n**6)
        d2 = (7./3)*(n**2) - (8./5)*(n**3) - (227./45)*(n**4) + (2704./315)*(n**5) + (2323./945)*(n**6)
        d3 = (56./15)*(n**3) - (136./35)*(n**4) - (1262./105)*(n**5) + (73814./2835)*(n**6)
        d4 = (4279./630)*(n**4) - (332./35)*(n**5) - (399572./14175)*(n**6)
        d5 = (4174./315)*(n**5) - (144838./6237)*(n**6)
        d6 = (601676./22275)*(n**6)
        return np.array([d0, d1, d2, d3, d4, d5, d6])

    m0 = 0.9999 
    a = 6378137.
    F = 298.257222101

    n = 1. / (2*F - 1)
    A_array = A_array(n)
    beta_array = beta_array(n)
    delta_array = delta_array(n)

    A_ = ( (m0*a)/(1.+n) )*A_array[0]
    S_ = ( (m0*a)/(1.+n) )*( A_array[0]*phi0_rad + np.dot(A_array[1:], np.sin(2*phi0_rad*np.arange(1,6))) )

    xi = (x + S_) / A_
    eta = y / A_

    xi2 = xi - np.sum(np.multiply(beta_array[1:], 
                                  np.multiply(np.sin(2*xi*np.arange(1,6)),
                                              np.cosh(2*eta*np.arange(1,6)))))
    eta2 = eta - np.sum(np.multiply(beta_array[1:],
                                   np.multiply(np.cos(2*xi*np.arange(1,6)),
                                               np.sinh(2*eta*np.arange(1,6)))))

    chi = np.arcsin( np.sin(xi2)/np.cosh(eta2) ) # [rad]
    latitude = chi + np.dot(delta_array[1:], np.sin(2*chi*np.arange(1, 7))) # [rad]

    longitude = lambda0_rad + np.arctan( np.sinh(eta2)/np.cos(xi2) ) # [rad]

    return np.rad2deg(latitude), np.rad2deg(longitude) # [deg]


def getZfactor(center_Lat):

    """
    Note: Only in case the unit of X,Y is decimal degrees AND the unit of meters.
    ref: https://pro.arcgis.com/ja/pro-app/latest/tool-reference/3d-analyst/applying-a-z-factor.htm
    """
    zfacters = np.array([
        [0, 0.00000898],
        [10,0.00000912],
        [20,0.00000956],
        [30,0.00001036],
        [40,0.00001171],
        [50,0.00001395],
        [60,0.00001792],
        [70,0.00002619],
        [80,0.00005156]
        ])
    
    fitfunc = fitfunc = interpolate.interp1d(zfacters[:,0],zfacters[:,1],fill_value='extrapolate')

    return fitfunc(center_Lat)

def calcHillshade(Z, center_Lat, azimuth_deg, altitude_deg):

    """
    ref: https://memomemokun.hateblo.jp/entry/2018/12/22/205419
    ref: https://pro.arcgis.com/ja/pro-app/latest/tool-reference/3d-analyst/how-hillshade-works.htm
    """

    # elevation mesh
    Zx, Zy = np.gradient(Z)
    Slope_rad = np.arctan(500 * getZfactor(center_Lat) * np.sqrt(Zx*Zx + Zy*Zy))
    Aspect_rad = np.arctan2(-Zx, Zy)

    # Sun
    Azimuth_rad = np.radians(360.0 - azimuth_deg + 90)
    Zenith_rad = np.radians(90 - altitude_deg)

    # hillshade
    shaded = np.cos(Zenith_rad) * np.cos(Slope_rad) + np.sin(Zenith_rad) * np.sin(Slope_rad) * np.cos(Azimuth_rad - Aspect_rad)
    
    return 255 * shaded

def tile2latlon(x,y,z):
    
    mapy = (y/2.0**z) * 2 * np.pi - np.pi
    lat = 2.*np.atan(np.exp(- mapy)) * 180./np.pi - 90.
    lon = (x/2.0**z) * 360. - 180.
    
    return lat,lon


def latlon2tile(lat,lon,z):
    
    x = int((lon/180.+1.) * 2.**z/2.)
    y = int(((-np.log(np.tan((45.+lat/2.) * np.pi/180.)) + np.pi) * 2.**z / (2.*np.pi)))
    
    return x,y