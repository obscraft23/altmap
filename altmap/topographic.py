import asyncio
import json

import numpy as np
import aiohttp
from scipy import interpolate

import matplotlib.style as mplstyle
mplstyle.use('fast')
import matplotlib.pyplot as plt
import japanize_matplotlib
from matplotlib.colors import Normalize

from altmap.calcutils import xy2latlon, getmagdec, latlon2xy

#lat0,lon0 = 35.42749067885798, 139.0390205383301
#lat1,lon1 = 35.44364499417055, 139.05309677124026

def _latlon2pq(lat,lon,z=18,L=85.05112827):
    p = 2**(z+7) * (lon/180. + 1)
    q = 2**(z+7)/np.pi * (-np.arctanh(np.sin(lat*np.pi/180.)) + np.arctanh(np.sin(L*np.pi/180.)))

    return int(p/256),int(q/256)

def _getXYlist(lat0,lon0,lat1,lon1):
    p0, q0 = _latlon2pq(lat0,lon0)
    p1, q1 = _latlon2pq(lat1,lon1)
    x_list = np.linspace(p0,p1,(p1-p0)+1).astype('int')
    y_list = np.linspace(q1,q0,(q0-q1)+1).astype('int')

    return x_list,y_list

async def _get_geojson(session,url):
    async with session.get(url) as resp:
        if resp.status == 200:
            res = await resp.read()
            geojsondata = json.loads(res)
            return [poi['geometry']['coordinates'] + [poi['properties']["alti"]] for poi in geojsondata["features"]]
        else:
            return []

async def getGeojson(x_list,y_list):
    points_3d = []

    async with aiohttp.ClientSession() as session:
        tasks_5a = []
        tasks_10b = []
        for x in x_list:
            for y in y_list:
                url_5a = f'https://cyberjapandata.gsi.go.jp/xyz/experimental_dem5a/18/{x}/{y}.geojson'
                tasks_5a.append(asyncio.ensure_future(_get_geojson(session, url_5a)))
                url_10b = f'https://cyberjapandata.gsi.go.jp/xyz/experimental_dem10b/18/{x}/{y}.geojson'
                tasks_10b.append(asyncio.ensure_future(_get_geojson(session, url_10b)))

        plist_5a = await asyncio.gather(*tasks_5a)
        plist_10b = await asyncio.gather(*tasks_10b)

    for pp_5a, pp_10b in zip(plist_5a,plist_10b):
        if pp_5a != []:
            points_3d += pp_5a
        else:
            points_3d += pp_10b
    
    return points_3d

def createTopographic(lat0,lon0,lat1,lon1,magnetic_north_line=True, max_length=2000):

    y_min0,x_min0=latlon2xy(lat0,lon0,(lat0+lat1)/2.,(lon0+lon1)/2.)
    y_max0,x_max0=latlon2xy(lat1,lon1,(lat0+lat1)/2.,(lon0+lon1)/2.)

    if (x_max0 - x_min0) > max_length:
        return None

    x_list,y_list = _getXYlist(lat0,lon0,lat1,lon1)
    plist = asyncio.run(getGeojson(x_list,y_list))
    #plist = await getGeojson(x_list,y_list)
    
    res = np.array(plist)
    nx = round((res[:,0].max() - res[:,0].min()) * 3600 * 5)
    ny = round((res[:,1].max() - res[:,1].min()) * 3600 * 5)
    nf=5

    lon_min = res[:,0].min()
    lon_max = res[:,0].min() + 0.2/3600.*(nx)
    lat_min = res[:,1].min()
    lat_max = res[:,1].min() + 0.2/3600.*(ny)
    
    y_min,x_min=latlon2xy(lat_min,lon_min,(lat_min+lat_max)/2.,(lon_min+lon_max)/2)
    y_max,x_max=latlon2xy(lat_max,lon_max,(lat_min+lat_max)/2.,(lon_min+lon_max)/2)

    x_coord = np.linspace(x_min,x_max,nf*nx+1)
    y_coord = np.linspace(y_min,y_max,nf*ny+1)
    _x_coord = np.linspace(lon_min,lon_max,nf*nx+1)
    _y_coord = np.linspace(lat_min,lat_max,nf*ny+1)
    __xx, __yy = np.meshgrid(_x_coord,_y_coord)
    
    _xx, _yy = np.meshgrid(x_coord, y_coord)
    
    mes = interpolate.griddata(points=res[:,0:2], values=res[:,2], xi=(__xx, __yy), method='cubic')

    Ex,Ey = np.gradient(mes,1.,1.)
    grad = np.sqrt(Ex**2 + Ey**2)
    grad[np.isnan(grad)]=0

    xyratio = (y_max0 - y_min0) / (x_max0 - x_min0)

    fig = plt.figure(figsize=(1, xyratio), dpi=2*round(x_max0 - x_min0))
    ax1 = fig.add_subplot(111)
    
    ax1.set_xlim((x_min0,x_max0))
    ax1.set_ylim((y_min0,y_max0))

    ax1.pcolormesh(_xx, _yy, mes, norm=Normalize(vmin=res[:,2].min(), vmax=res[:,2].max()), cmap='gray', shading='auto',alpha=.8)
    ax1.pcolormesh(_xx, _yy, grad, norm=Normalize(vmin=grad[grad>0].min(), vmax=min(2,grad[grad>0].max())), cmap='Reds', shading='auto',alpha=.3)

    sep = 5
    nn = int((res[:,2].max()-res[:,2].min()) / sep)
    ax1.contour(_xx, _yy, mes, levels=np.linspace(res[:,2].min(),sep*nn+res[:,2].min(),nn+1), colors=['black'], linewidths=0.0001,antialiased=False)
    ax1.contour(_xx, _yy, mes, levels=np.linspace(100,3800,75), colors=['black'], linewidths=0.1,antialiased=False)

    if magnetic_north_line:
        lat_c = (lat0+lat1)/2
        lon_c = (lon0+lon1)/2
        magrad = getmagdec(lat_c,lon_c) * np.pi/180.
        xx = np.linspace(x_min0,x_max0,100)
        max_maglines = int(x_max0/500)+1
        for dx in np.linspace(-max_maglines*500,max_maglines*500,2*max_maglines+1):
            ax1.plot(xx,np.tan(-np.pi/2+magrad)*(xx-lon_c-dx)+lat_c,c="blue",lw=0.1)

    """
    texts = [ ax1.text(float(p["lon"]), float(p["lat"]), p["name"], horizontalalignment='center', verticalalignment='center',fontsize=1) for p in tmp]

    [ax1.plot(np.array(route["coords"])[:,0],np.array(route["coords"])[:,1],c="r",lw=.3) for route in tqdm(data)]
    """

    """
    lat_c = (lat0+lat1)/2
    lon_c = (lon0+lon1)/2

    xx = np.linspace(res[:,0].min(),res[:,0].max(),100)
    magrad = getmagdec(lat_c,lon_c) * np.pi/180.
    for dx in [-1000,-500,0,500,1000]:
        y0, x0 = xy2latlon(0,dx,lat_c, lon_c)
        ax1.plot(xx,np.tan(-np.pi/2+magrad)*(xx-x0)+y0,c="gray",lw=0.2)
    """

    ax1.set_aspect('equal')
    ax1.axis("off")
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    
    return fig