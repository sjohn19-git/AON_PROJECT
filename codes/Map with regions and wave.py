import math
from scipy.integrate import trapz
from os.path import isfile, join
import obspy.signal.cross_correlation as cr
from shapely.geometry import Point, Polygon
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
import matplotlib.dates as mdates
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from medfilt import medfilt
import cdsapi
import matplotlib.pyplot as plt
import xarray as xr
import pygrib
import pickle
import numpy as np
import glob
import cv2
from obspy import UTCDateTime
import pandas as pd
import pygmt
import os
os.chdir("/Users/sebinjohn/AON_PROJECT/Data/codes")
# import dill


os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
with open("metadta.pkl", "rb") as f:
   long, lat, stationo, env = pickle.load(f)

datapath = "/Users/sebinjohn/AON_PROJECT/Data/*/"
st = UTCDateTime(2019, 11, 1)
et = UTCDateTime(2019, 12, 1)
windw = 1
tim_len = int((et-st)/(3600*windw))
smth_intg = np.zeros((len(stationo), tim_len))


freq = []
name = pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range(95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)

time_frame = []
for i in range(int((et-st)/(3600*windw))):
    time_frame.append(st)
    st = st+(3600*windw)
os.chdir("/Users/sebinjohn/AON_PROJECT/Data/median_Power_time_series")
# Python 3: open(..., 'wb')
with open(("median and mean"+str(time_frame[0])+"-"+str(time_frame[-1])+"Secondary.pkl"), 'rb') as f:
          mean_timeseries, median_timeseries = pickle.load(f)


def sort_counterclockwise(points, centre=None):
  if centre:
    centre_x, centre_y = centre
  else:
    centre_x, centre_y = sum([x for x, _ in points]) / \
                             len(points), sum(
                                 [y for _, y in points])/len(points)
  angles = [math.atan2(y - centre_y, x - centre_x) for x, y in points]
  counterclockwise_indices = sorted(
      range(len(points)), key=lambda i: angles[i])
  counterclockwise_points = [points[i] for i in counterclockwise_indices]
  return counterclockwise_points


def poly_seis(median_timeseries, mean_timeseries, *args):
    coords = []
    for point in args:
        coords.append(point)
    counter = sort_counterclockwise(coords, centre=None)
    poly = Polygon(counter)
    wista = []
    for i in range(np.shape(median_timeseries)[1]):
        p1 = Point(median_timeseries[1, i, 5], median_timeseries[2, i, 5])
        if p1.within(poly):
            wista.append(1)
        else:
            wista.append(0)
    c = wista.count(1)
    j = 0
    sum_array = np.empty([c, 2, len(median_timeseries[0, 0, 1:])])
    for i in range(len(wista)):
        if wista[i] == 1:
            sum_array[j, 0, :] = median_timeseries[0, i, 1:]
            sum_array[j, 1, :] = mean_timeseries[0, i, 1:]
            j += 1
    sum_tot = np.empty([2, len(median_timeseries[0, 0, 1:])])
    for i in range(len(sum_array[0, 0, :])):
        sum_tot[0, i] = np.sum(sum_array[:, 0, i]) / \
                               (c-np.count_nonzero(sum_array[:, 0, i] == 0))
        sum_tot[1, i] = np.sum(sum_array[:, 1, i]) / \
                               (c-np.count_nonzero(sum_array[:, 1, i] == 0))
    return wista, poly, counter, sum_array, sum_tot


def cross_correlate(time_frame, sum_tot, itera, *argv):
    c = 0
    if itera == 1:
        for i in range(len(time_frame)):
            tim = time_frame[i]
            print("searching for wave height for"+str(tim))
            grid, data_wave = wave(tim)
            # sum_w=wave_val(grid,wlat1,wlat2,wlon1,wlon2)
            wave_arr, counterw, lat, lon = wave_array(
                grid, (150, 30), (260, 73), (260, 30), (150, 73))
            if c == 0:
               a, b = np.shape(wave_arr.data)
               wave_all = np.zeros([a, b, len(time_frame)])
               c = 1
            wave_all[:, :, i] = wave_arr.data
    else:
        for i in range(1):
            tim = time_frame[i]
            print("searching for wave height for"+str(tim))
            grid, data_wave = wave(tim)
            # sum_w=wave_val(grid,wlat1,wlat2,wlon1,wlon2)
            wave_arr, counterw, lat, lon = wave_array(
                grid, (150, 30), (260, 73), (260, 30), (150, 73))
        wave_all = argv[0]
    corre_grid = np.zeros([np.shape(wave_all)[0], np.shape(wave_all)[1]])
    corre_inde = np.zeros(np.shape(corre_grid))
    for i in range(np.shape(wave_all)[0]):
        for j in range(np.shape(wave_all)[1]):
            k = cr.correlate(wave_all[i, j, :], sum_tot[0, :], 720)
            corre_grid[i, j] = k[720]
            corre_inde[i, j] = cr.xcorr_max(k)[0]
    return corre_grid, corre_inde, wave_all, counterw, lat, lon


regions = [[(-151.95, 70.3), (-151.6, 66.5), (-141.6, 66.5), (-141.88, 69.4), "north_east"], [(-165.9, 70), (-166.1, 66.5), (-154.3, 66.5), (-153.88, 70), "North_west"], [(-157.1, 59), (-152.8, 57), (-168.4, 52.3), (-170.4, 54.2), "south_west"],
         [(-140, 62), (-140.5, 59.5), (-131.5, 55), (-129, 56), "south_east"], [(-166.322, 62.951), (-153.788, 63.587), [-153.820, 60.768], (-163.582, 59.649), "central_west"], [(-149.853, 64.549), (-141.512, 64.760), (-142.462, 60.537), (-149.370, 60.446), "central_east"]]


def wave(tim):
    c = cdsapi.Client()
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/wave")
    files = os.listdir()
    if str(tim)[0:14]+"00"+".grib" not in files:
        print(str(tim)[0:14]+"00"+".grib not found Downloading...")
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': 'significant_height_of_combined_wind_waves_and_swell',
                'year': str(tim.year),
                'month': str(tim.month),
                'day': str(tim.day),
                'time': tim.ctime()[11:14]+"00",
                'format': 'grib',
                },
            str(tim)[0:14]+"00"+".grib")
    tim_wav = str(tim)[0:14]+"00"
    grib = str(tim)[0:14]+"00"+".grib"
    grbs = pygrib.open(grib)
    grb = grbs[1]
    data_wave = grb.values
    latg, long = grb.latlons()
    lat = (latg[:, 0])
    lon = (long[0, :])
    grid = xr.DataArray(
        data_wave, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
        )
    return grid, data_wave
 grid,data_wave=wave(tim)

def wave_array(grid,*args):
    # args=[(wlon1,wlat1),(wlon1,wlat2),(wlon2,wlat1),(wlon2,wlat2)]
    coordsw=[]
    for point in args:
        coordsw.append((point[0]+360,point[1]))
    lons=[]
    lats=[]
    for point in args:
        if point[0]<0:
            lons.append(point[0]+360)
        else:
            lons.append(point[0])
        lats.append(point[1])
    counterw=sort_counterclockwise(coordsw, centre = None)
    ilat1=np.where(grid.lat==min(lats))[0][0]
    ilat2=np.where(grid.lat==max(lats))[0][0]
    ilon1=np.where(grid.lon==min(lons))[0][0]
    ilon2=np.where(grid.lon==max(lons))[0][0]
    lat=grid.lat[ilat2:ilat1]
    lon=grid.lon[ilon1:ilon2]
    wave_arr=grid[ilat2:ilat1,ilon1:ilon2]
    return wave_arr,counterw,lat,lon    

median_timeseries[median_timeseries==0]=-180
wave_all=[]
corre_grid,corre_inde,wave_all,counterw,lat,lon=cross_correlate(time_frame,sum_tot,1,wave_all)
z=2
tim=time_frame[2]
grid1=pygmt.datasets.load_earth_relief(resolution='30s', region=[175,240, 50, 73])
os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
fig=pygmt.Figure()
os.chdir("/Users/sebinjohn/AON_PROJECT/Map_with_spectrogram")
pygmt.makecpt(cmap="bhw1_deadwood.cpt", series=[0,4000, 1],output="topo.cpt",reverse=True)
fig.grdimage(grid=grid1,region=[175,240, 50, 73],projection="Cyl_stere/200/60/20c", cmap="topo.cpt",frame="f",nan_transparent=True)
os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
# fig.basemap(region=[150,260, 30, 73], projection="L-159/35/33/45/22c")
fig.colorbar(cmap="psd.cpt",frame=["x+lPSD", "y+ldb"],projection="Cyl_stere/200/60/20c",position="n0/-0.1+w8c/0.5c+h") 
fig.grdimage(grid=grid,region=[175,240, 50, 73],projection="Cyl_stere/200/60/20c", cmap="expll1.cpt",frame="f",nan_transparent=True) 
fig.coast(region=[175,240, 50, 73], projection="Cyl_stere/200/60/20c",shorelines=True,frame="f",borders=["1/0.5p,black", "2/0.5p,red"],area_thresh='600')
fig.plot(x=median_timeseries[1,1:,z+1],y=median_timeseries[2,1:,z+1],color=median_timeseries[0,1:,z+1],style="c0.30c",cmap="psd.cpt",pen="black", projection="Cyl_stere/200/60/20c")
fig.text(text=str(tim),x=210,y=72,projection="Cyl_stere/200/60/20c")
fig.colorbar(projection="Cyl_stere/200/60/20c", cmap="expll1.cpt",frame=["x+lsignificant\twave\theight", "y+lm"],position="n0.57/-0.1+w8c/0.5c+h")
fig.text(text=["Gulf of Alaska", "Bering Sea","Chukchi Sea"], x=[-145, -178,-171], y=[58, 57,69.5])
fig.text(text="Bering Strait",x=-169.4,y=65.4,angle=90)
for i in range(len(regions)):
    p1=regions[i][0]
    p2=regions[i][1]
    p3=regions[i][2]
    p4=regions[i][3]
    reg=regions[i][-1]
    wista,poly,coords,sum_array,sum_tot=poly_seis(median_timeseries,mean_timeseries,p1,p2,p3,p4)
    xre=[]
    for ele in coords:
        xre.append(ele[0])
    xre.append(xre[0])
    yr=[]
    for ele in coords:
        yr.append(ele[1])
    yr.append(yr[0])

    xw=[]
    for ele in counterw:
        xw.append(ele[0])
    xw.append(xw[0])
    yw=[]
    for ele in counterw:
        yw.append(ele[1])
    yw.append(yw[0])
    fig.plot(x=xre,y=yr,pen="3p,Blue")
fig.show(method="external")
mattime=[]
for ele in time_frame:
    mattime.append(ele.matplotlib_date)


fig,axs=plt.subplots(nrows=1,ncols=1,figsize=(12,4))
axs.plot(mattime,wave_all[30,120])
axs.xaxis.set_major_locator(mdates.MonthLocator())
# axs.xaxis.set_major_formatter(mdates.DateFormatter())
axs.set_xticks([])
axs.set_xlim([min(mattime),max(mattime)])

p1=regions[0][0]
p2=regions[0][1]
p3=regions[0][2]
p4=regions[0][3]
reg=regions[0][-1]
wista,poly,coords,sum_array,sum_tot=poly_seis(median_timeseries,mean_timeseries,p1,p2,p3,p4)
fig,axs=plt.subplots(nrows=1,ncols=1,figsize=(12,4))
axs.plot(mattime,sum_tot[0,:])
axs.xaxis.set_major_locator(mdates.MonthLocator())
# axs.xaxis.set_major_formatter(mdates.DateFormatter())
axs.set_xticks([])
axs.set_xlim([min(mattime),max(mattime)])
