import os 
os.chdir("/Users/sebinjohn/AON_PROJECT/Data/codes")
import pygmt
#import dill
import pandas as pd
from obspy import UTCDateTime
import cv2
from os.path import isfile, join 
import glob
from scipy.integrate import trapz 
import numpy as np
import pickle
import pygrib
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cdsapi
from medfilt import medfilt
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import matplotlib.dates as mdates
from PIL import Image
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from shapely.geometry import Point, Polygon
import math 
import obspy.signal.cross_correlation as cr


os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
with open("metadta.pkl","rb") as f:
   long,lat,stationo ,env = pickle.load(f) 

datapath="/Users/sebinjohn/AON_PROJECT/Data/*/"
st=UTCDateTime(2019,11,1)
et=UTCDateTime(2019,12,1)
windw=1
tim_len=int((et-st)/(3600*windw))
smth_intg=np.zeros((len(stationo),tim_len))


freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)




time_frame=[]
for i in range (int((et-st)/(3600*windw))):
    time_frame.append(st)
    st=st+(3600*windw)



def wave(tim):
    c = cdsapi.Client()
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/wave")
    files=os.listdir()
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
    tim_wav=str(tim)[0:14]+"00"
    grib=str(tim)[0:14]+"00"+".grib"
    grbs=pygrib.open(grib)
    grb=grbs[1]
    data_wave = grb.values
    latg, long = grb.latlons()
    lat=(latg[:,0])
    lon=(long[0,:])
    grid = xr.DataArray(
        data_wave, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
        )
    return grid,data_wave,lat,lon
 
grid,data_wave=wave(tim)

def wave_grid_maker(data_wave,i):
    wave_array_comb[:,:,i]=data_wave
    return wave_array_comb
    
    
wave_array_comb=np.zeros((361,720,len(time_frame)))     
for i in range(len(time_frame)):
    print(i)
    tim=time_frame[i]
    grid,data_wave,lat,lon=wave(tim)
    wave_array_comb=wave_grid_maker(data_wave.filled(),i)
    wave_array_comb[wave_array_comb==9999.0]=np.nan
data_wave_mean=np.sum(wave_array_comb,axis=2)/np.shape(time_frame)[0]
grid_mean = xr.DataArray(
    data_wave_mean, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
    ) 

fig=pygmt.Figure()   
fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="expll1.cpt",frame="a",nan_transparent=True)
fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/22c",shorelines=True,frame="a")
fig.colorbar(cmap="expll1.cpt",frame='af+l"mean Significant wave height (m)"',projection="L-159/35/33/45/22c")  
#fig.plot(x=median_timeseries[1,1:,z+1],y=median_timeseries[2,1:,z+1],color=median_timeseries[0,1:,z+1],style="c0.30c",cmap="exp.cpt",pen="black", projection="L-159/35/33/45/22c")
#fig.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/22c")
fig.show(method="external")

data_wave.masked_values()
wave_array_comb
    
    