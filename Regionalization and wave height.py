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


mlat1=60
mlat2=50
mlon1=-132
mlon2=-145

def map_plo(st,et):
    mean_timeseries=np.zeros((len(stationo)+1)*3).reshape(3,len(stationo)+1)
    time_frame=[]
    for i in range (int((et-st)/(3600*windw))):
        time_frame.append(st)
        st=st+(3600*windw)
    for i in range(len(time_frame)):
        tim=time_frame[i]
        print(tim)
        mea_appen=np.zeros(3).reshape(3,1)
        integ_appen=np.zeros(3).reshape(3,1)
        mea_appen,integ_appen=mean_integ(tim,mea_appen,integ_appen)
        mean_timeseries=np.dstack((mean_timeseries,mea_appen))   
        mean_timeseries = mean_timeseries.astype('float64')
        median_timeseries=np.copy(mean_timeseries)
        #median_timeseries[median_timeseries==0]=np.nan
        np.count_nonzero(np.isnan(median_timeseries))
    for j in range(1,np.shape(mean_timeseries)[1]):
        median_timeseries[0,j,1:]=medfilt(median_timeseries[0,j,1:],7)
        print("interpolating station "+stationo[j-1][6:])
        interp_inde=np.array([])
        interp_x=median_timeseries[0,j,1:]
        datagap=np.where(interp_x==0)[0]
        data_x=np.where(interp_x!=0)[0]
        for i in range(len(data_x)-1):
            if data_x[i+1]-data_x[i]<12 and data_x[i+1]-data_x[i]>1:
                interp_inde=np.append(interp_inde,np.array( [i for i in range(int(data_x[i])+1,int(data_x[i+1]))]))
            else:
                continue
        if len(interp_inde)>1:
            interp=np.interp(interp_inde, data_x.reshape(np.shape(data_x)[0]), interp_x[data_x].reshape(np.shape(data_x)[0]))
            interp_inde=(interp_inde+1).astype("int32")
            median_timeseries[0,j,interp_inde]=interp
        else:
            continue
    # median_timeseries=np.nan_to_num(median_timeseries)
    #median_timeseries[median_timeseries==0]=np.nan  
    return mean_timeseries,median_timeseries,time_frame



        
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
    return grid,data_wave
    

def mean_integ(tim,mea_appen,integ_appen):
  for i in range(len(long)):
      os.chdir(join("/Users/sebinjohn/AON_PROJECT/Data",env[i].split("_")[-1].split(".")[-2]))
      sta=env[i].split("_")[-1].split(".")[-2]
      with open((str(sta)+".pkl"),"rb") as f:
         sta, starttimeta, endtimeta,starttimeak,endtimeak,sta,cha,loc = pickle.load(f) 
      j=int((tim-starttimeta)/(3600))
      with open((glob.glob(join(datapath+stationo[i]))[0]), 'rb') as g:
          final=np.load(g)
      res=final[[29:34,j]
      out_mean=mean(res)
      out_integ=integ(res) 
      mea_appen=np.append(mea_appen,np.array([out_mean,long[i],lat[i]]).reshape(3,1),axis=1)
      integ_appen=np.append(integ_appen,np.array([out_mean,long[i],lat[i]]).reshape(3,1),axis=1)   
  return mea_appen,integ_appen
     
def mean(res):
    su=0
    for i in range(len(res)):
        su+=res[i]
    mea=su/len(res)
    return mea

def integ(res):
    intg=trapz(res,freq[[29:34])
    return intg


def regionalization_seis(mlon1,mlon2,mlat1,mlat2,median_timeseries):
    lat_inde=[]
    final_inde=[]
    i=-1
    for ele in median_timeseries[2,:,:]:
        i+=1
        try:
            x=ele[ele>0][0]
        except:
            continue
        #print(x)
        if x<mlat1 and x>mlat2:
            lat_inde.append(i)
        else:
            pass
    for elem in lat_inde:
        lon=median_timeseries[1,elem,:]
        i=elem
        try:
            x=lon[lon<0][0]
        except:
            continue
            print(x)
        if x<mlon1 and x>mlon2:
            final_inde.append(i)
        else:
            pass
        
    sum_array=median_timeseries[0,final_inde[0],1:]
    for ele in final_inde[1:]:
        sum_array=sum_array+median_timeseries[0,ele,1:]
    return sum_array
        

def wave_val(grid,mlat1,mlat2,mlon1,mlon2):
    m_lon1=180+(180+mlon2)
    m_lon2=180+(180+mlon1)
    ilat1=np.where(grid.lat==mlat1)[0][0]
    ilat2=np.where(grid.lat==mlat2)[0][0]
    ilon1=np.where(grid.lon==m_lon1)[0][0]
    ilon2=np.where(grid.lon==m_lon2)[0][0]
    lat=grid.lat[ilat1:ilat2]
    lon=grid.lon[ilon1:ilon2]
    M1, M2 = np.meshgrid(lat, lon)
    a, b = M1.shape
    ngrid = a*b
    m1 = np.reshape(M1,(1,ngrid))
    m2 = np.reshape(M2,(1,ngrid))
    sum_wave=np.array([0.0])
    for pp in range(ngrid):
        mtry = np.array([m1[0,pp], m2[0,pp]])
        if not np.isnan(grid.loc[mtry[0],mtry[1]].data):
            sum_wave+=grid.loc[mtry[0],mtry[1]].data
    return sum_wave

  


mean_timeseries,median_timeseries,time_frame=map_plo(st,et)


sum_array=regionalization_seis(mlon1,mlon2,mlat1,mlat2,median_timeseries)
sum_wave=np.array([])  
sum_array=regionalization_seis(lon1,lon2,lat1,lat2,median_timeseries)
for i in range(len(time_frame)):
    tim=time_frame[i]
    print("searching for wave height for"+str(tim))
    grid,data_wave=wave(tim)
    sum_w=wave_val(grid,mlat1,mlat2,mlon1,mlon2)
    sum_wave=np.append(sum_wave,sum_w)

plt.plot(sum_wave)
plt.plot(sum_array)

