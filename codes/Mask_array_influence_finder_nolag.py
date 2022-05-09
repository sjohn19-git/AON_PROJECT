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
      res=final[34:42,j]
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
    intg=trapz(res,freq[34:42])
    return intg

def sort_counterclockwise(points, centre = None):
  if centre:
    centre_x, centre_y = centre
  else:
    centre_x, centre_y = sum([x for x,_ in points])/len(points), sum([y for _,y in points])/len(points)
  angles = [math.atan2(y - centre_y, x - centre_x) for x,y in points]
  counterclockwise_indices = sorted(range(len(points)), key=lambda i: angles[i])
  counterclockwise_points = [points[i] for i in counterclockwise_indices]
  return counterclockwise_points


def poly_seis(median_timeseries,*args):
    coords=[]
    for point in args:
        coords.append(point)
    counter=sort_counterclockwise(coords, centre = None)
    poly = Polygon(counter)
    wista=[]
    for i in range(np.shape(median_timeseries)[1]):
        p1 = Point(median_timeseries[1,i,5],median_timeseries[2,i,5] )
        if p1.within(poly):
            wista.append(1)
        else:
            wista.append(0)
    c=0
    for i in range(len(wista)):
        if wista[i]==1:
            median_timeseries[median_timeseries==0]=np.mean(median_timeseries[0,i,1:])
            mean_timeseries[mean_timeseries==0]=np.mean(mean_timeseries[0,i,1:])
            if c==0:
                sum_array=median_timeseries[0,i,1:]
                sum_mean=mean_timeseries[0,i,1:]
                c=1
            else:
                sum_array=sum_array+median_timeseries[0,i,1:]
                sum_mean=sum_mean+mean_timeseries[0,i,1:]
    return wista,poly,counter,sum_array,sum_mean
        
    
    


def poly_seis(median_timeseries,mean_timeseries,*args):
    coords=[]
    for point in args:
        coords.append(point)
    counter=sort_counterclockwise(coords, centre = None)
    poly = Polygon(counter)
    wista=[]
    for i in range(np.shape(median_timeseries)[1]):
        p1 = Point(median_timeseries[1,i,5],median_timeseries[2,i,5] )
        if p1.within(poly):
            wista.append(1)
        else:
            wista.append(0)
    c=wista.count(1)
    j=0
    sum_array=np.empty([c,2,len(median_timeseries[0,0,1:])])
    for i in range(len(wista)):
        if wista[i]==1:
            sum_array[j,0,:]=median_timeseries[0,i,1:]
            sum_array[j,1,:]=mean_timeseries[0,i,1:]
            j+=1
    sum_tot=np.empty([2,len(median_timeseries[0,0,1:])])
    for i in range(len(sum_array[0,0,:])):
        sum_tot[0,i]=np.sum(sum_array[:,0,i])/(c-np.count_nonzero(sum_array[:,0,i]==0))
        sum_tot[1,i]=np.sum(sum_array[:,1,i])/(c-np.count_nonzero(sum_array[:,1,i]==0))
    return wista,poly,counter,sum_array,sum_tot


os.chdir("/Users/sebinjohn/AON_PROJECT/Data/median_Power_time_series")
with open(("median and mean"+str(time_frame[0])+"-"+str(time_frame[-1])+"Secondary.pkl"), 'rb') as f:  # Python 3: open(..., 'wb')
          mean_timeseries,median_timeseries=pickle.load(f)


wista,poly,coords,sum_array,sum_tot=poly_seis(median_timeseries,mean_timeseries,(-151.95,70.3),(-151.6,66.5),(-141.6,66.5),(-141.88,69.4)) 

def wave_array(grid,*args):
    #args=[(wlon1,wlat1),(wlon1,wlat2),(wlon2,wlat1),(wlon2,wlat2)]
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

def cross_correlate(time_frame,sum_tot,itera,*argv):
    c=0
    if itera==1:
        for i in range(len(time_frame)):
            tim=time_frame[i]
            print("searching for wave height for"+str(tim))
            grid,data_wave=wave(tim)
            #sum_w=wave_val(grid,wlat1,wlat2,wlon1,wlon2)
            wave_arr,counterw,lat,lon=wave_array(grid,(150,30),(260,73),(260,30),(150,73))
            if c==0:
               a,b=np.shape(wave_arr.data)
               wave_all=np.zeros([a,b,len(time_frame)])
               c=1
            wave_all[:,:,i]=wave_arr.data
    else:
        for i in range(1):
            tim=time_frame[i]
            print("searching for wave height for"+str(tim))
            grid,data_wave=wave(tim)
            #sum_w=wave_val(grid,wlat1,wlat2,wlon1,wlon2)
            wave_arr,counterw,lat,lon=wave_array(grid,(150,30),(260,73),(260,30),(150,73))
        wave_all=argv[0]
    corre_grid=np.zeros([np.shape(wave_all)[0],np.shape(wave_all)[1]])
    corre_inde=np.zeros(np.shape(corre_grid))
    for i in range(np.shape(wave_all)[0]):
        for j in range(np.shape(wave_all)[1]):
            k=cr.correlate(wave_all[i,j,:],sum_tot[0,:],720)
            corre_grid[i,j]=cr.xcorr_max(k)[1]
            corre_inde[i,j]=cr.xcorr_max(k)[0]
    return corre_grid,corre_inde,wave_all,counterw,lat,lon
            
wave_all=[]

corre_grid,corre_inde,wave_all,counterw,lat,lon=cross_correlate(time_frame,sum_tot,1,wave_all)

def masked_corre_grid(corre_grid,corre_inde):
    bool_corre=corre_grid>0.6
    bool_inde=np.logical_and(corre_inde < 24, corre_inde>-1)
    comb_bool=np.logical_and(bool_corre,bool_inde)
    comb_corre=corre_grid.copy()
    comb_corre[~comb_bool] = np.nan
    return comb_corre

comb_corre=masked_corre_grid(corre_grid,corre_inde)



grid_corre = xr.DataArray(
    corre_grid, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
    )

corre_i=corre_inde.copy()
corre_i[corre_i<-720]=np.nan


grid_inde = xr.DataArray(
    corre_i, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
    )

grid_corre_mask=xr.DataArray(comb_corre, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
    )



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


os.chdir("/Users/sebinjohn/AON_PROJECT/cpts")


pygmt.makecpt(cmap="hot",output="hot1.cpt", series=[0,1,0.1],reverse=True)
fig=pygmt.Figure()


#fig.plot(x=xw,y=yw,pen="1p,blue")
#fig.plot(x=x,y=y,region=[150,260,30, 73], projection="L-159/35/33/45/22c",style="i0.5c", color="black",land=False)
#fig.contour(x=x,y=y,z=z,region=[150,260,30, 73], projection="L-159/35/33/45/22c",pen="0.5p",levels=0.2,annotation=0.2)
fig.grdimage(grid=grid_corre,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="hot1.cpt",frame="a",nan_transparent=True)
fig.coast(region=[150,260,30, 73], projection="L-159/35/33/45/22c",shorelines=True,frame="a",land="white")
fig.colorbar(cmap="hot1.cpt",frame='af+l"correlation_coeifficent"',projection="L-159/35/33/45/22c") 
fig.plot(x=xre,y=yr,pen="1p,red")
fig.show(method="external")




pygmt.makecpt(cmap="hot",output="hot1.cpt", series=[0.5,1,0.01],reverse=True)
fig=pygmt.Figure()


#fig.plot(x=xw,y=yw,pen="1p,blue")
#fig.plot(x=x,y=y,region=[150,260,30, 73], projection="L-159/35/33/45/22c",style="i0.5c", color="black",land=False)
#fig.contour(x=x,y=y,z=z,region=[150,260,30, 73], projection="L-159/35/33/45/22c",pen="0.5p",levels=0.2,annotation=0.2)
fig.grdimage(grid=grid_corre_mask,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="hot1.cpt",frame="a",nan_transparent=True)
fig.coast(region=[150,260,30, 73], projection="L-159/35/33/45/22c",shorelines=False,frame="a",land="200/200/200",borders=["1/0.5p,black", "2/0.5p,red"])
fig.colorbar(cmap="hot1.cpt",frame='af+l"correlation_coeifficent"',projection="L-159/35/33/45/22c") 
fig.plot(x=xre,y=yr,pen="1p,red")
fig.show(method="external")


#pygmt.makecpt(cmap="polar",output="polar1.cpt", series=[-100,100,1],reverse=False)
fig=pygmt.Figure()

#fig.plot(x=x,y=y,pen="1p,red")
#fig.plot(x=xw,y=yw,pen="1p,blue")
#fig.plot(x=x,y=y,region=[150,260,30, 73], projection="L-159/35/33/45/22c",style="i0.5c", color="black",land=False)
#fig.contour(x=x,y=y,z=z,region=[150,260,30, 73], projection="L-159/35/33/45/22c",pen="0.5p",levels=0.2,annotation=0.2)
fig.grdimage(grid=grid_inde,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="polar1.cpt",frame="a",nan_transparent=True)
fig.coast(region=[150,260,30, 73], projection="L-159/35/33/45/22c",shorelines=True,frame="a",land="white")
fig.colorbar(cmap="polar1.cpt",frame='af+l"correlation_index"',projection="L-159/35/33/45/22c") 
fig.plot(x=xre,y=yr,pen="1p,red")
 
fig.show(method="external")

