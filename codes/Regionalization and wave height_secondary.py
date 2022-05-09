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

mlat1=64
mlat2=61
mlon1=-142
mlon2=-151


wlat1=62
wlat2=60
wlon1=-142
wlon2=-151

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


def poly_wave(grid,*args):
    #args=[(wlon1,wlat1),(wlon1,wlat2),(wlon2,wlat1),(wlon2,wlat2)]
    coordsw=[]
    for point in args:
        coordsw.append((point[0]+360,point[1]))
    lons=[]
    lats=[]
    for point in args:
        lons.append(point[0]+360)
        lats.append(point[1])
    counterw=sort_counterclockwise(coordsw, centre = None)
    polyw = Polygon(counterw)
    ilat1=np.where(grid.lat==min(lats))[0][0]
    ilat2=np.where(grid.lat==max(lats))[0][0]
    ilon1=np.where(grid.lon==min(lons))[0][0]
    ilon2=np.where(grid.lon==max(lons))[0][0]
    lat=grid.lat[ilat2:ilat1]
    lon=grid.lon[ilon1:ilon2]
    pointsw=[]
    for i in range (len(lon)):
        for j in range (len(lat)):
            pointsw.append((lon[i],lat[j]))
    sum_wave=np.array([0.0])
    nopoints=0
    for ele in pointsw:
        p1 = Point(ele)
        if p1.within(polyw):
            if not np.isnan(grid.loc[ele[1],ele[0]].data):
                nopoints+=1
                sum_wave+=grid.loc[ele[1],ele[0]].data
    return sum_wave,counterw,nopoints
    
    


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
with open(("median and mean"+str(time_frame[0])+"-"+str(time_frame[-1])+"secondary.pkl"), 'rb') as f:  # Python 3: open(..., 'wb')
          mean_timeseries,median_timeseries=pickle.load(f)



wista,poly,coords,sum_array,sum_tot=poly_seis(median_timeseries,mean_timeseries,(-140,62),(-140.5,59.5),(-131.5,55),(-129,56))  

#mean_timeseries,median_timeseries,time_frame=map_plo(st,et)

#sum_array,sum_mean=regionalization_seis(mlon1,mlon2,mlat1,mlat2,median_timeseries)
sum_wave=np.array([])  
#sum_array=regionalization_seis(lon1,lon2,lat1,lat2,median_timeseries)
for i in range(len(time_frame)):
    tim=time_frame[i]
    print("searching for wave height for"+str(tim))
    grid,data_wave=wave(tim)
    #sum_w=wave_val(grid,wlat1,wlat2,wlon1,wlon2)
    sum_w,counterw,nopoints=poly_wave(grid,(-141.5,60),(-131,54.5),(-135,54),(-144,58.5))
    sum_wave=np.append(sum_wave,sum_w)
    sum_wave_n=sum_wave/nopoints

time_framem=[]
for i in range(len(time_frame)):
    time_framem.append(time_frame[i].matplotlib_date)
    
date_format=mdates.DateFormatter('%d,%b,%y')

mean_timeseries_plo=mean_timeseries.copy()
mean_timeseries_plo[mean_timeseries_plo==0]=np.nan







fig,axs=plt.subplots(4,1,figsize=(25,20))
axs[0].plot(time_framem,sum_wave_n)
axs[0].xaxis.set_major_formatter(date_format)
axs[0].set_ylabel("significant wave_height")
axs[0].set_title("mean wave Height over the region ")
axs[1].plot(time_framem,sum_tot[0,:])
axs[1].xaxis.set_major_formatter(date_format)
axs[1].set_ylabel("PSD")
axs[1].set_title("mean PSD over the stations in the region_7 hour median Filtered ")
axs[2].plot(time_framem,sum_tot[1,:])
axs[2].xaxis.set_major_formatter(date_format)
axs[2].set_title("mean raw-PSD over the stations in the region ")
axs[2].set_ylabel("PSD")
for i in range(len(wista)):
    if wista[i]==1:
        axs[3].plot(time_framem,mean_timeseries_plo[0,i,1:])
axs[3].xaxis.set_major_formatter(date_format)
axs[3].set_title("raw-PSD of stations in region-each line is a station ")
axs[3].set_ylabel("PSD")

for i in range(4):
    axs[i].set_xlim([18205,18214])
   # axs[i].set_xlim([min(time_framem),max(time_framem)])
for i in range(1,4):
    axs[i].set_ylim([-145,-100])


x=[]
for ele in coords:
    x.append(ele[0])
x.append(x[0])
y=[]
for ele in coords:
    y.append(ele[1])
y.append(y[0])

xw=[]
for ele in counterw:
    xw.append(ele[0])
xw.append(xw[0])
yw=[]
for ele in counterw:
    yw.append(ele[1])
yw.append(yw[0])


fig=pygmt.Figure()
#fig.basemap(region=[150,260, 30, 73],frame=True, projection="L-159/35/33/45/22c")
#fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="expll1.cpt",frame="a",nan_transparent=True)
fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/22c",shorelines=True,frame="a")
#fig.colorbar(cmap="exp.cpt",frame='af+l"PSD (db)"',projection="L-159/35/33/45/22c")  
#fig.plot(x=median_timeseries[1,1:,z+1],y=median_timeseries[2,1:,z+1],color=median_timeseries[0,1:,z+1],style="c0.30c",cmap="exp.cpt",pen="black", projection="L-159/35/33/45/22c")
#fig.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/22c")
fig.plot(x=x,y=y,pen="1p,red")
fig.plot(x=xw,y=yw,pen="1p,blue")

fig.show(method="external")

# def regionalization_seis(mlon1,mlon2,mlat1,mlat2,median_timeseries):
#     lat_inde=[]
#     final_inde=[]
#     i=-1
#     for ele in median_timeseries[2,:,:]:
#         i+=1
#         try:
#             x=ele[ele>0][0]
#         except:
#             continue
#         #print(x)
#         if x<mlat1 and x>mlat2:
#             lat_inde.append(i)
#         else:
#             pass
#     for elem in lat_inde:
#         lon=median_timeseries[1,elem,:]
#         i=elem
#         try:
#             x=lon[lon<0][0]
#         except:
#             continue
#             print(x)
#         if x<mlon1 and x>mlon2:
#             final_inde.append(i)
#         else:
#             pass
        
#     sum_array=median_timeseries[0,final_inde[0],1:]
#     sum_mean=mean_timeseries[0,final_inde[0],1:]
#     for ele in final_inde[1:]:
#         sum_array=sum_array+median_timeseries[0,ele,1:]
#         sum_mean=sum_mean+mean_timeseries[0,ele,1:]
#     return sum_array,sum_mean
        

# def wave_val(grid,mlat1,mlat2,mlon1,mlon2):
#     m_lon1=180+(180+mlon2)
#     m_lon2=180+(180+mlon1)
#     ilat1=np.where(grid.lat==mlat1)[0][0]
#     ilat2=np.where(grid.lat==mlat2)[0][0]
#     ilon1=np.where(grid.lon==m_lon1)[0][0]
#     ilon2=np.where(grid.lon==m_lon2)[0][0]
#     lat=grid.lat[ilat1:ilat2]
#     lon=grid.lon[ilon1:ilon2]
#     M1, M2 = np.meshgrid(lat, lon)
#     a, b = M1.shape
#     ngrid = a*b
#     m1 = np.reshape(M1,(1,ngrid))
#     m2 = np.reshape(M2,(1,ngrid))
#     sum_wave=np.array([0.0])
#     for pp in range(ngrid):
#         mtry = np.array([m1[0,pp], m2[0,pp]])
#         if not np.isnan(grid.loc[mtry[0],mtry[1]].data):
#             sum_wave+=grid.loc[mtry[0],mtry[1]].data
#     return sum_wave

  


# sum_array,sum_mean=regionalization_seis(mlon1,mlon2,mlat1,mlat2,median_timeseries)
# sum_wave=np.array([])  
# #sum_array=regionalization_seis(lon1,lon2,lat1,lat2,median_timeseries)
# for i in range(len(time_frame)):
#     tim=time_frame[i]
#     print("searching for wave height for"+str(tim))
#     grid,data_wave=wave(tim)
#     sum_w=wave_val(grid,wlat1,wlat2,wlon1,wlon2)
#     sum_wave=np.append(sum_wave,sum_w)


plt.plot(time_framem,median_timeseries[0,15,1:])



# os.chdir("/Users/sebinjohn/AON_PROJECT/Data/median_Power_time_series")
# with open(("median and mean"+str(time_frame[0])+"-"+str(time_frame[-1])+"secondary.pkl"), 'wb') as f:  # Python 3: open(..., 'wb')
#           pickle.dump([mean_timeseries,median_timeseries], f)