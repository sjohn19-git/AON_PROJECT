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

os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making/ice_video")


with open("metadta.pkl","rb") as f:
   long,lat,stationo ,env = pickle.load(f)

freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)
len(freq)


pygmt.makecpt(cmap="hot",output="exp.cpt", series=[-180,-90,0.005],reverse=True)
pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[0,14],reverse=True)
#grid = pygmt.datasets.load_earth_relief(resolution="10m", region=[-172, -135, 52, 73])
datapath="/Users/sebinjohn/AON_PROJECT/Data/*/"
st=UTCDateTime(2019,6,1)
et=UTCDateTime(2020,8,1)
windw=24
tim_len=int((et-st)/(3600*windw))
smth_intg=np.zeros((len(stationo),tim_len))


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
    median_timeseries[median_timeseries==0]=-200
    for i in range(len(time_frame)):
        tim=time_frame[i]
        print("searching for wave height for"+str(tim))
        grid,data_wave=wave(tim)
        print("plotting "+str(tim))
        ma_pic_generator(tim,grid,median_timeseries,i)
    return mean_timeseries,median_timeseries


        
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
    
    

def ma_pic_generator(tim,grid,median_timeseries,z):
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making/ice_video")
    fig=pygmt.Figure()
    #fig.basemap(region=[150,260, 30, 73],frame=True, projection="L-159/35/33/45/22c")
    fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="S200/65/20/20c", cmap="expll1.cpt",frame="a",nan_transparent=True)
    fig.coast(region=[150,260, 30, 73], projection="S200/65/20/20c",shorelines=True,frame="a")
    fig.colorbar(cmap="exp.cpt",frame='af+l"PSD (db)"',projection="S200/65/20/20c")  
    fig.plot(x=median_timeseries[1,1:,z+1],y=median_timeseries[2,1:,z+1],color=median_timeseries[0,1:,z+1],style="c0.45c",cmap="exp.cpt",pen="black", projection="S200/65/20/20c")
    fig.text(text=str(tim),x=200,y=80,projection="S200/65/20/20c")
    fig.savefig(str(tim)+"wave.jpg")



def mean_integ(tim,mea_appen,integ_appen):
  for i in range(len(long)):
      os.chdir(join("/Users/sebinjohn/AON_PROJECT/Data",env[i].split("_")[-1].split(".")[-2]))
      sta=env[i].split("_")[-1].split(".")[-2]
      with open((str(sta)+".pkl"),"rb") as f:
         sta, starttimeta, endtimeta,starttimeak,endtimeak,sta,cha,loc = pickle.load(f) 
      j=int((tim-starttimeta)/(3600))
      with open((glob.glob(join(datapath+stationo[i]))[0]), 'rb') as g:
          final=np.load(g)
      res=final[60:69,j]
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
    intg=trapz(res,freq[60:69])
    return intg
    




def convert_pictures_to_video(pathIn, pathOut, fps, time):
    ''' this function converts images to video'''
    frame_array=[]
    files=[f for f in glob.glob(pathIn) if isfile(join(pathIn,f))]
    files.sort()
    files
    for i in range(len(files)):
        print(i)
        '''reading images'''
        pil_image =Image.open(files[i]).convert('RGB') 
        height, width = pil_image.size
        open_cv_image = np.array(pil_image) 
        # Convert RGB to BGR 
        img = open_cv_image[:, :, ::-1].copy() 
        size=(height,width)
        for k in range (time):
            frame_array.append(img)
    out=cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'VP80'), fps,size)
    for i in range(len(frame_array)):
        out.write(frame_array[i])
    out.release()

# Example:
directory="//Users/sebinjohn/AON_PROJECT/Data/Video Making/ice_video"
pathIn=directory+'/*.jpg'
pathOut=directory+"/"+str(st)+"-"+str(et)+"2.avi"
fps=8
time=1 # the duration of each picture in the video




mean_timeseries,median_timeseries=map_plo(st,et)



convert_pictures_to_video(pathIn, pathOut, fps, time)



len(env)


#--------------------------------------------------------------------------------
pics=glob.glob(join(directory+"/*wave.jpg"))
len(pics)
for ele in glob.glob(join(directory+"/*plane.jpg")):
    pics.append(ele)
for ele in pics:
    os.remove(ele)
#--------------------------------------------------------------------------------
fig = pygmt.Figure()
fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/22c",shorelines=True)
fig.plot(
    region=[150,260, 30, 73],
    projection="L-159/35/33/45/22c",
    frame="ag",
    x=[-149.5432],
    y=[65.8251],
    style="v0.6c+e",
    direction=[[56], [-0.875]],
    pen="2p",
    color="red3",
)
fig.show(method="external")




vector_1 = [ -0.299093, -0.451724]
vector_2 = [-1, 0]

unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
dot_product = np.dot(unit_vector_1, unit_vector_2)
angle = np.arccos(dot_product)

print(angle*180/3.14)





