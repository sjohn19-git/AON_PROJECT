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

freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)
len(freq)

angle=np.array([])
wave_angle=np.array([])
datapath="/Users/sebinjohn/AON_PROJECT/Data/*/"
st=UTCDateTime(2019,11,1)
et=UTCDateTime(2019,12,1)
windw=1
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
        #print("interpolating station "+stationo[j-1][6:])
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
    median_timeseries[median_timeseries==0]=-180
    for i in range(len(time_frame)):
        tim=time_frame[i]
        print("searching for wave height for"+str(tim))
        grid,gridf,data_wave=wave(tim)
        print("plotting "+str(tim))
        angle,wave_angle=plane_fitting(tim,gridf,grid,median_timeseries,i,angle,wave_angle)
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
    gridf=xr.DataArray(
        data_wave, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
        )
    grid = xr.DataArray(data_wave[34:125,300:521], dims=["lat", "lon"], coords={"lat": lat[34:125], "lon": lon[300:521]})
    return grid,gridf,data_wave
    
   

def plane_fitting(tim,gridf,grid,median_timeseries,z,angle,wave_angle):
    try:
        xs=median_timeseries[1,1:,z+1].copy()
        ys=median_timeseries[2,1:,z+1].copy()
        zs=median_timeseries[0,1:,z+1].copy()
        dat_gap=np.where(np.logical_or(zs>=-90,zs<=-160))
        xs = np.delete(xs,dat_gap)
        ys=np.delete(ys,dat_gap )
        zs=np.delete(zs,dat_gap)
        # plot raw data
        # do fit
        tmp_A = []
        tmp_b = []
        for i in range(len(xs)):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_b.append(zs[i])
        b = np.matrix(tmp_b).T
        A = np.matrix(tmp_A)
        # Manual solution
        fit = (A.T * A).I * A.T * b
        errors = b - A * fit
        residual = np.linalg.norm(errors)
        print("solution: %f x + %f y + %f = z" % (fit[0], fit[1], fit[2]))
       # print("errors: \n", errors)
        print("residual:", residual)
        vector_1 = [fit[0,0],fit[1,0]]
        vector_2 = [-1, 0]
        unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
        unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        angl = np.arccos(dot_product)*(180/3.14159)
        angle=np.append(angle,angl)
        print("angle is =",angl)
        lonw=grid.where(grid==grid.max(), drop=True).lon.data
        latw=grid.where(grid==grid.max(), drop=True).lat.data
        m=(65.821-latw)/((360-149.5432)-lonw)
        wave_angl=180-(np.arctan(-m)*(180/3.14159))
        wave_angle=np.append(wave_angle,wave_angl)
        print("wave angle is =",wave_angl)
        fig2 = pygmt.Figure()
        fig2.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/15c", cmap="expll1.cpt",frame="a",nan_transparent=True)
        fig2.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/15c",shorelines=True)
        fig2.plot(
            region=[150,260, 30, 73],
            projection="L-159/35/33/45/15c",
            x=[-149.5432],
            y=[65.8251],
            frame="a",
            style="v0.6c+e",
            direction=[[angl], [-1.5*np.linalg.norm(vector_1)]],
            pen="2p",
            color="red3",
        )
        dire2=(grid.where(grid==grid.max(), drop=True).data[0][0])/10
        if m<0:
            dire2=-1*dire2
        fig2.plot(
            region=[150,260, 30, 73],
            projection="L-159/35/33/45/15c",
            x=[-149.5432],
            y=[65.8251],
            frame="a",
            style="v0.6c+e",
            direction=[[wave_angl][0], [dire2]],
            pen="2p",
            color="blue3",
        )
        fig2.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/22c")
        fig2.savefig("/Users/sebinjohn/AON_PROJECT/Data/Angle_comp/"+str(tim)+"angle_comp.jpg")
    except:
        print("passed "+str(tim))
        pass
    return angle,wave_angle

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
    





def convert_pictures_to_video(pathIn, pathOut, fps, time):
    ''' this function converts images to video'''
    pics=glob.glob(join(directory+"/*.jpg"))
    frame_array=[]
    files=[f for f in glob.glob(pathIn) if isfile(join(pathIn,f))]
    files.sort()
    files
    for i in range (len(files)):
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
directory="//Users/sebinjohn/AON_PROJECT/Data/Angle_comp"
pathIn=directory+'/*.jpg'
pathOut=directory+"/"+str(st)+"-"+str(et)+".avi"
fps=8
time=1 



mean_timeseries,median_timeseries=map_plo(st,et)


fig2 = pygmt.Figure()
fig2.coast(region="g", projection="Cyl_stere/180/-40/20c",shorelines=True)
fig2.grdimage(grid=grid,projection="Cyl_stere/180/-40/20c",frame="a", cmap="expll1.cpt",nan_transparent=True)
fig2.savefig("/Users/sebinjohn/AON_PROJECT/Data/Angle_comp/"+str(tim)+"angle_comdefjp.jpg")


