
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

pygmt.makecpt(cmap="hot",output="psd.cpt", series=[-150,-90,0.005],reverse=True)
pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[0,14],reverse=True)


#grid = pygmt.datasets.load_earth_relief(resolution="10m", region=[-172, -135, 52, 73])
datapath="/Users/sebinjohn/AON_PROJECT/Data/*/"
st=UTCDateTime(2019,11,1)
et=UTCDateTime(2019,12,1)
windw=1
tim_len=int((et-st)/(3600*windw))
smth_intg=np.zeros((len(stationo),tim_len))
os.chdir("/Users/sebinjohn/AON_PROJECT/Data/median_Power_time_series")
with open(("median and mean"+str(time_frame[0])+"-"+str(time_frame[-1])+"Secondary.pkl"), 'rb') as f:  # Python 3: open(..., 'wb')
          mean_timeseries,median_timeseries=pickle.load(f)
median_timeseries_dif=median_timeseries.copy()
median_timeseries[median_timeseries==0]=-180

os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")

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
    median_timeseries[median_timeseries==0]=-180
    for i in range(len(time_frame)):
        tim=time_frame[i]
        print("searching for wave height for"+str(tim))
        grid,data_wave=wave(tim)
        print("plotting "+str(tim))
        ma_pic_generator(tim,grid,median_timeseries,median_timeseries_dif,i)
        #fig_merge(tim)
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
    
    

def ma_pic_generator(tim,grid,median_timeseries,median_timeseries_dif,z):
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
    fig=pygmt.Figure()
    #fig.basemap(region=[150,260, 30, 73], projection="L-159/35/33/45/22c")
    fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/20c", cmap="expll1.cpt",frame="f",nan_transparent=True)
    fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/20c",shorelines=True,frame="f",borders=["1/0.5p,black", "2/0.5p,red"],area_thresh='600')
    fig.colorbar(cmap="psd.cpt",frame=["x+lPSD", "y+ldb"],projection="L-159/35/33/45/20c",position="n0/-0.1+w8c/0.5c+h")  
    fig.plot(x=median_timeseries[1,1:,z+1],y=median_timeseries[2,1:,z+1],color=median_timeseries[0,1:,z+1],style="c0.30c",cmap="psd.cpt",pen="black", projection="L-159/35/33/45/20c")
    fig.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/20c")
    angle,vector_1=plane_fitting(tim,grid,median_timeseries_dif,i)
    #fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/20c", cmap="expll1.cpt",frame="f",nan_transparent=True)
    #fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/20c",shorelines=True)
    fig.colorbar(projection="L-159/35/33/45/20c", cmap="expll1.cpt",frame=["x+lsignificant\twave\theight", "y+lm"],position="n0.57/-0.1+w8c/0.5c+h")
    fig.plot(
        region=[150,260, 30, 73],
        projection="L-159/35/33/45/20c",
        x=[-149.5432],
        y=[65.8251],
        frame="f",
        style="v0.6c+e",
        direction=[[angle], [-3*np.linalg.norm(vector_1)]],
        pen="2p",
        color="black",
    )
    fig.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/20c")
    #fig.show(method="external")
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making/one_panel")
    fig.savefig(str(tim)+"wave.jpg")



def plane_fitting(tim,grid,median_timeseries_dif,z):
        for i in range(median_timeseries_dif.shape[1]):
            mean_p=np.mean(median_timeseries_dif[0,i,:])
            median_timeseries_dif[0,i,:]=median_timeseries_dif[0,i,:]-mean_p
        xs=median_timeseries_dif[1,1:,z+1].copy()
        ys=median_timeseries_dif[2,1:,z+1].copy()
        zs=median_timeseries_dif[0,1:,z+1].copy()
        dat_gap=np.where(np.logical_or(zs==0,zs>40))
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
        print("errors: \n", errors)
        print("residual:", residual)
        vector_1 = [fit[0,0],fit[1,0]]
        vector_2 = [-1, 0]
        unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
        unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        angle = np.arccos(dot_product)*(180/3.14159)
        print("angle is =",angle)
        plt.figure(figsize=(12,12))
        ax = plt.subplot(111, projection='3d')
        ax.scatter(xs, ys, zs, color='b')
        # plot plane
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                          np.arange(ylim[0], ylim[1]))
        Z = np.zeros(X.shape)
        for r in range(X.shape[0]):
            for c in range(X.shape[1]):
                Z[r,c] = fit[0] * X[r,c] + fit[1] * Y[r,c] + fit[2]
        # ax.plot_wireframe(X,Y,Z, color='k')
        # ax.set_xlabel('lat')
        # ax.set_ylabel('lon')
        # ax.set_zlabel('power')
        # ax.set_title(str(tim))
        # #plt.savefig(str(tim)+"plane.jpg",dpi=700)
        # plt.close()
        return angle,vector_1



def fig_merge(tim):
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
    image1 = Image.open(str(tim)+"wave.jpg")
    image1_size = image1.size
    new_image = Image.new('RGB',(7423,3575), (250,250,250))
    new_image.paste(image1,(0,0))
    new_image.save(str(tim)+".jpg","JPEG")



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
    pics=glob.glob(join(directory+"/*wave.jpg"))
    frame_array=[]
    files=[f for f in glob.glob(pathIn) if isfile(join(pathIn,f))]
    files.sort()
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
directory="//Users/sebinjohn/AON_PROJECT/Data/Video Making/one_panel"
pathIn=directory+'/*.jpg'
pathOut=directory+"/"+str(st)+"-"+str(et)+".avi"
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
#--------------------------------------------------------------------------------
