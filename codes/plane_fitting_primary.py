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

os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making/primary_video")


with open("metadta.pkl","rb") as f:
   long,lat,stationo ,env = pickle.load(f)

freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)
len(freq)


pygmt.makecpt(cmap="hot",output="exp.cpt", series=[-165,-105,0.005],reverse=True)
pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[0,14],reverse=True)
#grid = pygmt.datasets.load_earth_relief(resolution="10m", region=[-172, -135, 52, 73])
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
        median_timeseries[0,j,1:]=medfilt(median_timeseries[0,j,1:],13)
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
        ma_pic_generator(tim,grid,median_timeseries,i)
        plane_fitting(tim,grid,median_timeseries,i)
        fig_merge(tim)
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
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making/primary_video")
    fig=pygmt.Figure()
    with fig.subplot(nrows=1, ncols=2, figsize=("45c", "60c"), frame="lrtb",clearance='2c',title=""):
        with fig.set_panel(panel=0): 
            #fig.basemap(region=[150,260, 30, 73],frame=True, projection="L-159/35/33/45/22c")
            fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="expll1.cpt",frame="a",nan_transparent=True)
            fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/22c",shorelines=True,frame="a")
            fig.colorbar(cmap="exp.cpt",frame='af+l"PSD (db)"',projection="L-159/35/33/45/22c")  
            fig.plot(x=median_timeseries[1,1:,z+1],y=median_timeseries[2,1:,z+1],color=median_timeseries[0,1:,z+1],style="c0.30c",cmap="exp.cpt",pen="black", projection="L-159/35/33/45/22c")
            fig.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/22c")
        with fig.set_panel(panel=1): 
            fig.coast(region="g", projection="Cyl_stere/180/-40/20c",shorelines=True)
            fig.grdimage(grid=grid,projection="Cyl_stere/180/-40/20c",frame="a", cmap="expll1.cpt",nan_transparent=True)
            fig.colorbar(projection="Cyl_stere/180/-40/20c", cmap="expll1.cpt",frame='af+l"Significant Wave Height (m)"')
            fig.text(text=str(tim),x=190,y=85,projection="Cyl_stere/180/-40/20c")
        fig.savefig(str(tim)+"wave.jpg")
 

def plane_fitting(tim,grid,median_timeseries,z):
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
        ax.plot_wireframe(X,Y,Z, color='k')
        ax.set_xlabel('lat')
        ax.set_ylabel('lon')
        ax.set_zlabel('power')
        ax.set_zlim(-170,-80)
        ax.set_title(str(tim))
        plt.savefig(str(tim)+"plane.jpg",dpi=700)
        plt.close()
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
            direction=[[angle], [-1.5*np.linalg.norm(vector_1)]],
            pen="2p",
            color="red3",
        )
        os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making/primary_video")
        fig2.savefig(str(tim)+"direction.jpg")


def fig_merge(tim):
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making/primary_video")
    image1 = Image.open(str(tim)+"wave.jpg")
    image2 = Image.open(str(tim)+"plane.jpg")
    image3 = Image.open(str(tim)+"direction.jpg")
    image1 = image1.resize((1280,500))
    image2 = image2.resize((640, 450))
    image3 = image3.resize((640, 360))
    image1_size = image1.size
    new_image = Image.new('RGB',(image1_size[0], 2*image1_size[1]), (250,250,250))
    new_image.paste(image1,(0,0))
    new_image.paste(image2,(0,10+image1_size[1]))
    new_image.paste(image3,(int(image1_size[0]/2),10+image1_size[1]))
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
      res=final[29:34,j]
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
    intg=trapz(res,freq[29:34])
    return intg
    




def convert_pictures_to_video(pathIn, pathOut, fps, time):
    ''' this function converts images to video'''
    pics=glob.glob(join(directory+"/*wave.jpg"))
    len(pics)
    for ele in glob.glob(join(directory+"/*plane.jpg")):
        pics.append(ele)
    for ele in glob.glob(join(directory+"/*direction.jpg")):
        pics.append(ele)
    for ele in pics:
        os.remove(ele)
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
directory="//Users/sebinjohn/AON_PROJECT/Data/Video Making/primary_video"
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





import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'grib',
        'variable': 'significant_height_of_combined_wind_waves_and_swell',
        'year': '2019',
        'day': '01',
        'time': '01:00',
        'month': '01',
    },
    'download.grib')