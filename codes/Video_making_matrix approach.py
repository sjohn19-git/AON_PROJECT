import os 
import pygmt
#import dill
import csv
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

long=["-161.2638","-133.121","-136.232","-154.4524",
      "-131.615","-147.4031","-156.4394","-158.9593",
      "-161.96","-156.6886","-158.2972","-159.4747",
      "-162.6829","-164.6555","-159.9126","-160.6514",
      "-153.483","-159.6514","-164.6483","-160.6027",
      "-163.0831","-143.7114","-144.9122","-150.6126",
      "-165.3436","-153.4196","-161.8016","-174.197495",
      "-154.9742","-156.6175","-161.0713","-159.0777",
      "-154.1467","-145.3697","-155.6214","-141.6153",
      "-149.3603","-143.1541","-151.9822","-147.8781",
      "-151.3773","-149.5432","-151.3773","-152.805",
      "-154.8808","-162.239","-141.878","-154.3915",
      "-147.7262","-142.0758","-145.7784","-154.07",
      "-153.8798","-152.1894","-150.9232","-141.6622",
      "-142.4931","-145.469894","-145.1568","-146.3399",
      "-155.208","-152.624","-155.8897","-157.5742",
      "-155.6225","-163.417603","-166.504501","-171.703598",
      "-157.1599","-159.5874","-154.7833","-148.937302",
      "-152.6821","-148.8233","-153.9721","-151.8132",
      "-148.4868","-145.568","-141.5951","-150.0239",
      "-147.4754","-141.6549","-151.5317","-134.5763",
      "-150.7703","-157.2316","-157.087","-160.0921",
      "-156.879196","-157.9906","-153.6446","-161.5308",
      "-150.223099","-145.7749","-148.860138","-150.741394",
      "-145.234207","-145.234207","-161.1943","-150.746704",
      "-139.6369","-153.1318","-143.2841","-146.3399",
      "-169.5484","-159.589493","-166.2011","-143.092606",
      "-175.103104","-177.1296","178.5112","-178.039",
      "-169.8951","-162.600006","-155.7251","-156.6132",
      "-139.5349","-137.8969","-136.7191","-139.8717",
      "-135.7863","-133.742","-138.0223","-136.2216",
      "-134.2708","-138.3689","-134.3426","-139.9355",
      "-138.3063","-136.3767","-138.2164","-136.3304",
      "-137.5201","-138.129","178.567","-138.4624",
      "-136.7935","-134.3906","-137.0885","-135.7796",
      "-138.7367","-140.1906","-138.5755","-136.0906",
      "-137.7381","-136.9598","-133.7147","-132.8174",
      "-132.2691","-130.9673","-131.1312","-130.2496",
      "-130.0257","-133.0818"
      ]
lat=["59.2533","56.1146","57.9616","57.5665",
     "55.3279","59.9979","59.1953","61.0224",
     "60.7515","62.2195","62.1344",'61.7105',
     "61.3416",'61.9331',"64.937","65.5011",
     "67.2221","66.6001","65.7077","67.4213",
     "67.6988","69.626","69.9175","69.836",
     "68.2746","70.34","69.3641","52.201599",
     "71.0033","71.3221","70.2043","63.3965",
     "64.1767","64.613","63.994","65.6035",
     "65.1479","65.3064","65.18","65.8371",
     "65.8937","65.8251","65.8937","65.6571",
     "65.4924","64.6379","62.3579","61.9037",
     "61.7929","64.0292","63.8036","63.3569",
     "62.4787","62.896198","63.5527","61.0595",
     "60.403","62.9698","62.3987","59.4296",
     "59.8542","60.0815","60.6801","57.5673",
     "55.8218","54.8564","53.8452","63.775799",
     "70.0079","69.1049","69.1565","63.7318",
     "68.8799","69.1532","68.4414","68.1343",
     "68.0748","68.1207","68.1861","66.7108",
     "66.7004","66.8088","61.98","57.4688",
     "62.5258","67.4572","66.1434","60.1686",
     "64.746101","59.0314","58.9287","62.2938",
     "60.5117","61.1292","63.405655","59.741699",
     "66.565697","66.565697","68.6483","61.4636",
     "59.9534","61.8823","60.1205","59.4296",
     "56.6011","55.831001","60.3849","60.7523",
     "52.730801","51.9303","51.9484","51.8339",
     "52.8235","66.895103","67.0486","68.7132",
     "68.6043","68.3889","66.3701","67.6136",
     "67.6106","67.441","66.9116","66.9808",
     "66.9227","66.2191","65.8052","65.4483",
     "65.3609","65.2225","64.4525","64.5753",
     "63.8433","63.109","51.932","62.4435",
     "62.5763","62.2024","61.4593","61.4817",
     "69.3286","60.7718","60.3024","60.7704",
     "59.6304","60.1218","59.5898","60.2114",
     "58.9601","59.3946","57.9128","56.9811",
     "55.9154","61.1512"
     ]
stationo=["final_O14K.npy","final_U33K.npy","final_S31K.npy","final_R18K.npy",
          "final_V35K.npy","final_P23K.npy","final_P17K.npy","final_M16K.npy",
          "final_M14K.npy","final_L18K.npy","final_L17K.npy","final_L16K.npy",
          "final_L14K.npy","final_K13K.npy","final_H17K.npy","final_G17K.npy",
          "final_F21K.npy","final_F18K.npy","final_F15K.npy","final_E18K.npy",
          "final_D17K.npy","final_C27K.npy","final_C26K.npy","final_C23K.npy",
          "final_C16K.npy","final_B22K.npy","final_B18K.npy","final_ATKA.npy",
          "final_A22K.npy","final_A21K.npy","final_A19K.npy","final_J17K.npy",
          "final_J20K.npy","final_J25K.npy","final_J19K.npy","final_I27K.npy",
          "final_I23K.npy","final_I26K.npy","final_I21K.npy","final_H24K.npy",
          "final_H22K.npy","final_H23K.npy","final_H22K.npy","final_H21K.npy",
          "final_H20K.npy","final_H16K.npy","final_M27K.npy","final_M19K.npy",
          "final_M23K.npy","final_K27K.npy","final_K24K.npy","final_K20K.npy",
          "final_L20K.npy","final_PPLA.npy","final_KTH.npy","final_BARN.npy",
          "final_BARK.npy","final_PAX.npy","final_HARP.npy","final_Q23K.npy",
          "final_O18K.npy","final_O20K.npy","final_N18K.npy","final_R16K.npy",
          "final_CHI.npy","final_FALS.npy","final_UNV.npy","final_GAMB.npy",
          "final_B20K.npy","final_C19K.npy","final_C21K.npy","final_MCK.npy", 
          "final_D22K.npy","final_D24K.npy","final_E21K.npy","final_E22K.npy",
          "final_E24K.npy","final_E25K.npy","final_E27K.npy","final_G23K.npy",
          "final_G24K.npy","final_G27K.npy","final_SKN.npy","final_S32K.npy",
          "final_L22K.npy","final_E19K.npy","final_G19K.npy","final_N15K.npy",
          "final_GCSA.npy","final_P16K.npy","final_Q19K.npy","final_K15K.npy",
          "final_SLK.npy","final_DIV.npy","final_RND.npy","final_BRSE.npy",
          "final_FYU.npy","final_TGL.npy","final_C18K.npy","final_SSN.npy",
          "final_BCP.npy","final_M20K.npy","final_BGLC.npy","final_Q23K.npy",
          "final_P08K.npy","final_CHN.npy","final_M11K.npy","final_CRQ.npy",
          "final_SMY.npy","final_KINC.npy","final_LSSA.npy","final_TASE.npy",
          "final_CLES.npy","final_KOTZ.npy","final_F20K.npy","final_D20K.npy",
          "final_E28M.npy","final_E29M.npy","final_EPYK.npy","final_F28M.npy",
          "final_F30M.npy","final_F31M.npy","final_G29M.npy","final_G30M.npy",
          "final_G31M.npy","final_H29M.npy","final_H31M.npy","final_I28M.npy",
          "final_I29M.npy","final_I30M.npy","final_J29N.npy","final_J30M.npy",
          "final_K29M.npy","final_L29M.npy","final_LSSE.npy","final_M29M.npy",
          "final_M30M.npy","final_M31M.npy","final_N30M.npy","final_N31M.npy",
          "final_D28M.npy","final_O28M.npy","final_O29M.npy","final_O30N.npy",
          "final_P29M.npy","final_P30M.npy","final_P32M.npy","final_P33M.npy",
          "final_Q32M.npy","final_R33M.npy","final_S34M.npy","final_T35M.npy",
          "final_U35K.npy","final_N32M.npy"
          ]
env=["notebook_env_O14K.db","notebook_env_U33K.db","notebook_env_S31K.db","notebook_env_R18K.db",
     "notebook_env_V35K.db","notebook_env_P23K.db","notebook_env_P17K.db","notebook_env_M16K.db",
     "notebook_env_M14K.db","notebook_env_L18K.db","notebook_env_L17K.db","notebook_env_L16K.db",
     "notebook_env_L14K.db","notebook_env_K13K.db","notebook_env_H17K.db","notebook_env_G17K.db",
     "notebook_env_F21K.db","notebook_env_F18K.db","notebook_env_F15K.db","notebook_env_E18K.db",
     "notebook_env_D17K.db","notebook_env_C27K.db","notebook_env_C26K.db","notebook_env_C23K.db",
     "notebook_env_C16K.db","notebook_env_B22K.db","notebook_env_B18K.db","notebook_env_ATKA.db",
     "notebook_env_A22K.db","notebook_env_A21K.db","notebook_env_A19K.db","notebook_env_J17K.db",
     "notebook_env_J20K.db","notebook_env_J25K.db","notebook_env_J19K.db","notebook_env_I27K.db",
     "notebook_env_I23K.db","notebook_env_I26K.db","notebook_env_I21K.db","notebook_env_H24K.db",
     "notebook_env_H22K.db","notebook_env_H23K.db","notebook_env_H22K.db","notebook_env_H21K.db",
     "notebook_env_H20K.db","notebook_env_H16K.db","notebook_env_M27K.db","notebook_env_M19K.db",
     "notebook_env_M23K.db","notebook_env_K27K.db","notebook_env_K24K.db","notebook_env_K20K.db",
     "notebook_env_L20K.db","notebook_env_PPLA.db","notebook_env_KTH.db","notebook_env_BARN.db",
     "notebook_env_BARK.db","notebook_env_PAX.db","notebook_env_HARP.db","notebook_env_Q23K.db",
     "notebook_env_O18K.db","notebook_env_O20K.db","notebook_env_N18K.db","notebook_env_R16K.db",
     "notebook_env_CHI.db","notebook_env_FALS.db","notebook_env_UNV.db","notebook_env_GAMB.db",
     "notebook_env_B20K.db","notebook_env_C19K.db","notebook_env_C21K.db","notebook_env_MCK.db",
     "notebook_env_D22K.db","notebook_env_D24K.db","notebook_env_E21K.db","notebook_env_E22K.db",
     "notebook_env_E24K.db","notebook_env_E25K.db","notebook_env_E27K.db","notebook_env_G23K.db",
     "notebook_env_G24K.db","notebook_env_G27K.db","notebook_env_SKN.db","notebook_env_S32K.db",
     "notebook_env_L22K.db","notebook_env_E19K.db","notebook_env_G19K.db","notebook_env_N15K.db",
     "notebook_env_GCSA.db","notebook_env_P16K.db","notebook_env_Q19K.db","notebook_env_K15K.db",
     "notebook_env_SLK.db","notebook_env_DIV.db","notebook_env_RND.db","notebook_env_BRSE.db",
     "notebook_env_FYU.db","notebook_env_TGL.db","notebook_env_C18K.db","notebook_env_SSN.db",
     "notebook_env_BCP.db","notebook_env_M20K.db","notebook_env_BGLC.db","notebook_env_Q23K.db",
     "notebook_env_P08K.db","notebook_env_CHN.db","notebook_env_M11K.db","notebook_env_CRQ.db",
     "notebook_env_SMY.db","notebook_env_KINC.db","notebook_env_LSSA.db","notebook_env_TASE.db",
     "notebook_env_CLES.db","notebook_env_KOTZ.db","notebook_env_F20K.db","notebook_env_D20K.db",
     "notebook_env_E28M.db","notebook_env_E29M.db","notebook_env_EPYK.db","notebook_env_F28M.db",
     "notebook_env_F30M.db","notebook_env_F31M.npy","notebook_env_G29M.npy","notebook_env_G30M.npy",
     "notebook_env_G31M.npy", "notebook_env_H29M.npy","notebook_env_H31M.npy","notebook_env_I28M.npy",
     "notebook_env_I29M.npy","notebook_env_I30M.npy","notebook_env_J29N.npy","notebook_env_J30M.npy",
     "notebook_env_K29M.npy", "notebook_env_L29M.npy","notebook_env_LSSE.npy","notebook_env_M29M.npy",
     "notebook_env_M30M.npy","notebook_env_M31M.npy","notebook_env_N30M.npy","notebook_env_N31M.npy",
     "notebook_env_D28M.npy", "notebook_env_O28M.npy","notebook_env_O29M.npy","notebook_env_O30N.npy",
     "notebook_env_P29M.npy","notebook_env_P30M.npy","notebook_env_P32M.npy","notebook_env_P33M.npy",
     "notebook_env_Q32M.npy","notebook_env_R33M.npy","notebook_env_S34M.npy","notebook_env_T35M.npy",
     "notebook_env_U35K.npy","notebook_env_N32M.npy"
     ]



os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
with open(("metadta.pkl"), 'wb') as f:  # Python 3: open(..., 'wb')
          pickle.dump([long,lat,stationo ,env], f)

freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)
len(freq)


pygmt.makecpt(cmap="hot",output="exp.cpt", series=[-150,-90,0.005],reverse=True)
pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[0,14],reverse=True)
grid = pygmt.datasets.load_earth_relief(resolution="10m", region=[-172, -135, 52, 73])
datapath="/Users/sebinjohn/AON_PROJECT/Data/*/"
st=UTCDateTime(2019,10,1)
et=UTCDateTime(2019,11,1)
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
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
    fig=pygmt.Figure()
    with fig.subplot(nrows=1, ncols=2, figsize=("45c", "60c"), frame="lrtb",clearance='2c',title=""):
        with fig.set_panel(panel=0): 
            #fig.basemap(frame=True,region="g", projection="L-159/35/33/45/22c")
            fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="expll1.cpt",frame="a",nan_transparent=True)
            fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/22c",shorelines=True)
            fig.colorbar(cmap="exp.cpt",frame='af+l"PSD (db)"',projection="L-159/35/33/45/22c")  
            fig.plot(x=median_timeseries[1,1:,z+1],y=median_timeseries[2,1:,z+1],color=median_timeseries[0,1:,z+1],style="c0.30c",cmap="exp.cpt",pen="black", projection="L-159/35/33/45/22c")
            fig.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/22c")
        with fig.set_panel(panel=1): 
            fig.coast(region="g", projection="Cyl_stere/180/-40/20c",shorelines=True)
            fig.grdimage(grid=grid,projection="Cyl_stere/180/-40/20c",frame="a", cmap="expll1.cpt",nan_transparent=True)
            fig.colorbar(projection="Cyl_stere/180/-40/20c", cmap="expll1.cpt",frame='af+l"Significant Wave Height (m)"')
            fig.text(text=str(tim),x=190,y=85,projection="Cyl_stere/180/-40/20c")
        fig.show()
        fig.savefig(str(tim)+".jpg")
        print("saved "+ str(tim)+".jpg")


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
    frame_array=[]
    files=[f for f in glob.glob(pathIn) if isfile(join(pathIn,f))]
    files.sort()
    files
    for i in range (len(files)):
        '''reading images'''
        img=cv2.imread(files[i])
        height, width, layers = img.shape
        size=(width,height)
        for k in range (time):
            frame_array.append(img)
    out=cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'VP80'), fps,size)
    for i in range(len(frame_array)):
        out.write(frame_array[i])
    out.release()

# Example:
directory="//Users/sebinjohn/AON_PROJECT/Data/Video Making"
pathIn=directory+'/*.jpg'
pathOut=directory+"/"+str(st)+"-"+str(et)+".avi"
fps=8
time=1 # the duration of each picture in the video




mean_timeseries,median_timeseries=map_plo(st,et)



convert_pictures_to_video(pathIn, pathOut, fps, time)



len(env)


#--------------------------------------------------------------------------------
pics=glob.glob(join(directory+"/*jpg"))
for ele in pics:
    os.remove(ele)
    os.remove(ele)
#--------------------------------------------------------------------------------
time=[]
for ele in time_frame:
    time.append(ele.matplotlib_date)
date_format=mdates.DateFormatter('%d,%b,%y')

fig,ax1=plt.subplots(3,4,figsize=(40,25))

for i in range(1,5):
    ax1[0][i-1].scatter(time,mean_timeseries[0,i,1:],color="blue",label="raw data")
    ax1[0][i-1].scatter(time,median_timeseries[0,i,1:],color="red",label="processed data")
    ax1[0][i-1].legend()
    ax1[0][i-1].xaxis.set_major_formatter(date_format)
    ax1[0][i-1].set_ylabel("power in db")
    ax1[0][i-1].set_title(stationo[i].split("_")[-1].split('.')[-2])
for i in range(5,9):
    ax1[1][i-5].scatter(time,mean_timeseries[0,i,1:],color="blue",label="raw data")
    ax1[1][i-5].scatter(time,median_timeseries[0,i,1:],color="red",label="processed data")
    ax1[1][i-5].legend()
    ax1[1][i-5].xaxis.set_major_formatter(date_format)
    ax1[1][i-5].set_ylabel("power in db")
    ax1[1][i-5].set_title(stationo[i].split("_")[-1].split('.')[-2])
for i in range(9,13):
    ax1[2][i-9].scatter(time,mean_timeseries[0,i,1:],color="blue",label="raw data")
    ax1[2][i-9].scatter(time,median_timeseries[0,i,1:],color="red",label="processed data")
    ax1[2][i-9].legend()
    ax1[2][i-9].xaxis.set_major_formatter(date_format)
    ax1[2][i-9].set_ylabel("power in db")
    ax1[2][i-9].set_title(stationo[i].split("_")[-1].split('.')[-2])

