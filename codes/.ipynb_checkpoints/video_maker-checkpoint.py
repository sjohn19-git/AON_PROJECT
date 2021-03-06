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
      "-169.8951"
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
     "52.8235"
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
          "final_CLES.npy"
          
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
     "notebook_env_CLES.db"
     ]
os.chdir("/home/sjohn/AON_PROJECT/Video Making")
with open(("metadta.pkl"), 'wb') as f:  # Python 3: open(..., 'wb')
          pickle.dump([long,lat,stationo ,env], f)

freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)
len(freq)


pygmt.makecpt(cmap="hot",output="exp.cpt", series=[-150,-90,0.005],reverse=True)
grid = pygmt.datasets.load_earth_relief(resolution="10m", region=[-172, -135, 52, 73])
datapath="/home/sjohn/Data/*/"
st=UTCDateTime(2019,11,1)
et=UTCDateTime(2019,12,1)
windw=2
tim_len=int((et-st)/(3600*windw))
smth_intg=np.zeros((len(stationo),tim_len))


def map_plo(st,et,co):
    time_frame=[]
    for i in range (int((et-st)/(3600*windw))):
        time_frame.append(st)
        st=st+(3600*windw)
    for i in range(len(time_frame)):
        tim=time_frame[i]
        print(tim)
        integral=text_maker(tim,co)
        co+=1
        grid,data_wave=wave(tim)
        ma_pic_generator(tim,grid,data_wave)

        
def wave(tim):
    c = cdsapi.Client()
    os.chdir("/home/sjohn/Data/wave")
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
    
    

def ma_pic_generator(tim,grid,data_wave):
    os.chdir("/home/sjohn/AON_PROJECT/Video Making")
    data= pd.read_csv("mapplo.csv")
    fig=pygmt.Figure()
    with fig.subplot(nrows=1, ncols=2, figsize=("45c", "60c"), frame="lrtb",clearance='2c',title=""):
        with fig.set_panel(panel=0): 
            #pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[np.amin(data_wave),np.amax(data_wave)],reverse=True)
            fig.basemap(frame=True,region=[150,260,30,73], projection="L-159/35/33/45/22c")
            fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/45/22c",shorelines=True)
           # fig.grdimage(grid=grid,projection="S200/90/20c", cmap="expll1.cpt",nan_transparent=True)
            fig.colorbar(cmap="exp.cpt",frame='af+l"PSD (db)"',projection="L-159/35/33/45/22c")
            fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="L-159/35/33/45/22c", cmap="expll1.cpt",frame="a",nan_transparent=True)
            fig.plot(x=data['long'],y=data['lat'],color=data['powe'],style="c0.35c",cmap="exp.cpt",pen="black", projection="L-159/35/33/45/22c")
            fig.text(text=str(tim),x=-160,y=72,projection="L-159/35/33/45/22c")
        with fig.set_panel(panel=1): 
            pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[np.amin(data_wave),np.amax(data_wave)],reverse=True)
            fig.coast(region="g", projection="Cyl_stere/180/-40/20c",shorelines=True)
            fig.grdimage(grid=grid,projection="Cyl_stere/180/-40/20c",frame="a", cmap="expll1.cpt",nan_transparent=True)
            fig.colorbar(projection="Cyl_stere/180/-40/20c", cmap="expll1.cpt",frame='af+l"Significant Wave Height (m)"')
            fig.text(text=str(tim),x=190,y=85,projection="Cyl_stere/180/-40/20c")
        fig.show()
        #fig.show(method="external")
        fig.savefig(str(tim)+".jpg")


def text_maker(tim,co):
  fields=["long","lat","powe"]
  rows=[]
  for i in range(len(long)):
      os.chdir(join("/home/sjohn/Data",env[i].split("_")[-1].split(".")[-2]))
      sta=env[i].split("_")[-1].split(".")[-2]
      with open((str(sta)+".pkl"),"rb") as f:
         sta, starttimeta, endtimeta,starttimeak,endtimeak,sta,cha,loc = pickle.load(f)
      j=int((tim-starttimeta)/(3600))
      with open((glob.glob(join(datapath+stationo[i]))[0]), 'rb') as g:
          final=np.load(g)
      for z in range(3): 
          res=final[34:42,j+z]
          mea=mean(res)
          intg=integ(res)
          if z==0:
              smooth_d1=mea
              intg_d1=intg
          elif z==1:
              smooth_d2=mea
              intg_d2=intg
          else:
              smooth_d3=mea
              intg_d3=intg
      smooth=min(smooth_d1,smooth_d2,smooth_d3) 
      smthd_intgeral=min(intg_d1,intg_d2,intg_d3)
      if smooth ==0:
          smooth=-160
          print(stationo[i])
      if smthd_intgeral==0:
          smthd_intgeral=float("Nan")
      smth_intg[i,co]=smthd_intgeral
      rows.append([str(long[i]),str(lat[i]),str(smooth)])
  os.chdir("/home/sjohn/AON_PROJECT/Video Making")
  with open("mapplo.csv", 'w') as csvfile:
      csvwriter = csv.writer(csvfile) 
      csvwriter.writerow(fields)
      csvwriter.writerows(rows)
  return smth_intg
     
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
            print(frame_array)
    out=cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'VP80'), fps,size)
    for i in range(len(frame_array)):
        out.write(frame_array[i])
    out.release()

# Example:
directory="/home/sjohn/AON_PROJECT/Video Making"
pathIn=directory+'/*.jpg'
pathOut=directory+"/"+'video_EX8.mp4'
fps=8
time=1 # the duration of each picture in the video




map_plo(st,et,0)



convert_pictures_to_video(pathIn, pathOut, fps, time)



len(env)


#--------------------------------------------------------------------------------
pics=glob.glob(join(directory+"/*jpg"))
for ele in pics:
    os.remove(ele)
time_frame=[]
for i in range (int((et-st)/(3600*windw))):
    time_frame.append(st)
    st=st+(3600*windw)  
 tim=time_frame[2]
c = cdsapi.Client()
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

