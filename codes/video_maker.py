import os 
import pygmt
import dill
import csv
import pandas as pd
from obspy import UTCDateTime
import cv2
from os.path import isfile, join 
import glob

long=["-161.2638","-133.121","-136.232","-154.4524",
      "-131.615","-147.4031","-156.4394","-158.9593",
      "-161.96","-156.6886","-158.2972","-159.4747",
      "-162.6829","-164.6555","-159.9126","-160.6514",
      "-153.483","-159.6514","-164.6483","-160.6027",
      "-163.0831"
      ]
lat=["59.2533","56.1146","57.9616","57.5665",
     "55.3279","59.9979","59.1953","61.0224",
     "60.7515","62.2195","62.1344",'61.7105',
     "61.3416",'61.9331',"64.937","65.5011",
     "67.2221","66.6001","65.7077","67.4213",
     "-163.0831"
     ]
stationo=["final_O14K.npy","final_U33K.npy","final_S31K.npy","final_R18K.npy",
          "final_V35K.npy","final_P23K.npy","final_P17K.npy","final_M16K.npy",
          "final_M14K.npy","final_L18K.npy","final_L17K.npy","final_L16K.npy",
          "final_L14K.npy","final_K13K.npy","final_H17K.npy","final_G17K.npy",
          "final_F21K.npy","final_F18K.npy","final_F15K.npy","final_E18K.npy",
          "final_D17K.npy"
          
          ]
env=["notebook_env_O14K.db","notebook_env_U33K.db","notebook_env_S31K.db","notebook_env_R18K.db",
     "notebook_env_V35K.db","notebook_env_P23K.db","notebook_env_P17K.db","notebook_env_M16K.db",
     "notebook_env_M14K.db","notebook_env_L18K.db","notebook_env_L17K.db","notebook_env_L16K.db",
     "notebook_env_L14K.db","notebook_env_K13K.db","notebook_env_H17K.db","notebook_env_G17K.db",
     "notebook_env_F21K.db","notebook_env_F18K.db","notebook_env_F15K.db","notebook_env_E18K.db",
     "notebook_env_D17K.db"
     ]
pygmt.makecpt(cmap="rainbow",output="exp.cpt", series=[-135,-100,0.005])
grid = pygmt.datasets.load_earth_relief(resolution="10m", region=[-172, -135, 52, 73])
os.chdir("/home/sjohn/AON_PROJECT/Video Making")

st=UTCDateTime(2019,1,1)
et=st+(3600*600)


def map_plo(st,et):
    time_frame=[]
    for i in range (int((et-st)/(3600*6))):
        time_frame.append(st)
        st=st+(3600*6)
    for i in range(len(time_frame)):
        tim=time_frame[i]
        text_maker(tim)
        ma_pic_generator(tim)

def ma_pic_generator(tim):
    fig=pygmt.Figure()
    fig.basemap(frame=True,region=[-190,-130, 48, 73], projection="S200/90/20c")
    #fig.grdimage(grid='@earth_relief_01m.grd ',region=[-190,-130, 48, 73],projection="S200/90/20c", cmap="light.cpt",frame="a")
    fig.coast(region=[-190,-130, 48, 73], projection="S200/90/20c",shorelines=True)
    data= pd.read_csv("mapplo.csv")
    fig.plot(x=data['long'],y=data['lat'],color=data['powe'],style="c0.4c",cmap="exp.cpt",pen="black")
    fig.text(text=str(tim),x=-160,y=72)
    fig.colorbar(cmap="exp.cpt",frame='af+l"PSD (db)"')
    fig.show()
    fig.savefig(str(tim)+".jpg")



def text_maker(tim):
  fields=["long","lat","powe"]
  rows=[]
  for i in range(len(long)):
      dill.load_session(env[i])
      j=int((tim-starttimeta)/(3600))
      with open(stationo[i], 'rb') as g:
          final=np.load(g)
      rows.append([str(long[i]),str(lat[i]),str(final[50,j])])
  with open("mapplo.csv", 'w') as csvfile:
      csvwriter = csv.writer(csvfile) 
      csvwriter.writerow(fields)
      csvwriter.writerows(rows)

map_plo(st,et)


def convert_pictures_to_video(pathIn, pathOut, fps, time):
    ''' this function converts images to video'''
    frame_array=[]
    files=[f for f in glob.glob(pathIn) if isfile(join(pathIn,f))]
    files.sort()
    files
    for i in range (len(files)):
        '''reading images'''
        img=cv2.imread(files[i])
        img=cv2.resize(img,(1400,1000))
        height, width, layers = img.shape
        size=(width,height)
        for k in range (time):
            frame_array.append(img)
            print(frame_array)
    out=cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'avc1'), fps,size)
    for i in range(len(frame_array)):
        out.write(frame_array[i])
    out.release()

# Example:
directory="/Users/sebinjohn/AON_PROJECT/Video Making"
pathIn=directory+'/*.jpg'
pathOut=directory+"/"+'video_EX9.mp4'
fps=3
time=2# the duration of each picture in the video
convert_pictures_to_video(pathIn, pathOut, fps, time)
files
