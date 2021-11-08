import os 
import pygmt
import dill
import csv
import pandas as pd
from obspy import UTCDateTime
import cv2
from os.path import isfile, join 
import glob

long=["-161.2638","-133.121",
      "-131.615"
      ]
lat=["59.2533","56.1146",
     "55.3279"
     ]
stationo=["final_O14K.npy","final_U33K.npy",
          "final_V35K.npy"
          ]
env=["notebook_env_O14K.db","notebook_env_U33K.db","notebook_env_V35K.db"]
pygmt.makecpt(cmap="rainbow",output="exp.cpt", series=[-135,-100,0.005])
grid = pygmt.datasets.load_earth_relief(resolution="10m", region=[-172, -135, 52, 73])
os.chdir("/Users/sebinjohn/AON_PROJECT/video making")

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
