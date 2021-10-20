import datetime as dt
import matplotlib.dates as mdates
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import dill
from datetime import datetime
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
import dill
import os
#----------------------------------------------------------------------------------------
os.chdir("/home/sjohn/AON_PROJECT/O14K")
dill.load_session("notebook_env.db")

with open('final.npy', 'rb') as g:
    final=np.load(g)
num_rows, num_cols =final.shape

startt=UTCDateTime(starttimes[0])
endt=UTCDateTime(endtimes[-1])
l=((endt-startt)/num_cols)
print(l)
time=np.arange(startt,endt,l)
np.shape(time)
pmin=(np.amin(final))
pmax=np.amax(final)

date_format=mdates.DateFormatter('%d,%b,%y')
date_format1=mdates.DateFormatter('%d,%b,%y')
date_format

freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)

with open('finalp.npy', 'rb') as g:
    finalp=np.load(g)
with open('finalw.npy', 'rb') as g:
    finalw=np.load(g)

fint_frames=np.arange(starttimeta,endtimeak,3600)
for i in range(len(fint_frames)):
    fint_frames[i]=mdates.date2num(fint_frames[i])
mintp=min(fint_frames)
maxtp=max(fint_frames)
#-----------------------------------------------------------------------------------------
j=53
xmin=18400
xmax=18600


tmin=mdates.num2date(xmin)
tmax=mdates.num2date(xmax)
ans=UTCDateTime(k)-UTCDateTime(startt)
startcol=int(ans/(60*60))
endcol=int(startcol+(UTCDateTime(tmax)-UTCDateTime(tmin))/(60*60))
print(startcol,endcol)
powe=final[j,startcol:endcol]
startcolw=(np.where(finalw==xmin))
endcolw=(np.where(finalw==xmax))
wspd=finalw[(int(startcolw[0])):(int(endcolw[0]))]
startcolp=(np.where(finalp==xmin))
endcolp=(np.where(finalp==xmax))
prsr=finalp[(int(startcolp[0])):(int(endcolp[0]))]
len(prsr)

#----------------------------------------------------------------------------


UTCDateTime(tmax)






