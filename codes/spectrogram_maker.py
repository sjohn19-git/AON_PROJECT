import datetime as dt
import matplotlib.dates as mdates
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob
import re
import matplotlib.cm as cm
import matplotlib.colors as cl
from matplotlib import gridspec
import dill
from datetime import datetime
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import requests

os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")


with open("metadta.pkl","rb") as f:
   long,lat,stationo ,env = pickle.load(f)



sta="O14K"
os.chdir(r"/Users/sebinjohn/AON_PROJECT/Data/"+sta)
with open((str(sta)+".pkl"),"rb") as f:
    sta, starttimeta, endtimeta,starttimeak,endtimeak,sta,cha,loc=pickle.load(f)

fint_frames=np.arange(starttimeta,endtimeak,3600)
for i in range(len(fint_frames)):
    fint_frames[i]=mdates.date2num(fint_frames[i])
    
with open("final_"+sta+".npy", 'rb') as g:
    final=np.load(g)

xmin=18262
xmax=18628
xlmi=xmin
xlma=xmax

freq=[]
name= pd.read_xml("/Users/sebinjohn/AON_PROJECT/Data/Video Making/pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)
len(freq)

date_format=mdates.DateFormatter('%d-%B-%Y')
date_format1=mdates.DateFormatter('%m/%d/%yT%H')
date_format

cmap = plt.get_cmap('nipy_spectral').copy()
cmap.set_over(color = 'w')
fig = plt.figure(constrained_layout=True)
fig = plt.figure()
fig.set_figheight(6)
fig.set_figwidth(18)
spec = gridspec.GridSpec(ncols=2, nrows=1,width_ratios=[10, 1], wspace=-0.1,hspace=0.2, height_ratios=[1])
ax1 = fig.add_subplot(spec[0])
c = ax1.pcolormesh(fint_frames,freq,final,cmap=cmap,vmin=-195,vmax=-95,shading='auto',rasterized=True)
# ax1.axhline(y=0.1,c="r",linestyle="--")
# ax1.axhline(y=0.2,c="r",linestyle="--")
ax1.axhline(y=0.2,c="k",linestyle="--")
#ax1.axhline(y=0.07,c="k",linestyle="--")
nom=cl.Normalize(vmin=-195,vmax=-95)
#ax1.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA, SU)))
ax1.xaxis.set_major_locator(mdates.MonthLocator(1))
#ax1.xaxis.set_major_locator(mdates.HourLocator(interval=4))
ax1.xaxis.set_major_formatter(date_format)
ax1.xaxis.set_minor_formatter(date_format1)
#ax1.set_ylim(freq[29],0.1)
#ax1.set_ylabel("Frequency(Hz)")
ax1.set_xlim([xmin,xmax])
# for label in ax1.get_xticklabels():
#     label.set_rotation(40)
#     label.set_horizontalalignment('right')
ax1.set_yscale('log')
#ax1.title.set_text(sta+"_"+cha+" "+str(UTCDateTime((mdates.num2date(xmin)).strftime('%Y-%m-%d')))[0:10]+" - "+str(UTCDateTime((mdates.num2date(xmax)).strftime('%Y-%m-%d')))[0:10])
#ax1.title.set_text(sta+"_"+cha+" ")
ax0 = fig.add_subplot(spec[1])
ax0.axis('off')
fig.colorbar(cm.ScalarMappable(norm=nom, cmap="nipy_spectral"), ax=ax0)


