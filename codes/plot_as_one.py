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
#_--------------------------------------------------------------------------------------------
dill.load_session("notebook_env.db")
fint_frames=np.arange(starttimeta,endtimeak,3600)
for i in range(len(fint_frames)):
    fint_frames[i]=mdates.date2num(fint_frames[i])
#-----------------------------------------------------------------------------------------------------
xmin=18383
xmax=18628
xlmi=xmin
xlma=xmax
j=42

cmap = plt.get_cmap('nipy_spectral').copy()
cmap.set_over(color = 'w')
fig = plt.figure(constrained_layout=True)
fig.set_figheight(20)
fig.set_figwidth(12)
spec = gridspec.GridSpec(5,1,figure=fig)
ax1 = fig.add_subplot(spec[0])
c = ax1.pcolormesh(fint_frames,freq,final,cmap=cmap,vmin=-195,vmax=-95,shading='auto',rasterized=True)

nom=cl.Normalize(vmin=-150,vmax=-95)
ax1.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
#ax1.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=(WE)))
ax1.xaxis.set_major_formatter(date_format)
#ax1.xaxis.set_minor_formatter(date_format1)

ax1.set_ylabel("Frequency(Hz)")
ax1.set_xlim([xmin,xmax])
for label in ax1.get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')
ax1.set_yscale('log')
ax1.title.set_text(sta+"_"+cha+" "+str(UTCDateTime((mdates.num2date(xmin)).strftime('%Y-%m-%d')))+" - "+str(UTCDateTime((mdates.num2date(xmax)).strftime('%Y-%m-%d'))))

#---------------------------------------------------------------------------------------
ax2 = fig.add_subplot(spec[1])


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
len(freq)
with open('finalp.npy', 'rb') as g:
    finalp=np.load(g)
with open('finalw.npy', 'rb') as g:
    finalw=np.load(g)
fint_frames=np.arange(starttimeta,endtimeak,3600)
for i in range(len(fint_frames)):
    fint_frames[i]=mdates.date2num(fint_frames[i])
mintp=min(fint_frames)
maxtp=max(fint_frames)

ax2.set_ylim(-80,-160)

p1 = ax2.scatter(fint_frames[:],final[j,:],s=0.5,label="power (db)",c="red")
ax2.xaxis.set_major_formatter(date_format)
ax2.xaxis.set_minor_formatter(date_format1)
ax2.set_xlim([xmin,xmax])
ax2.set_ylabel("power ("+str(freq[j])+"Hz)")
ax2.title.set_text(sta+"_"+cha+" "+str(UTCDateTime((mdates.num2date(xmin)).strftime('%Y-%m-%d')))+" - "+str(UTCDateTime((mdates.num2date(xmax)).strftime('%Y-%m-%d'))))
ax2.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
ax2.legend()
for label in ax2.get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')

#---------------------------------------------------------------------------------------

ax3=fig.add_subplot(spec[2])

p1 = ax3.scatter(finalp[:,0],finalp[:,1], s=0.5,label="pressure (atm)",c="blue")
ax3.xaxis.set_major_formatter(date_format)
ax3.xaxis.set_minor_formatter(date_format1)
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([0.8,1])
ax3.set_ylabel("Pressure (atm)"+ str(str(freq[j])+"Hz"))
ax3.title.set_text(sta+"_"+cha+" "+str(UTCDateTime((mdates.num2date(xmin)).strftime('%Y-%m-%d')))+" - "+str(UTCDateTime((mdates.num2date(xmax)).strftime('%Y-%m-%d'))))
ax3.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
ax3.legend()
for label in ax3.get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')
#--------------------------------------------------------------------------------------

import obspy
from obspy import read 
from medfilt import medfilt
import numpy as np
import obspy.signal.cross_correlation as cr
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import datetime
from nan_help import nan_helper


xmin=18375
xmax=18634

shift=1
window=170
file_name=input("File_name?")
freuer=input("freq")
st0=read("/home/sjohn/AON_PROJECT/O14K/"+file_name)
time=xmin+(window/48)
tr1=st0[0]
tr2=st0[2]
arr=tr1.data
arr[arr==0]=np.nan
nans, x= nan_helper(arr)
arr[nans]= np.interp(x(nans), x(~nans), arr[~nans])
tr2starttime=tr2.stats.starttime
arr2=tr2.data
arr2[arr2==0]=np.nan
arr2[arr2<0.85]=np.nan
arr2[arr2>0.96]=np.nan
nans, x= nan_helper(arr2)
arr2[nans]= np.interp(x(nans), x(~nans), arr2[~nans])
arr2=medfilt(arr2,11)
tr2.data=arr2
st=st0.copy()
st3=st0.copy()







dt=tr1.stats.starttime
starttime=tr1.stats.starttime
endtime=tr1.stats.endtime
starttimep=tr2.stats.starttime
endtimep=tr2.stats.endtime
#tr1.plot()

fint_frames=np.arange(starttime,endtime,3600)
for i in range(len(fint_frames)):
    fint_frames[i]=mdates.date2num(fint_frames[i])

window/2

k=window
corre=np.array([])
index=np.array([])
lag=np.array([])
for i in range(len(tr1)):
    try:
        st.trim(dt, dt+3600*k)
        tr2=st[2]
        tr1=st[0]
        tr2.detrend(type="demean")
        tr1.detrend(type="demean")
        tr2.detrend(type="linear")
        tr1.detrend(type="linear")
        tr1.taper(0.05, type='cosine',side='both')
        tr2.taper(0.05,type='cosine',side='both')
        tr1par=tr1.data
        tr2par=tr2.data
        tr1par=medfilt(tr1par,11)
        coeif=cr.correlate(tr1par, tr2par,int(window/2))
        corre=np.append(corre,np.amin(coeif))
        index=np.append(index,np.array((dt+1800*k)))
        lag=np.append(lag,np.where(coeif==np.amin(coeif))[0][0]-int((window/2)+1))
        st=st0.copy() 
        dt+=3600
    except:
        break
        print("error")

np.amin(coeif)
len(lag)
np.where(coeif==np.amin(coeif))[0][0]



for i in range(len(index)):
    index[i]=index[i].matplotlib_date

#--------------------------------------------------------------------------------------
ax4=fig.add_subplot(spec[3])
ax4.scatter(index,corre,s=0.2)
ax4.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
date_format=mdates.DateFormatter('%d,%b,%y')
ax4.xaxis.set_major_formatter(date_format)
for label in ax4.get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')


ax4.set_xlabel("Time")
ax4.set_xlim([xlmi,xlma])
ax4.set_ylabel("cross correlation coeifficient")
ax4.title.set_text("cross correlation in window length "+str(round((window/24),2))+"days "+str(UTCDateTime(mdates.num2date(xlmi)))+" - "+str(UTCDateTime(mdates.num2date(xlma)))+" of "+freuer+"Hz")
#-----------------------------------------------------------------------------------------------------
ax5=fig.add_subplot(spec[4])

ax5.scatter(index,lag,s=0.2)

ax5.set_xlim([xlmi,xlma])
ax5.title.set_text("cross correlation lag time in hours  "+str(UTCDateTime(mdates.num2date(xlmi)))+" - "+str(UTCDateTime(mdates.num2date(xlma)))+" of "+freuer+"Hz")
ax5.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
date_format=mdates.DateFormatter('%d,%b,%y')
plt.ylabel("Lag in hours")
ax5.xaxis.set_major_formatter(date_format)
for label in ax5.get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')
fig.savefig(freuer+"full_plot.png")
