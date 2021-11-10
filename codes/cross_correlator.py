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
freq=input("freq")
st0=read("/home/sjohn/AON_PROJECT/O14K/"+file_name)
time=xmin+(window/48)
tr1=st0[0]
tr2=st0[2]
st0[0].plot()
st0[2].plot()
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
st0[2].plot()
st0[0].plot()
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
        if i==2400:
            tr1.plot()
            tr2.plot()
            plt.figure(1)
            ax0=plt.gca()
            ax0.set_xlim([0,window])
            plt.plot(tr2par)
            plt.figure(2)
            ax0=plt.gca()
            plt.title("Median filtering detrended and tapered PSD")
            ax0.set_xlim([0,window])
            plt.plot(tr1.data,label="data")
            plt.plot(tr1par,label="median_filt")
            plt.xlabel("Time in hour")
            plt.legend()
        coeif=cr.correlate(tr1par, tr2par,int(window/2))
        corre=np.append(corre,np.amin(coeif))
        index=np.append(index,np.array((dt+1800*k)))
        lag=np.append(lag,np.where(coeif==np.amin(coeif))[0][0]-int((window/2)+1))
        st=st0.copy() 
        dt+=3600
    except:
        break
        print("error")
    
for ele in lag:
    print(ele)
np.amin(coeif)
len(lag)
np.where(coeif==np.amin(coeif))[0][0]



for i in range(len(index)):
    index[i]=index[i].matplotlib_date

#--------------------------------------------------------------------------------------

xlmi=18383
xlma=18628

plt.figure(3)
plt.scatter(index,corre,s=0.2)
fig=plt.gcf()
fig.set_figwidth(18)
fig.set_figheight(4)
ax = plt.gca()
ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
date_format=mdates.DateFormatter('%d,%b,%y')
ax.xaxis.set_major_formatter(date_format)
for label in ax.get_xticklabels():
    label.set_rotation(20)
    label.set_horizontalalignment('right')


plt.xlabel("Time")
ax.set_xlim([xlmi,xlma])
plt.ylabel("cross correlation coeifficient")
plt.title("cross correlation in window length "+str(round((window/24),2))+"days "+str(UTCDateTime(mdates.num2date(xlmi)))+" - "+str(UTCDateTime(mdates.num2date(xlma)))+" of "+freq+"Hz")

plt.figure(4)

plt.scatter(index,lag,s=0.2)
fig=plt.gcf()
fig.set_figwidth(18)
fig.set_figheight(4)
ax = plt.gca()
ax.set_xlim([xlmi,xlma])
plt.title("cross correlation lag time in hours  "+str(UTCDateTime(mdates.num2date(xlmi)))+" - "+str(UTCDateTime(mdates.num2date(xlma)))+" of "+freq+"Hz")
ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
date_format=mdates.DateFormatter('%d,%b,%y')
plt.ylabel("Lag in hours")
ax.xaxis.set_major_formatter(date_format)
for label in ax.get_xticklabels():
    label.set_rotation(20)
    label.set_horizontalalignment('right')
