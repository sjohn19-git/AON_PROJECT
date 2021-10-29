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




shift=1
window=170
st0=read("/home/sjohn/AON_PROJECT/O14K/O14K_2020-05-18T00:00:00.000000Z2020-12-04T00:00:00.000000Z.MSEED")
time=18400+(window/48)
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
        if k>4772:
            break
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


plt.figure(3)
plt.scatter(index,corre,s=0.2)
fig=plt.gcf()
fig.set_size_inches(18.5, 10.5)
ax = plt.gca()
ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
date_format=mdates.DateFormatter('%d,%b,%y')
ax.xaxis.set_major_formatter(date_format)
for label in ax.get_xticklabels():
    label.set_rotation(20)
    label.set_horizontalalignment('right')

plt.xlabel("Time")
ax.set_xlim([18410,18590])
plt.ylabel("cross correlation coeifficient")
plt.title("cross correlation in window length "+str(round((window/24),2))+" "+str(UTCDateTime(mdates.num2date(18410)))+" - "+str(UTCDateTime(mdates.num2date(18590))))

plt.figure(4)

plt.scatter(index,lag,s=0.2)
fig=plt.gcf()
fig.set_size_inches(18.5, 10.5)
ax = plt.gca()
ax.set_xlim([18410,18590])
ax.set_ylim([-10,10])
plt.title("cross correlation lag time in hours  "+str(UTCDateTime(mdates.num2date(18410)))+" - "+str(UTCDateTime(mdates.num2date(18590))))
ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
date_format=mdates.DateFormatter('%d,%b,%y')
plt.ylabel("Lag in hours")
ax.xaxis.set_major_formatter(date_format)
for label in ax.get_xticklabels():
    label.set_rotation(20)
    label.set_horizontalalignment('right')
# fig = plt.figure()
# fig.set_size_inches(10,6)
# ax = fig.add_axes([0,0,1,1])
# ax.bar([])
# plt.show()

# import pylab as p
# p.plot(tr2.data)
# p.plot(medfilt(tr2.data,7))
# p.scatter(tr1.data)
# p.scatter(medfilt(tr1.data,4))
# p.show ()

# st.normalize()
#tr1.filter("bandpass",freqmin=0.0000005,freqmax=2)
#tr2.filter("bandpass",freqmin=0.0000005,freqmax=2)
# dn= tr2.stats.endtime
#tr2.plot()
#tr2.plot(starttime=dt , endtime=dt+6000000)
#tr1.plot()

# tr1d=medfilt(tr1d,25)
# tr2d=medfilt(tr2d,25)

# tr2s=tr2.slice(starttime,starttime+168*3600)
# tr2s.plot()
# len(tr2d)