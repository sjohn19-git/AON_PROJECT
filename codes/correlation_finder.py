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
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes

#--------------------------------------------------------------------------------
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
#---------------------------------------------------------------------
j=53
xmin=18400
xmax=18600

fig=plt.figure()

fig.set_figheight(10)
fig.set_figwidth(16)

host = fig.add_axes([2, 0.1, 0.7, 0.3], axes_class=HostAxes)
host.axis["right"].set_visible(False)
par1 = ParasiteAxes(host)
host.parasites.append(par1)

par1.axis["right"].set_visible(True)
par1.axis["right"].major_ticklabels.set_visible(True)
par1.axis["right"].label.set_visible(True)

par1.axis["bottom2"] = par1.new_fixed_axis(loc="bottom", offset=(0, -25))
par1.axis["bottom2"].set_visible(False)
host.set_xlim(xmin,xmax)
host.set_ylim(-80,-160)
par1.set_xlim(xmin,xmax)
par1.set_ylim([0,18])
par1.set_ylabel("Wind (m/s)")
host.set_ylabel("power ("+str(freq[j])+"Hz)")
host.title.set_text(sta+"_"+cha+" "+str(UTCDateTime((mdates.num2date(xmin)).strftime('%Y-%m-%d')))+" - "+str(UTCDateTime((mdates.num2date(xmax)).strftime('%Y-%m-%d'))))
x=(fint_frames[1])
y=(final[45,1])
p1 = host.scatter(fint_frames[:],finalp[j,:],s=0.5,label="power (db)",c="red")
p2 =par1.scatter(finalw[:,0],finalw[:,1], s=0.5,label="power (db)",c="green")
host.xaxis.set_major_formatter(date_format)
host.xaxis.set_minor_formatter(date_format1)
host.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(SU)))
plt.setp(host.axis["bottom"].major_ticklabels, rotation=40,ha="right")
host.legend()
fig.savefig('pressure&powerO14K.eps', bbox_inches='tight')