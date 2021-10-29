import datetime as dt
import matplotlib.dates as mdates
from obspy import UTCDateTime
from obspy import read
import matplotlib.pyplot as plt
import dill
from datetime import datetime as datetime
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
import dill
import os
import re
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
ans=UTCDateTime(tmin)-UTCDateTime(startt)
startcol=int(ans/(60*60))
endcol=int(startcol+(UTCDateTime(tmax)-UTCDateTime(tmin))/(60*60))
print(startcol,endcol)
powe=final[j,startcol:endcol]
fint_frames=fint_frames[startcol:endcol]
startcolw=(np.where(finalw==xmin))
endcolw=(np.where(finalw==xmax))
wspd=finalw[(int(startcolw[0])):(int(endcolw[0]))]
startcolp=(np.where(finalp==xmin))
endcolp=(np.where(finalp==xmax))
prsr=finalp[(int(startcolp[0])):(int(endcolp[0]))]
len(prsr)

tmin
tmax
#----------------------------------------------------------------------------

net1="TIMESERIES "+net[-1]+"_"
sta1=sta+"Power"+"__"
cha1=cha+"_M,"
samp1=" "+str(len(powe))+" samples, "
sampr1="0.0002777777777777 sps, "
startt1=str(UTCDateTime(tmin))[:-1]+",  TSPAIR, FLOAT, COUNTS\n"
header=net1+sta1+cha1+samp1+sampr1+startt1
startt1
header
with open (net[-1]+"_"+sta+"_"+cha+"_"+str(UTCDateTime(tmin))[:-1]+"power.txt",'w') as f:
    f.write(header)
    for i in range(len(powe)):
        ti=str(UTCDateTime(mdates.num2date(fint_frames[i])))[:-1]
        f.write(ti+" "+str(powe[i])+"\n")
st = read("/home/sjohn/AON_PROJECT/O14K/AK_O14K_BHZ_2020-05-18T00:00:00.000000power.txt")
#--------------------------------------------------------------------------------------
sta1=sta+"wnd"+"__"
finalw
header=net1+sta1+cha1+samp1+sampr1+startt1
tminw=mdates.num2date(wspd[0,0])
tmaxw=mdates.num2date(wspd[-1,0])
samp1=" "+str(len(wspd))+" samples, "
namew=net[-1]+"_"+sta+"_"+cha+"_"+str(UTCDateTime(tminw))[:-1]+"wind.txt"
with open (net[-1]+"_"+sta+"_"+cha+"_"+str(UTCDateTime(tminw))[:-1]+"wind.txt",'w') as f:
    f.write(header)
    for i in range(len(wspd)):
        ti=str(UTCDateTime(mdates.num2date(wspd[i,0])))[:-1]
        f.write(ti+" "+str(wspd[i,1])+"\n")
st+=read(os.path.join("/home/sjohn/AON_PROJECT/O14K/",namew))


#----------------------------------------------------------------------------------------
sta1=sta+"prsr"+"__"
samp1=" "+str(len(powe))+" samples, "
sampr="0.0002777777777777 sps, "
cha1="LDV_M,"
startt1=str(UTCDateTime(tmin))[:-1]+",  SLIST, FLOAT, COUNTS\n"
header=net1+sta1+cha1+samp1+sampr+startt1
print(header)
tminp=mdates.num2date(prsr[0,0])
tmaxp=mdates.num2date(prsr[-1,0])
namep=net[-1]+"_"+sta+"_"+cha+"_"+str(UTCDateTime(tminp))[:-1]+"pressure.txt"
with open (net[-1]+"_"+sta+"_"+cha+"_"+str(UTCDateTime(tminp))[:-1]+"pressure.txt",'w') as f:
    f.write(header)
    for i in range(len(prsr)):
        ti=str(UTCDateTime(mdates.num2date(prsr[i,0])))[:-1]
        f.write(ti+" "+str(prsr[i,1])+"\n")

plines=[]
with open (os.path.join("/home/sjohn/AON_PROJECT/O14K/",namep),'r') as infile:
    for line in infile:
        plines.append(line)

#------------------------------------------------------------------------------------------
startt=UTCDateTime(mdates.num2date(xmin))
endt=UTCDateTime(mdates.num2date(xmax))
pressr=np.zeros(len(st[0]))
for i in range(len(st[0])):
    pattern=re.compile("(%s) ([\d\.]+)"%str(startt)[:-1])
    for ele in plines:   
        matches=re.finditer(pattern, ele)
        for match in matches:
            pressr[i]=match.group(2)
 
    startt+=3600

with open (net[-1]+"_"+sta+"_"+cha+"_"+str(UTCDateTime(tminp))[:-1]+"pressuref.txt",'w') as f:
    f.write(header)
    for i in range(len(pressr)):
        f.write(str(pressr[i])+"\n")
 
os.remove(os.path.join("/home/sjohn/AON_PROJECT/O14K/",namep)) 
k=net[-1]+"_"+sta+"_"+cha+"_"+str(UTCDateTime(tminp))[:-1]+"pressuref.txt"
k
st+=read("/home/sjohn/AON_PROJECT/O14K/"+k)
#----------------------------------------------------------------------------------------
st
st.write("O14K_"+str(UTCDateTime(tmin))+str(UTCDateTime(tmax))+".MSEED", format="MSEED")  


