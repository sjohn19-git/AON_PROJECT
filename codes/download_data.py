import os
import glob
from os.path import isfile, join 
import datetime as dt
import matplotlib.dates as mdates
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import matplotlib.cm as cm
import matplotlib.colors as cl
from matplotlib import gridspec
import shutil
from datetime import datetime
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import requests
import pickle


def pickle_maker():
    sta=input("station?")
    try:
        os.mkdir(join("/home/sjohn/spectrogram/"+sta))
    except:
        pass
    starttimeta=UTCDateTime(input("start time TA?"))
    endtimeta=UTCDateTime(input("endtime time TA?"))
    starttimeak=UTCDateTime(input("start time AK?"))
    endtimeak=UTCDateTime(input("endtime time AK?"))
    loc=input("location?")
    cha=input("channel?")
    stepta=(starttimeta-endtimeta)/76
    stepak=(starttimeak-endtimeak)/100
    timesta=np.arange(starttimeta,endtimeta,-stepta)
    timesak=np.arange(starttimeak,endtimeak,-stepak)
    starttimes=[]
    endtimes=[]
    net=[]
    for i in range(75):
        k=str(timesta[i])
        starttimes.append((k[:10])+"T00:00:00")
    starttimes.append((str(timesta[-1])[:10])+"T00:00:00")
    for i in range(99):
        starttimes.append((str(timesak[i])[:10])+"T00:00:00")
    for i in range(1,76):
        endtimes.append((str(timesta[i])[:10])+"T00:00:00")
    endtimes.append((str(endtimeta)[:10])+"T00:00:00")
    for i in range(1,100):
        endtimes.append((str(timesak[i])[:10])+"T00:00:00")
    starttimes.append(endtimes[-1])
    endtimes.append((str(endtimeak)[:10])+"T00:00:00")
    for i in range(76):
        net.append("TA")
    for i in range(100):
        net.append("AK")
    return sta, starttimeta, endtimeta,starttimeak,endtimeak,cha,loc,net,starttimes,endtimes




def unit_sort(ele):
    return float(ele.split("/")[-1][3:].split(".")[-2])

def download(net,starttimes,endtimes,sta):
    for z in range (len(starttimes)):
        os.chdir(join("/home/sjohn/spectrogram/"+sta))
        url1='https://service.iris.edu/mustang/noise-psd/1/query?target='+net[z]+"."+sta+"."+loc+"."+cha+".M&starttime="+str(starttimes[z])+"&endtime="+str(endtimes[z])+"&format=xml"
        url2='https://service.iris.edu/irisws/timeseries/1/query?net='+net[z]+"&sta="+sta+"&cha=LDV&start="+str(starttimes[z])+"&end="+str(endtimes[z])+"&deci=0.0002777777777777&format=ascii2&loc=EP"
        url3='https://service.iris.edu/irisws/timeseries/1/query?net='+net[z]+"&sta="+sta+"&cha=LWS&start="+str(starttimes[z])+"&end="+str(endtimes[z])+"&deci=0.0002777777777777&format=ascii2&loc=EP"
        try:
            print("urlpdf",url1)
            resp1= requests.get(url1)
            print("urlpressure",url2)
            resp2=requests.get(url2)
            print("urlwind",url3)
            resp3=requests.get(url3)
            with open('pdf'+str(z)+".xml", 'w') as foutput:
                foutput.write(resp1.content.decode('utf-8'))
            with open('pressure'+str(z)+".txt", 'w') as foutput:
                foutput.write(resp2.content.decode('utf-8'))
            with open('wind'+str(z)+".txt", 'w') as foutput:
                foutput.write(resp3.content.decode('utf-8'))
        except:
            print("check")
            break
    


sta, starttimeta, endtimeta,starttimeak,endtimeak,cha,loc,net,starttimes,endtimes=pickle_maker()

os.chdir(join("/home/sjohn/spectrogram/"+sta))

with open((str(sta)+".pkl"),"wb") as f:
    pickle.dump([sta, starttimeta, endtimeta,starttimeak,endtimeak,sta,cha,loc],f)


download(net,starttimes,endtimes,sta)


pattern=re.compile("<Psd target=\"([A-z]*)\.[\.A-z\d\"\s]+start=\"(\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d\.\d\d\dZ)\" end=\"(\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d\.\d\d\dZ)\"")
pattern1=re.compile("<value freq=\"([1]?[0-9]?\.[\d]+)\"[\s,a-z,=]*\"(\-?[\d]+\.?[\d]+)\"")



#-------------------------POWER----------------------------------------------

filespdf=[]
for name in glob.glob(os.path.join(os.getcwd(),"pdf*.xml")):
    filespdf.append(name)
unit_sort(filespdf[0])
filespdf = sorted(filespdf, key=unit_sort)
lines=[]
for names in filespdf:
    with open(names,'r') as infile:
        for line in infile:
            lines.append(line)

time_frames=np.arange(starttimeta,endtimeak,1800)
for i in range(len(time_frames)):
    time_frames[i]=mdates.date2num(time_frames[i])
count=0
avail_time=[]
for line in lines:
    result=re.finditer(pattern,line)
    for match in result:
        count+=1
        avail_time.append(mdates.date2num(UTCDateTime(match.group(2))))  
both = set(avail_time).intersection(time_frames)
result=[]
for ele in both:
    result.append(np.where(time_frames==ele))
result.sort()
ma=np.zeros((96,len(time_frames)))
r=-1
co=0
nmatch=0
for line in lines:
    result1=re.finditer(pattern1,line)
    for match in result1:
        nmatch+=1
        try:
            if r==95:
                r=0
                co+=1
            else:
                r+=1 
            ma[r,result[co][0][0]]=(match.group(2))
        except:
            pass
ma=np.delete(ma, range(1, ma.shape[1], 2), axis=1)
np.save("final_"+sta+'.npy',ma)
print(ma)


#--------------------------------pressure--------------------------------------
          
files=[]
for name in glob.glob(os.path.join(os.getcwd(),"pressure*.txt")):
    files.append(name)
plines=[]
for names in files:
    with open(names,'r') as infile:
        for line in infile:
            plines.append(line)
len(plines)

pattern=re.compile("(\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d.\d\d\d\d\d\d)  ([\d\.]+)")
r=-1
nmat=0
for line in plines:
    result2=re.finditer(pattern,line)
    for match in result2:
        nmat+=1
pma=np.zeros((nmat,2))
for line in plines:
    result2=re.finditer(pattern,line)
    for match in result2:
        r+=1
        #pma[r,1]=(match.group(1))
        pma[r,0]=mdates.date2num(UTCDateTime(match.group(1)))
        pma[r,1]=((float(match.group(2))/10132.5))
np.save('finalp_'+sta+'.npy',pma)
print(pma)


#-----------------------------------------------------wind---------------------------

filesw=[]
for name in glob.glob(os.path.join(os.getcwd(),"wind*.txt")):
    filesw.append(name)
wlines=[]
for names in files:
    with open(names,'r') as infile:
        for line in infile:
            wlines.append(line)


pattern=re.compile("(\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d.\d\d\d\d\d\d)  ([\d\.]+)")

r=-1
nmat=0
for line in wlines:
    result2=re.finditer(pattern,line)
    for match in result2:
        nmat+=1
wma=np.zeros((nmat,2))
for line in wlines:
    result2=re.finditer(pattern,line)
    for match in result2:
        r+=1
        #pma[r,1]=(match.group(1))
        wma[r,0]=mdates.date2num(UTCDateTime(match.group(1)))
        wma[r,1]=((float(match.group(2))))
np.save('finalw_'+sta+'.npy',wma)
print(wma)


tocopy=["final_","finalp_","finalw_"]

os.mkdir(join("/tmp/.x2go-sjohn/media/disk/_Users_sebinjohn_AON_PROJECT/"+sta))
os.mkdir("/home/sjohn/Data/"+sta)
for i in range(len(tocopy)):
    shutil.copy(join(os.getcwd(),tocopy[i]+sta+".npy"),join("/tmp/.x2go-sjohn/media/disk/_Users_sebinjohn_AON_PROJECT",sta))
    shutil.copy(join(os.getcwd(),tocopy[i]+sta+".npy"),join('/home/sjohn/Data',sta))

shutil.copy(join(os.getcwd(),str(sta)+".pkl"),join("/tmp/.x2go-sjohn/media/disk/_Users_sebinjohn_AON_PROJECT",sta))
shutil.copy(join(os.getcwd(),str(sta)+".pkl"),join('/home/sjohn/Data',sta))


