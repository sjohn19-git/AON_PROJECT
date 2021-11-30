import os
import glob
from os.path import isfile, join 
import dill
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

def unit_sort(ele):
    return float(ele.split("/")[-1][3:].split(".")[-2])
def download():
    stepta=(starttimeta-endtimeta)/76
    stepak=(starttimeak-endtimeak)/26
    timesta=np.arange(starttimeta,endtimeta,-stepta)
    timesak=np.arange(starttimeak,endtimeak,-stepak)
    starttimes=[]
    endtimes=[]
    net=[]
    for i in range(75):
        k=str(timesta[i])
        starttimes.append((k[:10])+"T00:00:00")
    starttimes.append((str(timesta[-1])[:10])+"T00:00:00")
    for i in range(25):
        starttimes.append((str(timesak[i])[:10])+"T00:00:00")
    for i in range(1,76):
        endtimes.append((str(timesta[i])[:10])+"T00:00:00")
    endtimes.append((str(endtimeta)[:10])+"T00:00:00")
    for i in range(1,26):
        endtimes.append((str(timesak[i])[:10])+"T00:00:00")
    starttimes.append(endtimes[-1])
    endtimes.append((str(endtimeak)[:10])+"T00:00:00")
    for i in range(76):
        net.append("TA")
    for i in range(26):
        net.append("AK")
    for z in range (len(starttimes)):
        url1='https://service.iris.edu/mustang/noise-psd/1/query?target='+net[z]+"."+sta+"."+loc+"."+cha+".M&starttime="+str(starttimes[z])+"&endtime="+str(endtimes[z])+"&format=xml"
        url2='https://service.iris.edu/irisws/timeseries/1/query?net='+net[z]+"&sta="+sta+"&cha=LDV&start="+str(starttimes[z])+"&end="+str(endtimes[z])+"&deci=0.0002777777777777&format=ascii2&loc=EP"
        try:
            print("urlpdf",url1)
            resp1= requests.get(url1)
            print("urlpressure",url2)
            resp2=requests.get(url2)
            with open('pdf'+str(z)+".xml", 'w') as foutput:
                foutput.write(resp1.content.decode('utf-8'))
            with open('pressure'+str(z)+".txt", 'w') as foutput:
                foutput.write(resp2.content.decode('utf-8'))
        except:
            print("check")
            break

pattern=re.compile("<Psd target=\"([A-z]*)\.[\.A-z\d\"\s]+start=\"(\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d\.\d\d\dZ)\" end=\"(\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d\.\d\d\dZ)\"")
pattern1=re.compile("<value freq=\"([1]?[0-9]?\.[\d]+)\"[\s,a-z,=]*\"(\-?[\d]+\.?[\d]+)\"")

folders=[]
for path in glob.glob("/home/sjohn/spectrogram/P16K"):
    folders.append(path)


for i in range(3,len(folders)):
    path=folders[i]
    print(path)
    os.chdir(path)
    try:
        dill.load_session(glob.glob(join(os.getcwd()+'/*v_*.db'))[0])
    except:
         dill.load_session(glob.glob(join(os.getcwd()+'/*.db'))[0])
    download()
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
    print(path)
                
for i in range(len(folders)):
    sta=folders[i].split("/")[-1]
    shutil.copy(join(folders[i],"final_"+sta+".npy"),join("/tmp/.x2go-sjohn/media/disk/_Users_sebinjohn_AON_PROJECT",sta))
    shutil.copy(join(folders[i],"final_"+sta+".npy"),join('/home/sjohn/Data',sta))
    print(sta)

    
  
