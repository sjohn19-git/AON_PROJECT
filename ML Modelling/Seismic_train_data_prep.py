#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 17:04:55 2022

@author: sebinjohn
"""
os.chdir("/Users/sebinjohn/AON_PROJECT/Data/codes")
import matplotlib.pyplot as plt
import os
import numpy as np
import glob
import pandas as pd
from obspy import UTCDateTime
from medfilt import medfilt


st=UTCDateTime(2018,1,1)
et=UTCDateTime(2022,1,1)
windw=1

time_frame=[]
for i in range (int((et-st)/(3600*windw))):
    time_frame.append(st)
    st=st+(3600*windw)




os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
freq=[]
name= pd.read_xml("pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
for i in range (95):
    freq.append(name.iloc[i]['freq'])
freq.append(19.740300000000000)
len(freq)

final_seis_ML=np.zeros((len(time_frame),10))

inde=-1

for files in glob.glob("/Users/sebinjohn/AON_PROJECT/Data/ML_seismic_train/*/*.npy"):
    print(files)
    inde+=1
    seis=np.load(files)
    res=seis[53:61,:]
    out_mean=mean(res)
    out_mean_med=medfilt(out_mean,7)
    interp_inde=np.array([])
    interp_x=out_mean_med.copy()
    datagap=np.where(interp_x==0)[0]
    data_x=np.where(interp_x!=0)[0]
    for i in range(len(data_x)-1):
        if data_x[i+1]-data_x[i]<12 and data_x[i+1]-data_x[i]>1:
            interp_inde=np.append(interp_inde,np.array( [i for i in range(int(data_x[i])+1,int(data_x[i+1]))]))
        else:
            continue
    if len(interp_inde)>1:
        interp=np.interp(interp_inde, data_x.reshape(np.shape(data_x)[0]), interp_x[data_x].reshape(np.shape(data_x)[0]))
        interp_inde=(interp_inde+1).astype("int32")
        interp_x[interp_inde]=interp
    else:
        pass
    final_seis_ML[:,inde]=interp_x







def mean(res):
    mea=np.mean(res,axis=0)
    return mea

os.chdir("/Users/sebinjohn/AON_PROJECT/Data/ML_seismic_train")

np.save("final_seis_ML.npy", final_seis_ML)
