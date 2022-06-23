#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 21:17:52 2022

@author: sebinjohn
"""

import os 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

st=UTCDateTime(2018,1,1)
et=UTCDateTime(2022,1,1)
windw=1
time_frame=[]
for i in range (int((et-st)/(3600*windw))):
    time_frame.append(st)
    st=st+(3600*windw)

os.chdir("/Users/sebinjohn/AON_PROJECT/Data/ML_seismic_train/C27K")

seismic_data=np.load("final_MLC27K.npy")

res=seismic_data[53:61,:]
out_mean=np.mean(res,axis=0)
seismic_dat_h=np.copy(out_mean).reshape((1461,24))
seismic_dat_d=np.zeros((seismic_dat_h.shape[0]))

for i in range(seismic_dat_h.shape[0]-1400):
    k=seismic_dat_h[i,:]
    sumk=np.sum(k)
    lent=np.count_nonzero(k)
    seismic_dat_d[i]=sumk/lent
    print(sumk)
    print(lent)

mattime=[]
for i in range(len(time_frame)):
    if i%24==0:
        mattime.append(time_frame[i].matplotlib_date)





fig,ax=plt.subplots(figsize=(12,3))
date_format=mdates.DateFormatter('%d,%b,%y') 
ax.plot(mattime,seismic_dat_d)
ax.set_xlim(17532+600,17532+610)
ax.set_ylim(-160,-120)
ax.xaxis.set_major_formatter(date_format)