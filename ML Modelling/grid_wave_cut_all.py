#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 02:03:45 2022

@author: sebinjohn
"""

import numpy as np
import matplotlib.pyplot as plt
from global_land_mask import globe
import pygmt
from scipy.interpolate import RegularGridInterpolator
import os

seis=np.load("/Users/sebinjohn/AON_PROJECT/Data/ML_seismic_train/final_seis_ML.npy")
wave=np.load("/Users/sebinjohn/AON_PROJECT/Data/wave/wav_arr.npy")
   



wave_ML=np.zeros((86,220,5844))
j=-1
for i in range(35064):
    if i%6==0:
        j+=1
        wave_ML[:,:,j]=wave[:,:,i]
 


lat_wave_ML=np.arange(73,30,-0.5)
lon_wave_ML=np.arange(150,180,0.5)
lon_wave_ML=np.append(lon_wave_ML,np.arange(-180,-100,0.5))

for ele in lon_wave_ML:
    if ele<0:
        lon_wave_ML[np.where(lon_wave_ML==ele)]=ele+360

sort_lat=np.sort(lat_wave_ML)


lat_wave_rs=np.arange(73,30,-1.5)
lon_wave_rs=np.arange(150,260,1.5)
X,Y=np.meshgrid(lat_wave_rs,lon_wave_rs)
wave_ml=np.zeros((74,29,5844))
for i in range(5844):
    my_interpolating_function = RegularGridInterpolator((sort_lat, lon_wave_ML),wave_ML[:,:,i])
    c=my_interpolating_function((X,Y))
    c=np.fliplr(c)
    wave_ml[:,:,i]=c
    




plt.pcolormesh(Y,X,wave_ml[:,:,10])
plt.pcolormesh(lon_wave_ML,lat_wave_ML,wave_ML[:,:,0])



def grid_downgrade(X,Y):
    points=[]
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
                points.append((X[i,j],Y[i,j]))
    indices=[]
    for point in points:
        if point[1]<180:
            if globe.is_ocean(point[0],point[1]):
                indices.append((np.where(lat_wave_rs==point[0])[0][0],np.where(lon_wave_rs==point[1])[0][0]))
        else:
            if globe.is_ocean(point[0],point[1]-360):
                indices.append((np.where(lat_wave_rs==point[0])[0][0],np.where(lon_wave_rs==point[1])[0][0]))
    return indices
 
def Reverse(tuples):
    new_tup = ()
    for k in reversed(tuples):
        new_tup = new_tup + (k,)
    return new_tup

    
wave_mlo=np.zeros((1339,5844))
for i in range(1339):
    wave_mlo[i,:]=wave_ml[Reverse(indices[i])]

    
jx=[]
jy=[]

for ele in indices:
    jx.append(lon_wave_rs[ele[1]])
    jy.append(lat_wave_rs[ele[0]])


os.chdir("/Users/sebinjohn/AON_PROJECT/Data/wave")
np.save("wave_mlo_5844.npy",wave_mlo)

fig=pygmt.Figure()
fig.coast(region=[150,260,30, 73], projection="L-159/35/33/45/22c",shorelines=False,frame="a",land="200/200/200",borders=["1/0.5p,black", "2/0.5p,red"])
fig.plot(x=jx,y=jy,style="i0.1c", color="black")
fig.show(method="external")