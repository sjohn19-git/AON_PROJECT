#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 21 22:22:39 2022
@author: sebinjohn
"""

import numpy as np
import scipy.io 
import scipy
import matplotlib.pyplot as plt 
import random
from scipy import optimize
from obspy import UTCDateTime
import matplotlib.dates as mdates

seis=np.load("/Users/sebinjohn/AON_PROJECT/Data/ML_seismic_train/final_seis_ML.npy")
sea_ice=np.load("/Users/sebinjohn/AON_PROJECT/Data/sea_ice_con/sea_ice_grdcut_5844.npy")
wave=np.load("/Users/sebinjohn/AON_PROJECT/Data/wave/wave_mlo_5844.npy")

sea_ice[sea_ice>1]=0



   
seis_ML=np.zeros((10,5844))
j=-1
for i in range(35064):
    if i%6==0:
        j+=1
        seis_ML[:,j]=seis[i,:]

seis_ML[seis_ML==0]=np.nan       
wave[np.isnan(wave)]=0
      
mean=-np.mean(wave)*np.ones(wave.shape)
div=(1/(np.amax(wave)-np.amin(wave)))*np.ones(wave.shape)

wave=wave+mean
wave=wave*div

wave=wave+-(np.amin(wave))



def feature_scaling(y):
    mean=-np.nanmean(y)*np.ones(y.shape)
    div=(1/(np.nanmax(y)-np.nanmin(y)))*np.ones(y.shape)
    y=y+mean
    y=y*div
    return y



seis_ML_f=feature_scaling(seis_ML)  

seis_ML_f=seis_ML_f+-(np.nanmin(seis_ML_f))




def X_make(sea_ice,wave,seis):
    X=np.zeros((5844,3607))
    for i in range(5844):
        X[i,:]=np.hstack((sea_ice[:,i].flatten(),wave[:,i].flatten()))
        print(i)
    return X

avail_dat=[]
for i in range(5844):
    if any(np.isnan(seis_ML_f[:,i])):
        pass
    else:
        avail_dat.append(i)



input_layer_size  = 3607
hidden_layer_size = 3607
num_labels = 1
x=X[avail_dat[:3000],:]
y=seis_ML_f[0,avail_dat[:3000]]
lam=1


def sigmoid(z):
    ex=1/((1+np.exp(-z)))
    return ex

def sigmoid_grd(z):
    ex=sigmoid(z)*(1-sigmoid(z))
    return ex

def random_init(Lin,Lout):
    epsilon_init=0.12
    ex=np.random.rand(Lout,Lin+1)*(2*epsilon_init)-epsilon_init
    return ex




def nnCostFunction(nn_params,input_layer_size, hidden_layer_size,num_labels,X, y, lam):
    theta1=nn_params[:(hidden_layer_size*(input_layer_size+1))].reshape((hidden_layer_size,input_layer_size+1),order="F")
    theta2=nn_params[(hidden_layer_size*(input_layer_size+1)):].reshape((num_labels,hidden_layer_size+1),order="F")
    theta1_grad = np.zeros(theta1.shape)
    theta2_grad = np.zeros(theta2.shape)
    m=X.shape[0]
    a1=X
    le=(np.shape(a1)[0])
    a1=np.hstack((np.ones((le, 1)),a1 ))
    z2=a1@theta1.T
    a2=sigmoid(z2)
    le=(np.shape(a2)[0])
    a2=np.hstack((np.ones((le, 1)),a2 ))
    z3=a2@theta2.T
    a3=sigmoid(z3)
    J=0
    del3=np.zeros((num_labels,m))
    for i in range(m):
        #yi=y[:,i]
        yi=y[i]
        Ji=np.abs(np.square(yi)-np.square(a3))
        J+=1/m*np.sum(Ji)
        del3[:,i]=(a3[i].reshape(num_labels,1)-yi.reshape((num_labels,1))).flatten()
    J=J+(lam/(2*m))*(np.sum(np.square(theta1[:,1:]))+np.sum(np.square(theta2[:,1:])))
    z2=np.hstack((np.zeros((le, 1)),z2))
    del2=(theta2.T@del3)*sigmoid_grd(z2).T
    Del2=del3@a2
    Del1=del2[1:,:]@a1
    theta1_grad=(1/m)*(theta1_grad+Del1)
    reg1=(lam/m)*(theta1_grad)
    reg1[:,0]=0
    theta1_grad=theta1_grad+reg1
    theta2_grad=(1/m)*(theta2_grad+Del2)
    reg2=(lam/m)*(theta2_grad)
    reg2[:,0]=0
    theta2_grad=theta2_grad+reg2
    grd=np.concatenate((theta1_grad.flatten(order="F"),theta2_grad.flatten(order="F")),axis=0)  
    return J,grd



J,grd=nnCostFunction(initial_nn_params,input_layer_size, hidden_layer_size,num_labels,x, y, lam)

initial_Theta1 = random_init(input_layer_size, hidden_layer_size);
initial_Theta2 = random_init(hidden_layer_size, num_labels);

initial_nn_params =np.concatenate((initial_Theta1.flatten(order="F"),initial_Theta2 .flatten(order="F")),axis=0) 



def optimization(initial_nn_params,input_layer_size, hidden_layer_size,num_labels,X, y,lam):
    args=(input_layer_size, hidden_layer_size,num_labels,X, y, lam)
    res=scipy.optimize.minimize(nnCostFunction,initial_nn_params,args=args,method="Newton-CG",jac=True,options={'disp': True,"maxiter":50})
    theta1=res.x[:(hidden_layer_size*(input_layer_size+1))].reshape((hidden_layer_size,input_layer_size+1),order="F")
    theta2=res.x[(hidden_layer_size*(input_layer_size+1)):].reshape((num_labels,hidden_layer_size+1),order="F")
    return theta1,theta2

theta1o,theta2o=optimization(initial_nn_params,input_layer_size, hidden_layer_size,num_labels,x, y,lam)

xt=X[avail_dat[:],:]
yt=seis_ML_f[0,avail_dat[:]]



def predict(theta1o,theta2o,xt):
    m=xt.shape[0]
    num_labels = theta2o.shape[0]
    a1=xt
    le=(np.shape(a1)[0])
    a1=np.hstack((np.ones((le, 1)),a1 ))
    h1=sigmoid(a1@theta1o.T)
    h1s=np.hstack((np.ones((le, 1)),h1 ))
    h2=sigmoid(h1s@theta2o.T)
    return h2

h2=predict(theta1o,theta2o,xt)


np.save("theta1_sta0",theta1o)
np.save("theta2_sta0",theta2o)

st=UTCDateTime(2018,1,1)
et=UTCDateTime(2022,1,1)
windw=1*6
time_frame=np.array([])
for i in range (int((et-st)/(3600*windw))):
    time_frame=np.append(time_frame,st)
    st=st+(3600*windw)

avail_time_f=time_frame[avail_dat[:]]
plot_time=[]
for ele in avail_time_f:
    plot_time.append(ele.matplotlib_date)
    

date_format=mdates.DateFormatter('%d-%b-%Y')

fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(18,6))
axes.scatter(plot_time[:],yt[:],s=1.5,c="blue")
axes.scatter(plot_time[:],h2[:],s=1.5,c="orange")
axes.set_xlim([min(plot_time[:]),max(plot_time[:])])
axes.xaxis.set_major_formatter(date_format)
axes.set_ylim([0,0.6])
axes.vlines(plot_time[3000],0,1)
fig.suptitle("A19K")
axes.set_ylabel("power_normalized")
