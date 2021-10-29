import math
import obspy 
from obspy import read 
import numpy as np
import obspy.signal.cross_correlation as cr
import matplotlib.pyplot as plt

x=np.linspace(0,5*math.pi,1000)
x
    
t=np.array([])
t1=np.array([])

for ele in x:
    t=np.append(t,math.sin(ele))
    t1=np.append(t1,2*math.exp(-ele/10)*math.sin(ele+(math.pi)/4))
    coeif=cr.xcorr(t, t1,0)
coeif
plt.plot(t)
plt.plot(t1)
