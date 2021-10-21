import obspy
from obspy import read 
import numpy as np
import obspy.signal.cross_correlation as cr


st=read("/home/sjohn/AON_PROJECT/O14K/O14K_2020-05-18T00:00:00.000000Z2020-12-04T00:00:00.000000Z.MSEED")
tr1=st[0]
tr2=st[1]
tr1=tr1.data
tr2=tr.data
len(tr1)
j=0
corre=np.array([])
for i in range(len(tr1)):
    k=j+24
    tr1par=tr1[j:k]
    tr2par=tr2[j:k]
    j=k
    coeif=cr.xcorr(tr1par, tr2par,1)
    corre=np.append(corre,np.array(coeif[1]))
    
tr1
