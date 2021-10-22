import obspy
from obspy import read 
import numpy as np
import obspy.signal.cross_correlation as cr
import matplotlib.pyplot as plt

shift=1
window=480
st=read("/home/sjohn/AON_PROJECT/O14K/O14K_2020-05-18T00:00:00.000000Z2020-12-04T00:00:00.000000Z.MSEED")

tr1=st[0]
tr2=st[2]
tr2.plot()
tr1=tr1.data
tr2=tr2.data
len(tr1)
j=0
k=window
corre=np.array([])
index=np.array([])
for i in range((int(len(tr1)/shift))):
#    print(j,k)
    tr1par=tr1[j:k]
    tr2par=tr2[j:k]
    coeif=cr.xcorr(tr1par, tr2par,1)
    corre=np.append(corre,np.array(coeif[1]))
    index=np.append(index,np.array(coeif[0]))
    k+=shift
    j+=shift
    if k>4800:
        break



val=[]



for i in range(len(corre)):
    val.append(format(corre[i]))


cun,cuo,cup,cuq=0,0,0,0
for ele in val:
    if float(ele)<=-0.5:
        cun+=1
    elif float(ele)<=0 and float(ele)>-0.5:
        cuo+=1
    elif float(ele)>0 and float(ele)<=0.5:
        cup+=1
    elif float(ele)>0.5 and float(ele)<=1:
        cuq+=1
    else:
        print(ele)
x=len(corre)
cun
cuo
cup
cuq

fig = plt.figure()
fig.set_size_inches(10,6)
ax = fig.add_axes([0,0,1,1])
ax.bar(["xcorr<-0.5",'-0.5<xcorr<0',"0<xcorr<0.5","xcorr>0.5)"],[cun,cuo,cup,cuq])
ax.set_ylabel("counts")
ax.set_title('correlation between pressure and power of frequency (0.5187 Hz), window length='+str(window)+"Hour")
plt.show()
