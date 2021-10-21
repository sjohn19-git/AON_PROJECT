import obspy
from obspy import read 

st=read("/home/sjohn/AON_PROJECT/O14K/O14K_2020-05-18T00:00:00.000000Z2020-12-04T00:00:00.000000Z.MSEED")
st
st[0].plot()
