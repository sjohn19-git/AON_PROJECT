import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
import pickle
import matplotlib.dates as mdates


st=UTCDateTime(2015,1,1)
et=UTCDateTime(2021,6,1)
windw=24
tim_len=int((et-st)/(3600*windw))
time_frame=[]
for i in range (int((et-st)/(3600*windw))):
    time_frame.append(st)
    st=st+(3600*windw)


with open(("/Users/sebinjohn/AON_PROJECT/Data/Timeseris_PSD/timeseries_2015-2021-1day.pkl"),"rb") as f:
    mean_timeseries,median_timeseries = pickle.load(f) 
    
mat_t=[]
for ele in time_frame:
    mat_t.append(ele.matplotlib_date)
    
    
median_timeseries[median_timeseries==-180]=np.nan
date_format=mdates.DateFormatter('%d,%b,%y') 
 
f,axs=plt.subplots(5,1,sharex=True,figsize=(20,24))

axs[0].plot(mat_t,median_timeseries[0,23,1:])
axs[0].set_title("Power (PSD) in ice band-2015/1/1 - 2021/06/1\n\nC26K")
axs[2].plot(mat_t,median_timeseries[0,27,1:])
axs[2].set_title("B18K")
axs[1].plot(mat_t,median_timeseries[0,30,1:])
axs[1].set_title("A21K")
axs[3].plot(mat_t,median_timeseries[0,25,1:])
axs[3].set_title("C16K")
axs[4].plot(mat_t,median_timeseries[0,19,1:])
axs[4].set_title("F15K")
axs[4].set_xlim([min(mat_t),max(mat_t)])
axs[4].xaxis.set_major_locator(mdates.MonthLocator(bymonth=[3,6,9,12]))
axs[4].xaxis.set_major_formatter(date_format)
for label in axs[4].get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')