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
import pygmt
import pickle





def plotter(stations):
    os.chdir("/Users/sebinjohn/AON_PROJECT/Data/Video Making")
    with open("metadta.pkl","rb") as f:
       long,lat,stationo ,env = pickle.load(f)
    for ele in stations:
        sta=ele
        os.chdir("/Users/sebinjohn/AON_PROJECT/Map_with_spectrogram")
        if not os.path.exists(sta+".png"):
            os.chdir(r"/Users/sebinjohn/AON_PROJECT/Data/"+sta)
            with open((str(sta)+".pkl"),"rb") as f:
                sta, starttimeta, endtimeta,starttimeak,endtimeak,sta,cha,loc=pickle.load(f)
            
            fint_frames=np.arange(starttimeta,endtimeak,3600)
            for i in range(len(fint_frames)):
                fint_frames[i]=mdates.date2num(fint_frames[i])
                
            with open("final_"+sta+".npy", 'rb') as g:
                final=np.load(g)
            
            xmin=17532
            xmax=18682
            xlmi=xmin
            xlma=xmax
            
            freq=[]
            name= pd.read_xml("/Users/sebinjohn/AON_PROJECT/Data/Video Making/pdf0.xml", xpath="/PsdRoot/Psds[1]/Psd[1]/value[@freq]")
            for i in range (95):
                freq.append(name.iloc[i]['freq'])
            freq.append(19.740300000000000)
            len(freq)
            
            date_format=mdates.DateFormatter('%d-%B-%Y')
            date_format1=mdates.DateFormatter('%m/%d/%yT%H')
            date_format
            
            cmap = plt.get_cmap('nipy_spectral').copy()
            cmap.set_over(color = 'w')
            fig = plt.figure(constrained_layout=True)
            fig = plt.figure(facecolor=None)
            fig.patch.set_alpha(0)
            fig.set_figheight(6)
            fig.set_figwidth(18)
            spec = gridspec.GridSpec(ncols=2, nrows=1,width_ratios=[10, 1], wspace=-0.1,hspace=0.2, height_ratios=[1])
            ax1 = fig.add_subplot(spec[0])
            c = ax1.pcolormesh(fint_frames,freq,final,cmap=cmap,vmin=-195,vmax=-95,shading='auto',rasterized=True)
            # ax1.axhline(y=0.1,c="r",linestyle="--")
            # ax1.axhline(y=0.2,c="r",linestyle="--")
            #ax1.axhline(y=0.1,c="k",linestyle="--")
            #ax1.axhline(y=0.07,c="k",linestyle="--")
            nom=cl.Normalize(vmin=-195,vmax=-95)
            #ax1.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(MO, TU, WE, TH, FR, SA, SU)))
            ax1.xaxis.set_major_locator(mdates.MonthLocator(1))
            #ax1.xaxis.set_major_locator(mdates.HourLocator(interval=4))
            ax1.xaxis.set_major_formatter(date_format)
            ax1.xaxis.set_minor_formatter(date_format1)
            #ax1.set_ylim(freq[29],0.1)
            #ax1.set_ylabel("Frequency(Hz)")
            ax1.set_xlim([xmin,xmax])
            ax1.set_facecolor(color=None)
            # for label in ax1.get_xticklabels():
            #     label.set_rotation(40)
            #     label.set_horizontalalignment('right')
            ax1.set_yscale('log')
            #ax1.title.set_text(sta+"_"+cha+" "+str(UTCDateTime((mdates.num2date(xmin)).strftime('%Y-%m-%d')))[0:10]+" - "+str(UTCDateTime((mdates.num2date(xmax)).strftime('%Y-%m-%d')))[0:10])
            ax1.title.set_text(sta+"_"+cha)
            ax0 = fig.add_subplot(spec[1])
            ax0.axis('off')
            ax0.set_facecolor(color=None)
            fig.colorbar(cm.ScalarMappable(norm=nom, cmap="nipy_spectral"), ax=ax0)
            os.chdir("/Users/sebinjohn/AON_PROJECT/Map_with_spectrogram")
            fig.savefig(sta+".png")
    map_maker(stations)
    #grid=pygmt.datasets.load_earth_relief(resolution='30s', region=[175,240, 50, 73])
    

def map_maker(stations):
    os.chdir("/Users/sebinjohn/AON_PROJECT/Map_with_spectrogram")
    lon=[]
    lati=[]
    pics=[]
    for j in range(len(stations)):
        sta=stations[j]
        for i in range(len(stationo)):
            if stationo[i]=="final_"+sta+".npy":
                lon.append(long[i])
                lati.append(lat[i])
                pics.append(glob.glob(sta+".jpg"))
    fig=pygmt.Figure()
    #pygmt.makecpt(cmap="abyss.cpt", series=[-7977,0, 1],output="topo.cpt")
    #pygmt.makecpt(cmap="bhw1_deadwood.cpt", series=[0,4000, 1],output="topo.cpt",reverse=True)
    fig.grdimage(grid=grid,region=[175,240, 50, 73],projection="Cyl_stere/200/60/20c",cmap="final_cpt_topo_bathy_handmade.cpt")
    fig.coast(region=[175,240, 50, 73], projection="Cyl_stere/200/60/20c",shorelines=True,frame="f",borders=["1/0.5p,black", "2/0.5p,red"],area_thresh='600')
    for i in range(len(stations)):
        fig.image("legend-01.png", region=[175,240, 50, 73], projection="Cyl_stere/200/60/20c", position="g"+str(180+(180+float(lon[i])))+"/"+lati[i]+"+w1c+jCM")
        fig.image(stations[i]+".png", region=[175,240, 50, 73], projection="Cyl_stere/200/60/20c", position="g"+str(180+(180+float(lon[i])))+"/"+lati[i]+"+w5c+jRM")
    #fig.coast(region=[180,260, 50, 73], projection="D220/60/65/73/12c",shorelines=True,frame="f",borders=["1/0.5p,black", "2/0.5p,red"],area_thresh='600')
    fig.show(method="external")
    #fig.savefig("map_with_spect.png")
        

        
stations=["K13K","A22K","C27K","F21K","S31K","R16K","M27K","I27K"]
plotter(stations)
