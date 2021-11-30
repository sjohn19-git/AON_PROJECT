#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 00:34:28 2021

@author: sebinjohn
"""

import pygmt
import xarray as xr
import pygrib
import numpy as np
import matplotlib.pyplot as plt


grid = pygmt.datasets.load_earth_relief(resolution="01d")
lat=grid["lat"]
lon=grid["lon"]



grib = "/Users/sebinjohn/Downloads/multi_1.ak_4m.hs.201609.grb2"
grib1="/Users/sebinjohn/Downloads/multi_1.ak_10m.hs.201903.grb2"
grbs = pygrib.open(grib)
grb=grbs[1]
grbs1=pygrib.open(grib1)
grb1=grbs1[1]
grb1


data = grb.values
latg, long = grb.latlons()
np.shape(data)
lat=latg[:,0]
lon=long[0,:]
np.shape(lat)
np.shape(lon)
np.shape(data)
np.amax(data)
np.amin(data)



data1 = grb1.values
latg1, long1 = grb1.latlons()
np.shape(data1)
lat1=(latg1[:,0])
lon1=(long1[0,:])



grid1 = xr.DataArray(
    data1, dims=["lat", "lon"], coords={"lat": lat1, "lon": lon1}
)

grid = xr.DataArray(
    data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
)

pygmt.makecpt(cmap="hot",output="expll.cpt", series=[np.amin(data),np.amax(data)],reverse=True)
pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[np.amin(data1),np.amax(data1)],reverse=True)


fig = pygmt.Figure()
fig.coast(region=[-190,-130, 48, 73], projection="S200/90/20c",shorelines=True)
fig.grdimage(grid=grid1,projection="S200/90/20c",frame="a", cmap="expll1.cpt",nan_transparent=True)
# fig.grdimage(grid=grid,projection="S200/90/20c", frame="a", cmap="expll.cpt",nan_transparent=True)


fig.show()

