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
import cdsapi

c = cdsapi.Client()

os.chdir("")




c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': 'significant_height_of_combined_wind_waves_and_swell',
        'year': '2018',
        'month': '01',
        'day': '01',
        'time': '00:00',
        'format': 'grib',
    },
    'download.grib')



grib="/Users/sebinjohn/Downloads/adaptor.mars.internal-1638789059.0528095-23607-6-cac177c2-831e-45e6-8366-424c3cd6dd2a.grib"
grbs=pygrib.open(grib)
grb=grbs[1]
grb



data = grb.values
latg, long = grb.latlons()
np.shape(data)
lat=(latg[:,0])
lon=(long[0,:])



grid = xr.DataArray(
    data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon}
)


pygmt.makecpt(cmap="bathy",output="expll1.cpt", series=[np.amin(data),np.amax(data)],reverse=True)


fig = pygmt.Figure()
fig.coast(region="g", projection="X6i/6i",shorelines=True)
fig.grdimage(grid=grid,projection="X6i/6i",frame="a", cmap="expll1.cpt",nan_transparent=True)
fig.colorbar(projection="X6i/6i",frame="a", cmap="expll1.cpt")
#fig.grdimage(grid=grid,projection="S200/90/20c", frame="a", cmap="expll.cpt",nan_transparent=True)
fig.show(method="external")


