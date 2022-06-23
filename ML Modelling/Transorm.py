#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 29 08:37:54 2022

@author: sebinjohn
"""
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os
import pyproj 
import pygmt
import xarray as xr
from mpl_toolkits import mplot3d
from scipy.interpolate import griddata
from global_land_mask import globe


sea_ice_con=np.load("/Users/sebinjohn/AON_PROJECT/Data/sea_ice_con/sea_ice_con.npy")

#sea_ice_con = ma.masked_greater(sea_ice_con, 1.0)
ice=sea_ice_con[:,:,0]


dx = dy = 25000

x = np.arange(-3850000, +3750000, +dx)
y = np.arange(+5850000, -5350000, -dy)

x.shape, y.shape, ice.shape


fig = plt.figure(figsize=(9, 9))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))

cs = ax.coastlines(resolution='110m', linewidth=0.5)

ax.gridlines()
ax.set_extent([-180, 180, 40, 90], crs=ccrs.PlateCarree())

kw = dict(central_latitude=90, central_longitude=-45, true_scale_latitude=70)
cs = ax.pcolormesh(x, y, ice, cmap=plt.cm.Blues,
                   transform=ccrs.Stereographic(**kw))



X,Y=np.meshgrid(x,y)

NorthPolar_WGS=pyproj.Transformer.from_crs(3411,4326)
WGSvalues=NorthPolar_WGS.transform(X,Y)

lat=WGSvalues[0]
lon=WGSvalues[1]

points=np.array((lon.flatten(),lat.flatten())).T

x_s=np.concatenate((np.arange(170,180,0.5),np.arange(-180,-135,0.5)))
y_s=np.arange(55,73,0.5)
grid_x, grid_y = np.meshgrid(x_s,y_s)

sea_ice_ml=np.zeros((36,110,1461))
for i in range(1461):
    values=sea_ice_con[:,:,i].flatten()
    grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
    sea_ice_ml[:,:,i]=grid_z0 
    print(i)

grid_z0.shape

points.shape
values.shape

np.save("sea_ice_downgrid_oandl.npy",sea_ice_ml)

sea_ice_ml=np.load("sea_ice_downgrid_oandl.npy")
def grid_downgrade(grid_x, grid_y):
    points=[]
    for i in range(grid_x.shape[0]):
        for j in range(grid_x.shape[1]):
            points.append((grid_y[i,j],grid_x[i,j]))
    indices=[]
    inpoint=[]
    for point in points:
        if globe.is_ocean(point[0],point[1]):
            inpoint.append(point)
            indices.append(points.index(point))
    return indices
 
ice_grdcut_ocean=np.zeros((len(indices),1461)) 


for i in range(1461):
    ice_grdcut_ocean[:,i]=sea_ice_ml[:,:,i].flatten()[indices]




#ice_grdcut_ocean=np.load("sea_ice_grdcut.npy")
ice_grdcut_ocean_stacked=np.zeros((2268,11688))

xi=np.arange(0,11688)
xp=np.arange(0,11688,8)
for j in range(2268):
    fp=ice_grdcut_ocean[j,:]
    zi=np.interp(xi,xp,fp)
    ice_grdcut_ocean_stacked[j,:]=zi


np.save("sea_ice_grdcut.npy",ice_grdcut_ocean_stacked)
jx=[]
jy=[]

for point in inpoint:
    jx.append(point[1])
    jy.append(point[0])


fig=pygmt.Figure()
fig.coast(region=[150,260,30, 73], projection="L-159/35/33/45/22c",shorelines=False,frame="a",land="200/200/200",borders=["1/0.5p,black", "2/0.5p,red"])
fig.plot(x=jx,y=jy,style="i0.1c", color="black")
fig.show(method="external")






pygmt.makecpt(cmap="haxby",output="ice.cpt", series=[0,1.5,0.005])

fig=pygmt.Figure()
#fig.coast(region=[150,260, 30, 73], projection="L-159/35/33/70/20c",shorelines=True,frame="f",borders=["1/0.5p,black", "2/0.5p,red"],area_thresh='600')
fig.plot(x=jx,y=jy,cmap="ice.cpt",color=ice_grdcut_ocean[:,0].flatten(),style="c0.07c", projection="Cyl_stere/30/-20/12c",region=[170,230, 50, 80])
#fig.grdimage(grid=grid,region=[150,260, 30, 73],projection="Cyl_stere/30/-20/12c", cmap="ice.cpt",frame="f",nan_transparent=True)
fig.coast(region=[170,230, 50, 80], projection="Cyl_stere/30/-20/12c",shorelines=True,frame="f",borders=["1/0.5p,black", "2/0.5p,red"],area_thresh='600',land="grey")
fig.colorbar(projection="Cyl_stere/30/-20/12c", cmap="ice.cpt",frame=["x+lIce\tConcenteration", "y+lm"],position="n0.57/-0.1+w8c/0.5c+h")
fig.show(method="external")




# f = interpolate.interp2d(latA, lonA, ice_A, kind='linear')

# lat_g=np.arange(55,75,0.05)
# lon_g=np.arange(170,230,0.05)
# ice_new = f(lat_g.flatten(), lon_g.flatten())
# ice_new= ma.masked_greater(ice_new, 1.0)
# ice_new=ma.masked_less(ice_new, 0)

#grid = xr.DataArray(ice_new.T, dims=["lat", "lon"], coords={"lat":lat_g, "lon":  lon_g})