#--------------------------------------------------------------------------------#
# WRF4PALM to process data from WRF to PALM v6.0
# Output of this script is the NetCDF dynamic driver for PALM following
# PALM Input Data Standard (PIDS) v1.9
#
# Users must provide PALM domain configuration first and run the create_cfg.py script 
#
# Users must provide their own WRF output file.
#
# Users must define the start/end date and time step for the dynamic driver.
#
# Users now can define streched vertical grid spacing
#
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
# Acknowledgement: The author would like to acknowledge Ricardo Faria for his initial
# contribution of WRF2PALM https://github.com/ricardo88faria/WRF2PALM.

#sreenath updated for HRRR data
#updating the code for lod=1, except for soil data
#--------------------------------------------------------------------------------#


#%% clear environment
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#%%

import os
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master')

import gc
import numpy as np
from netCDF4 import Dataset, num2date
from wrf import getvar, destagger, interplevel
import time
import pandas as pd
import xarray as xr
from tqdm import tqdm
from util.geostrophic import *
from util.nearest import nearest
from datetime import datetime
from util.surface_nan_solver import *
from util.interp_array import *
start = datetime.now()
from metpy.interpolate import interpolate_1d, log_interpolate_1d

#%%
###############################################################################
##-------------------------------- User INPUT -------------------------------##
###############################################################################

os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master')

case_name = 'ches_90m_nz180'

#wrf_file = 'wrf_output/wrfout_d04_2017-02-10_00.nc'
wrf_file = 'wrf_output/hrrr_check2_20210421.nc'


interp_mode = 'linear'

# !!! give start and end time to interpolate WRF output here !!!
# this depends on
# 1) WRF output time frequency
# 2) the desired PALM input update frequency
# Time in UTC
dt_start = datetime(2019, 9, 24,22,)
dt_end = datetime(2019, 9, 26, 5,)
interval = 1
ts = '1hour'

# layers for soil temperature and moisture calculation in the LES model
# this shall be changed depending on different cases
# setting it to PALM default now
dz_soil = np.array([0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86])
# define z_origin in meters
z_origin = 0 
###############################################################################
##--------------------------- stretch vertically ----------------------------##
###############################################################################
# stretch factor for a vertically stretched grid
# set this to 1 if no streching required
#dz_stretch_factor = 1.02
dz_stretch_factor = 1.08

# Height level above which the grid is to be stretched vertically (in m)
#dz_stretch_level = 1200
dz_stretch_level = 1800

# aloowed maximum vertical grid spacing (in m)
dz_max = 60

def calc_stretch(z, dz):
    dz_lvl = np.zeros_like(z)
    dz_lvl[:] = dz
    z_stretch = np.copy(z)
    zw_stretch = np.copy(zw)
    for idz, height in enumerate(z):
        if height>dz_stretch_level:
            dz_lvl[idz] = dz_lvl[idz-1]*dz_stretch_factor
            if dz_lvl[idz]<=dz_max:
                z_stretch[idz] = z_stretch[idz-1]+dz_lvl[idz]
            else:
                z_stretch[idz] = z_stretch[idz-1]+dz_max
    for i in range(0, zw.shape[0]):
        zw_stretch[i] = (z_stretch[i]+z_stretch[i+1])*0.5
    return(z_stretch, zw_stretch)

###############################################################################
##--------------------------------- Read CFG --------------------------------##
###############################################################################

cfg = pd.read_csv('cfg_input/'+case_name + '.cfg')
dx = cfg.dx.values[0]
dy = cfg.dy.values[0]
dz = cfg.dz.values[0]
nx = cfg.nx.values[0]
ny = cfg.ny.values[0]
nz = cfg.nz.values[0]
north = cfg.north.values[0]
south = cfg.south.values[0]
east = cfg.east.values[0]
west = cfg.west.values[0]
lat_palm = cfg.centlat[0]
lon_palm = cfg.centlon[0]

y = np.arange(dy/2,dy*ny+dy/2,dy)
x = np.arange(dx/2,dx*nx+dx/2,dx)
z = np.arange(dz/2, dz*nz, dz) + z_origin
xu = x + np.gradient(x)/2
xu = xu[:-1]
yv = y + np.gradient(y)/2
yv = yv[:-1]
zw = z + np.gradient(z)/2
zw = zw[:-1]

if dz_stretch_factor>1:
    z, zw = calc_stretch(z, dz)


# lat_s = ds_wrf['gridlat_0'][:,0].data
# lon_s = ds_wrf['gridlon_0'][0,:].data

###############################################################################
##---------------------------- read WRF variables ---------------------------##
###############################################################################


print(f'Loading WRF netCDF: {wrf_file}' )
#ds_wrf = xr.open_dataset(wrf_file)
#create the check_ds dataset using the combine_nc_v2 code
#check_ds = check_ds.sortby('lv_ISBL0',ascending=False).copy()

ds_wrf = check_ds.copy() #
#nc_wrf = Dataset(wrf_file, 'r')

# ds_wrf.close()
# nc_wrf.close()

# lat_s = ds_wrf['XLAT'][0,:,0].data
# lon_s = ds_wrf['XLONG'][0,0,:].data
###########################################################3
lat_s = ds_wrf['gridlat_0'][:,0].data
lon_s = ds_wrf['gridlon_0'][0,:].data

south_idx, north_idx = nearest(lat_s, south)[1], nearest(lat_s, north)[1]
west_idx, east_idx = nearest(lon_s, west)[1], nearest(lon_s, east)[1]

#lat_v = ds_wrf['XLAT_V'][0,:,0].data

lat_v = ds_wrf['gridlat_0'][:,0].data
southv_idx, northv_idx = nearest(lat_v, south)[1], nearest(lat_v, north)[1]

#lon_u = ds_wrf['XLONG_U'][0,0,:].data
lon_u = ds_wrf['gridlon_0'][0,:].data
westu_idx, eastu_idx = nearest(lon_u, west)[1], nearest(lon_u, east)[1]

# lat_wrf = ds_wrf['XLAT'][0,south_idx:north_idx,0].data
# lon_wrf = ds_wrf['XLONG'][0,0,west_idx:east_idx].data

lat_wrf = ds_wrf['gridlat_0'][south_idx:north_idx,0].data
lon_wrf = ds_wrf['gridlon_0'][0,west_idx:east_idx].data

# If PALM domain smaller than one WRF grid
if north_idx-south_idx<=1:
    north_idx = south_idx+2
if east_idx-west_idx<=1:
    east_idx = west_idx+2
if northv_idx-southv_idx<=1:
    northv_idx = southv_idx+2
if eastu_idx-westu_idx<=1:
    eastu_idx = westu_idx+2

#get soil layer depths and thicknesses from the hrrr data

# depth of the center of soil layer
#zs = ds_wrf['ZS'][0, :].data
#ds_wrf.lv_DBLL7_l1[0].data

##Note, the hrrr soil data start with depth 0 ~ skin properties?
a = 0
zs = []
for i, d in enumerate(ds_wrf.lv_DBLL7_l1[0].values):
    zs.append(np.round((a + d)/2,3))
    a = d
    print(zs)
del a
zs = np.array(zs)    

# thickness of soil layer
#dzs = ds_wrf['DZS'][0, :].data
a = 0
dzs = []
for i, d in enumerate(ds_wrf.lv_DBLL7_l1[0].values):
    dzs.append(np.round((d-a),3))
    a = d
    print(dzs)
del a
dzs = np.array(dzs)    

##Note HRRR data, first soil level is at 0.0

# landmask - 1 is land, 0 is water
#landmask = ds_wrf['LANDMASK'][0, south_idx:north_idx, west_idx:east_idx].data
landmask = ds_wrf['LAND_P0_L1_GLC0'][0,south_idx:north_idx, west_idx:east_idx].data


# TMN - soil temperature at lower boundary
#tmn = ds_wrf['TMN'][0, south_idx:north_idx, west_idx:east_idx].data
#no sol temp at the lower boundary, use soil temp at soil level 0 
tmn = ds_wrf['TSOIL_P0_2L106_GLC0'][0,0,south_idx:north_idx, west_idx:east_idx].data


#incoming radiation at the surface
rad_lw_in = ds_wrf['DLWRF_P0_L1_GLC0'][:,south_idx:north_idx, west_idx:east_idx].data
rad_sw_in = ds_wrf['DSWRF_P0_L1_GLC0'][:,south_idx:north_idx, west_idx:east_idx].data

# # Staggered ASL taking directly from WRF
# PH = ds_wrf['PH']
# PHB = ds_wrf['PHB']
# H_stag = (PH + PHB) / 9.81
# H = destagger(H_stag, stagger_dim=1)

# # De-Staggered ASL taking directly from HRRR
# H = ds_wrf['HGT_P0_L100_GLC0']
# #we don't have H_stag in HRRR so setting it equal to H
# H_stag =  H


#wrf_time = nc_wrf.variables['XTIME']


# time_var = num2date(nc_wrf.variables['XTIME'][:], nc_wrf.variables['XTIME'].units)
# time_var = np.array(time_var).astype('datetime64[s]')

#adding the time var directly from the combined nc file


time_var = ds_wrf.initial_time0_hours.values.astype('datetime64[s]')
for i in range(0,time_var.shape[0]):
    if time_var[i] == np.datetime64(dt_start):
        start_idx = i
        
    if time_var[i] == np.datetime64(dt_end):
        end_idx = i

time_idx = np.arange(start_idx,end_idx+1,interval)

# round up the end time index so that PALM doesn't crash when the final time step is not given
input_lag = (dt_end-dt_start).total_seconds()
tmp_lag = (time_var[time_idx[-1]]-time_var[time_idx[0]]).astype('float')
if input_lag-tmp_lag > 0:
     time_idx= np.append(time_idx,end_idx)


#%%

ds1 = check_ds.copy()
#ds1 = ds_wrf.copy()
#ds1 = ds1.mean('xgrid_0').mean('ygrid_0')

#drop the first two pressure levels and update the third one with the surface values
#because the bottom two pressure values are below surface at the Ches site

ds2 = ds1.where( (ds1.lv_ISBL0!=(101320)) & (ds1.lv_ISBL0!=(100000)), drop=True)

#update surface height and offset

#getting the height of the first pressure level at all times into an array
#note: this first pressure level IS NOT surface pressure, it is at 97500 Pa
Ht_surf = ds2.HGT_P0_L100_GLC0[:,0].values 

#creating a dummy array to fill in with updated heights from the surface
#so that the first data point's value is at 0
H = np.empty_like(ds2.HGT_P0_L100_GLC0.values)

for i in tqdm(range(0,len(ds2.initial_time0_hours)),ascii=True):
    
    #updating surface height values
    ds2.HGT_P0_L100_GLC0[i,0].values = Ht_surf[i]
    for j in range(ds2.HGT_P0_L100_GLC0.values.shape[1]):
        #Updating the actual height by subtracting the surface values
        #creating a 1D array of model heights, starting at the surface for all times
        H[i,j] = ds2.HGT_P0_L100_GLC0[i,j].values - Ht_surf[i]
        
#%%
#now need to update the profiles of the geostrophic winds
#need 3d pres and temp data for it :(

#pull in pressure and temp data to arrays
#run the previous code from create dynamic
pres = np.empty((time_idx.shape[0],H.shape[1]+2,north_idx-south_idx,east_idx-west_idx))
tk = np.empty((time_idx.shape[0],H.shape[1]+2,north_idx-south_idx,east_idx-west_idx))

z_wrf = np.empty((time_idx.shape[0],H.shape[1],north_idx-south_idx,east_idx-west_idx))
times = []


for t, wrf_t in enumerate(time_idx):
    #wrf_t=1
    #t = 1
    
    #pres[t,:,:,:] = getvar(nc_wrf, 'pres', timeidx = wrf_t, units='Pa')[:,south_idx:north_idx, west_idx:east_idx]
    pres[t,:,:,:] = ds_wrf.Pres[wrf_t,:,south_idx:north_idx, west_idx:east_idx]
    # #temp in kelvin
    # #tk[t,:,:,:] = getvar(nc_wrf, 'tk', timeidx = wrf_t)[:,south_idx:north_idx, west_idx:east_idx]
    tk[t,:,:,:] = ds_wrf.TMP_P0_L100_GLC0[wrf_t,:,south_idx:north_idx, west_idx:east_idx]
    z_wrf[t,:,:,:] = H[wrf_t,:, south_idx:north_idx, west_idx:east_idx]
    
   
    times = np.append(times, time_var[t].astype(datetime))
#nc_wrf.close()
print("WRF output reading done.",flush=True)

time_step_sec = ((times[1]-times[0])).total_seconds()
times_sec = np.zeros(time_idx.shape[0])
for t, wrf_t in enumerate(time_idx):
    times_sec[t] = (time_var[wrf_t]-time_var[time_idx[0]]).astype('timedelta64[ns]')
    
#%% 

pres_hint = np.empty((times.shape[0], pres.shape[1]-2,y.shape[0],x.shape[0]))
tk_hint = np.empty((times.shape[0], tk.shape[1]-2,y.shape[0],x.shape[0]))
z_wrf_hint = np.empty((times.shape[0], z_wrf.shape[1],y.shape[0],x.shape[0]))

#have to update the surface pressure and temp values
#don't consider the first two entries from the arrays and 
#sub surface values for the third
for t in tqdm(range(pres.shape[0]),ascii=True):
    pres[t,2] = check_ds.PRES_P0_L1_GLC0[t,south_idx:north_idx, west_idx:east_idx].values
    tk[t,2] = check_ds.TMP_P0_L1_GLC0[t,south_idx:north_idx, west_idx:east_idx].values

#for pressure and temp, only use values from index 2 onwards
# horizontal interpolation
for t in tqdm(range(pres.shape[0]),ascii=True):
    for z_idx in range(2,pres.shape[1]):
        #z_idx =2
        pres_hint[t,z_idx-2,:,:] = interp_array_2d(pres[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        tk_hint[t,z_idx-2,:,:] = interp_array_2d(tk[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        z_wrf_hint[t,z_idx-2,:,:] =  interp_array_2d(z_wrf[t,z_idx-2,:,:], x.shape[0], y.shape[0], interp_mode)

#%%        
pres_tmp = np.empty((times.shape[0], z.shape[0],pres_hint.shape[2], pres_hint.shape[3]))
tk_tmp = np.empty((times.shape[0], z.shape[0],tk_hint.shape[2], tk_hint.shape[3]))

for l_idx, l in tqdm(enumerate(z), desc="Interpolating unstaggered vertical levels"):
    for t in range(0, times.shape[0]):
        pres_tmp[t, int(l_idx), :, :] = interplevel(pres_hint[t, :, :, :], z_wrf_hint[t,:,:,:], l).data
        tk_tmp[t, int(l_idx), :, :] = interplevel(tk_hint[t, :, :, :], z_wrf_hint[t,:,:,:], l).data
        
        
#%%


# calculate geostrophic winds at every levels
# latitudes and longitudes are still required here
def rolling_mean(var, window):
    roll_mean = []
    for i in range(0, var.shape[0] - window, window):
        roll_mean.append(np.nansum(var[i:i + window]) / window)
    return (np.array(roll_mean))

if dz > 10:
    lat_wrf_f = interp_array_1d(lat_wrf,y.shape[0])
    lon_wrf_f = interp_array_1d(lon_wrf,x.shape[0])
    
    geo_wind_u = np.zeros((pres_tmp.shape[0], pres_tmp.shape[1]))
    geo_wind_v = np.zeros((pres_tmp.shape[0], pres_tmp.shape[1]))
    geo_wind_u_f = np.zeros((pres_tmp.shape[0], z.shape[0]))
    geo_wind_v_f = np.zeros((pres_tmp.shape[0], z.shape[0]))
    for t in tqdm(range(pres_tmp.shape[0]),ascii=True, desc="Calculating geostropihc winds"):
        for h in range(0, pres_tmp.shape[1]):
            geo_wind = geostr(pres_tmp[t, h, :, :], tk_tmp[t, h, :, :], lat_wrf_f[:], lon_wrf_f[:], dy, dx)
            geo_wind_u[t, h] = geo_wind[0]
            geo_wind_v[t, h] = geo_wind[1]
         # "smooth" the geostrophic winds after calculation by taking rolling mean
        geo_wind_u_f[t, :] = interp_array_1d(rolling_mean(geo_wind_u[t, :], 10), z.shape[0])
        geo_wind_v_f[t, :] = interp_array_1d(rolling_mean(geo_wind_v[t, :], 10), z.shape[0])
    
    print(flush=True)
    print("Geostrophic wind calculation done.",flush=True)

#check the interpolated data
if False :  
    #Uearth_geo
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,geo_wind_u_f.shape[0],1): 
        ax.plot(geo_wind_u_f[t],z,'-',label=str(t))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    #plt.xlim((-0.02, 0.02))
    #check the interpolated data
    #Vearth_geo
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,geo_wind_v_f.shape[0],1): 
        ax.plot(geo_wind_v_f[t],z,'-',label=str(t))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
#plt.xlim((-0.02, 0.02))

#%%

# vertical interpolation to PALM grid
#linearly interpolating the wind values from ~ 250 m to the surface


Ugearth_int = np.zeros((geo_wind_u_f.shape[0],geo_wind_u_f.shape[1]))

for t in range(0,Ugearth_int.shape[0]): 
    fx = np.polyfit(z[22:63], geo_wind_u_f[t,22:63],2)
    f = np.poly1d(fx)
    Ugearth_int[t,:22] = f(z[:22])   # use interpolation function returned by `interp1d`
    Ugearth_int[t,22:] = geo_wind_u_f[t,22:]   # use interpolation function returned by `interp1d`
    
#try plotting as a check
if False :
    #Uearth_geo
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,Ugearth_int.shape[0],1): 
        #ax.plot(interp(z, z[22:], geo_wind_u_f[t,22:]),z,'-',label=str(t))#backup
        ax.plot(Ugearth_int[t],z,'-',label=str(t))#backup
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 6000.))
    #plt.xlim((-0.02, 0.02))


# vertical interpolation to PALM grid
Vgearth_int = np.zeros((geo_wind_u_f.shape[0],geo_wind_u_f.shape[1]))

for t in range(0,Ugearth_int.shape[0]): 
    fx = np.polyfit(z[22:63], geo_wind_v_f[t,22:63],2)
    f = np.poly1d(fx)
    Vgearth_int[t,:22] = f(z[:22])   # use interpolation function returned by `interp1d`
    Vgearth_int[t,22:] = geo_wind_v_f[t,22:]   # use interpolation function returned by `interp1d`
    
#try plotting as a check
if False:
    #Vearth_geo
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,Vgearth_int.shape[0]): 
        #ax.plot(interp(z, z[22:], geo_wind_u_f[t,22:]),z,'-',label=str(t))#backup
        ax.plot(Vgearth_int[t],z,'-',label=str(t))#backup
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 6000.))
    #plt.xlim((-0.02, 0.02))

#%%
#delete all the 3D data now
del pres_tmp, pres_hint, tk_hint, tk_tmp, z_wrf, z_wrf_hint, tk, pres, H, Ht_surf
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master\wrf_output')


#after calculating the geostrophic wind, average the forcings data for L1 dynamic file

#go on with the rest of the data processing
#subset the HRRR data over PALM domain
check_ds = check_ds.isel(ygrid_0=slice(south_idx,north_idx),xgrid_0=slice(west_idx,east_idx)).copy()

#average it for 1D profiles
check_ds = check_ds.mean('xgrid_0').mean('ygrid_0')


#after averaging
check_ds.to_netcdf('hrrr_v20210722_201909242200-2019092605UTC.nc')
check_ds.close()
ds1.close()
ds2.close()
# test data
#check_ds = xr.open_dataset('hrrr_check2.nc')
check_ds = xr.Dataset()
ds1 = xr.Dataset()
ds2= xr.Dataset()

#reading data back again to convert from dask arrays to xarrays
ds1 = xr.open_dataset('hrrr_v20210722_201909242200-2019092605UTC.nc')

#%% updating and calculating the pt and q 1D arrays

#drop the first two pressure levels and update the third one with the surface values

ds2 = ds1.where( (ds1.lv_ISBL0!=(101320)) & (ds1.lv_ISBL0!=(100000)), drop=True)

#update the surface values for pt,q, pres
P0 = 100000 #Pa
Pt_surface = ds2.TMP_P0_L1_GLC0[:,0]*((P0/ds2.PRES_P0_L1_GLC0[:,0])**(0.286))

for i in range(0,len(ds2.initial_time0_hours)):
    ds2.Pt[i,0] = Pt_surface[i]

#q
for i in range(0,len(ds2.initial_time0_hours)):
    ds2.SPFH_P0_L100_GLC0[i,0] = ds2.SPFH_P0_L103_GLC0[i,0]

    
# #update first pressure level, don't do will have to add time as a coordinate!
# for i in range(0,len(ds2.initial_time0_hours)):
#     ds2.lv_ISBL0[i,0] = Pt_surface[i]


#update surface height and offset
Ht_surf = ds2.HGT_P0_L100_GLC0[:,0].values
H = np.empty_like(ds2.HGT_P0_L100_GLC0.values)

for i in range(0,len(ds2.initial_time0_hours)):
    ds2.HGT_P0_L100_GLC0[i,0].values = Ht_surf[i]
    for j in range(ds2.HGT_P0_L100_GLC0.values.shape[1]):
        H[i,j] = ds2.HGT_P0_L100_GLC0[i,j].values - Ht_surf[i]

#check the plots if needed
if False :
    #check the plots

    #pt
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for i in range(0,ds2.initial_time0_hours.shape[0],2): 
        ax.plot(np.asarray(ds2.Pt.sel(initial_time0_hours=ds2.initial_time0_hours[i])),H[i],'-',label=str(i))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    ax.set_xlabel('theta updated (m)',fontsize=18)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    plt.xlim((275., 310.))
    
    
    #pt old
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for i in range(0,ds1.initial_time0_hours.shape[0],2): 
        ax.plot(np.asarray(ds1.Pt.sel(initial_time0_hours=ds1.initial_time0_hours[i])),np.asarray(ds1.HGT_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),'-',label=str(i))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    ax.set_xlabel('theta wrong',fontsize=18)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    plt.xlim((275., 310.))
    
    #q
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for i in range(1,len(ds1.initial_time0_hours),1): 
        ax.plot(np.asarray(ds2.SPFH_P0_L100_GLC0.sel(initial_time0_hours=ds2.initial_time0_hours[i])),H[i],'-',label=str(pd.to_datetime(ds2.initial_time0_hours[i].values).hour))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    ax.set_xlabel('q updated',fontsize=18)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    #plt.xlim((275., 310.))
    
    #q old
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for i in range(1,len(ds1.initial_time0_hours),1): 
        ax.plot(np.asarray(ds1.SPFH_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),np.asarray(ds1.HGT_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),'-',label=str(pd.to_datetime(ds1.initial_time0_hours[i].values).hour))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    ax.set_xlabel('q wrong',fontsize=18)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    #plt.xlim((275., 310.))


#interpolate the pt and q data to the dynamic file vertical levels
#get the height data to a variable called H
#H = ds2['HGT_P0_L100_GLC0']

#pt interpolation
pt_int = np.zeros((ds2.Pt.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,pt_int.shape[0]): 
    pt_int[t] = np.interp(z, H[t], ds2.Pt[t])
    
#check the interpolated data
if False :
    #pt
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,pt_int.shape[0],2): 
        ax.plot(pt_int[t],z,'-',label=str(t))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    plt.xlim((275., 310.))
    #looks cool!   

#qv interpolation
qv_int = np.zeros((ds2.SPFH_P0_L100_GLC0.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,qv_int.shape[0]): 
    qv_int[t] = np.interp(z, H[t], ds2.SPFH_P0_L100_GLC0[t])
    
#check the interpolated data
#qv
if False:
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,pt_int.shape[0],2): 
        ax.plot(qv_int[t],z,'-',label=str(t))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    #plt.xlim((275., 310.))
    #looks cool!   
        
#%% #interpolate Uearth and Vearth to the PALM grid

#Uearth interpolation
Uearth_int = np.zeros((ds2.Uearth.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,qv_int.shape[0]): 
    Uearth_int[t] = np.interp(z, H[t], ds2.Uearth[t])


#repeat the same for Vearth

#Vearth interpolation
Vearth_int = np.zeros((ds2.Vearth.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,qv_int.shape[0]): 
    Vearth_int[t] = np.interp(z, H[t], ds2.Vearth[t])


#check the interpolated data
if False:
    
    #u
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,pt_int.shape[0],1): 
        ax.plot(Uearth_int[t],z,'-',label=str(t))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    #plt.xlim((275., 310.))
    #looks cool!   

    #v
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,pt_int.shape[0],1): 
        ax.plot(Vearth_int[t],z,'-',label=str(t))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    #plt.xlim((275., 310.))
    #looks cool!   
    
#for vertical velocity, pulling in the updated surface values and keeping it at that

# vertical interpolation to PALM grid
WVEL_int = np.zeros((H.shape[0],zw.shape[0]))

#define zw from the creat_dynamic code
#interpolate to zw
for t in range(0,Uearth_int.shape[0]): 
    WVEL_int[t] = np.interp(zw, H[t], ds2.WVEL[t])

#check the interpolated data
#WVEL
if False:
    fig = plt.figure(figsize=(8,8))
    ax = fig.subplots(1,1)
    for t in range(0,WVEL_int.shape[0],2): 
        ax.plot(WVEL_int[t],zw,'-',label=str(t))
        ax.set_ylabel('Geopotential height (m)',fontsize=18)
    fig.legend(fontsize=12)
    #plt.gca().invert_yaxis()
    plt.ylim((0., 2000.))
    #plt.xlim((-0.02, 0.02))
    #looks cool!   

#%% read in the 1D radiation arrays

#incoming radiation at the surface
rad_lw_in = np.asarray(ds_wrf['DLWRF_P0_L1_GLC0'][:,south_idx:north_idx, west_idx:east_idx].mean('xgrid_0').mean('ygrid_0'))
rad_sw_in = np.asarray(ds_wrf['DSWRF_P0_L1_GLC0'][:,south_idx:north_idx, west_idx:east_idx].mean('xgrid_0').mean('ygrid_0'))



#%%##############################################################################

# Write to NetCDF file
# Based on INIFOR format
print('Writing NetCDF file',flush=True)
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master')

nc_output = xr.Dataset()
res_origin = str(dx) + 'x' + str(dy) + ' m'
nc_output.attrs['description'] = f'Contains dynamic data from NCEP HRRR mesoscale. HRRR output file: {wrf_file}'
nc_output.attrs['author'] = 'Sreenath Paleri (paleri@wisc.edu)'
nc_output.attrs['history'] = 'Created at ' + time.ctime(time.time())
nc_output.attrs['source']= 'netCDF4 python'
nc_output.attrs['origin_lat'] = np.float(lat_palm)
nc_output.attrs['origin_lon'] = np.float(lon_palm)
nc_output.attrs['z'] = np.float(0)
nc_output.attrs['x'] = np.float(0)
nc_output.attrs['y'] = np.float(0)
nc_output.attrs['rotation_angle'] = np.float(0)
nc_output.attrs['origin_time'] =  str(times[0]) + ' UTC'
nc_output.attrs['end_time'] =  str(times[-1]) + ' UTC'


nc_output['x'] = xr.DataArray(x, dims=['x'], attrs={'units':'m'})
nc_output['y'] = xr.DataArray(y, dims=['y'], attrs={'units':'m'})
nc_output['z'] = xr.DataArray(z-z_origin, dims=['z'], attrs={'units':'m'})
nc_output['zsoil'] = xr.DataArray(dz_soil, dims=['zsoil'], attrs={'units':'m'})
nc_output['xu'] = xr.DataArray(xu, dims=['xu'], attrs={'units':'m'})
nc_output['yv'] = xr.DataArray(yv, dims=['yv'], attrs={'units':'m'})
nc_output['zw'] = xr.DataArray(zw-z_origin, dims=['zw'], attrs={'units':'m'})
nc_output['time'] = xr.DataArray(times_sec, dims=['time'], attrs={'units':'seconds'})
nc_output['time_rad'] = xr.DataArray(times_sec, dims=['time_rad'], attrs={'units':'seconds'})

#remember to change this output file name accordingly at each run
#nc_output.to_netcdf(f'dynamic_files/{case_name}_dynamic_{ts}')

#not adding the code for initial profiles and soil moisture or soil temp now.
#add those in for the other test cases later

nc_output['rad_lw_in'] = xr.DataArray(rad_lw_in, dims=['time_rad'], 
         attrs={'units':'W/m^2','lod':np.int32(1), 'source':'WRF', 'long_name':'Downward long-wave rad. flux'}) 

nc_output['rad_sw_in'] = xr.DataArray(rad_sw_in, dims=['time_rad'], 
         attrs={'units':'W/m^2','lod':np.int32(1), 'source':'WRF', 'long_name':'Downward short-wave rad. flux'}) 


nc_output['init_atmosphere_pt'] = xr.DataArray(pt_int[0],dims=['z'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_pt'] = xr.DataArray(pt_int,dims=['time', 'z'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_pt'] = xr.DataArray(pt_int,dims=['time', 'z'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_pt'] = xr.DataArray(pt_int,dims=['time', 'z'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_pt'] = xr.DataArray(pt_int,dims=['time', 'z'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_pt'] = xr.DataArray(pt_int[:,-1],dims=['time'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_qv'] = xr.DataArray(qv_int[0],dims=['z'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_qv'] = xr.DataArray(qv_int,dims=['time', 'z'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_qv'] = xr.DataArray(qv_int,dims=['time', 'z'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_qv'] = xr.DataArray(qv_int,dims=['time', 'z'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_qv'] = xr.DataArray(qv_int,dims=['time', 'z'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_qv'] = xr.DataArray(qv_int[:,-1],dims=['time'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})


nc_output['init_atmosphere_u'] = xr.DataArray(Uearth_int[0],dims=['z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_u'] = xr.DataArray(Uearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_u'] = xr.DataArray(Uearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_u'] = xr.DataArray(Uearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_u'] = xr.DataArray(Uearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_u'] = xr.DataArray(Uearth_int[:,-1],dims=['time'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})


nc_output['init_atmosphere_v'] = xr.DataArray(Vearth_int[0],dims=['z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_v'] = xr.DataArray(Vearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_v'] = xr.DataArray(Vearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_v'] = xr.DataArray(Vearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_v'] = xr.DataArray(Vearth_int,dims=['time', 'z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_v'] = xr.DataArray(Vearth_int[:,-1],dims=['time'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_w'] = xr.DataArray(WVEL_int[0],dims=['zw'],
          attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_w'] = xr.DataArray(WVEL_int,dims=['time', 'zw'],
          attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_w'] = xr.DataArray(WVEL_int,dims=['time', 'zw'],
          attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_w'] = xr.DataArray(WVEL_int,dims=['time', 'zw'],
          attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_w'] = xr.DataArray(WVEL_int,dims=['time', 'zw'],
          attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_w'] = xr.DataArray(WVEL_int[:,-1],dims=['time'],
          attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
    
nc_output['surface_forcing_surface_pressure'] = xr.DataArray(ds2.PRES_P0_L1_GLC0[:,0].values,dims=['time'],
         attrs={'units':'Pa', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})

nc_output['ls_forcing_ug'] = xr.DataArray(Ugearth_int,dims=['time','z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'long_name':'u wind component geostrophic', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_vg'] = xr.DataArray(Vgearth_int,dims=['time','z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'long_name':'v wind component geostrophic', 'source':'WRF', 'res_origin':res_origin})

#remember to change this output file name accordingly at each run
nc_output.to_netcdf(f'dynamic_files/{case_name}_dynamictest_{ts}')

for var in nc_output.data_vars:
    encoding = {var: {'dtype': 'float32', '_FillValue': -9999, 'zlib':True}}
    nc_output[var].to_netcdf(f'dynamic_files/{case_name}_dynamictest_{ts}', encoding=encoding, mode='a')


print('PALM dynamic input file is ready. Script duration: {}'.format(end - start))
print('Start time: '+str(times[0]))
print('End time: '+str(times[-1]))
print('Time step: '+str(time_step_sec)+' seconds')

ds_wrf.close()
gc.collect()


#%% update the time delay if necessary

source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master\dynamic_files'

os.chdir(source)
#dataset for t = 0
ds1 = xr.Dataset()
ds1 = xr.open_dataset('ches_90m_dynamic_1hour.nc')
#just for IOP03 for this time removing data from 22:00 UTC Sep 24th
ds1 = ds1.isel(time=slice(0,17),time_rad=slice(0,17)).copy()

#ds1
#time is both a coordinate and a dimension
#we need to update the values of the coordinate
#create an array with new values
time_new = (np.asarray(range(0,ds1.time.shape[0]))*3600).astype('float64')
#pull in those new values as a dummy coordinate
ds1.coords['time_new'] = time_new
#change the time dimension to a dummy dimension
#following code automatically creates a dummy dimension of integers of the same length
ds1 = ds1.swap_dims({'time': 'dummy'})
#ds1.ls_forcing_left_pt
#drop the time coordinate 
ds1 = ds1.drop('time')
#rename the dummy dimension as time
ds1 = ds1.rename({'dummy': 'time'})
#rename the new time coordinate as time
ds1 = ds1.rename({'time_new': 'time'})


ds1.coords['time_new'] = time_new
#ds1 = ds1.swap_dims({'time': 'time_new'})
ds1 = ds1.swap_dims({'time_rad': 'dummy'})
#ds1.ls_forcing_left_pt
#ds1.coords['time_updated'] = 
ds1 = ds1.drop('time_rad')
#rename the dimension
ds1 = ds1.rename({'dummy': 'time_rad'})
#rename the coordinate
ds1 = ds1.rename({'time_new': 'time_rad'})



#all the other data

ds = xr.open_dataset(f'{case_name}_dynamictest_{ts}')
ds.close()
#ds = nc_output.copy()
#ds.ls_forcing_left_pt
delay = 57600+3600

time_new = (np.asarray(range(0,ds.time.shape[0]))*3600 + delay).astype('float64')

#ds1
#time is both a coordinate and a dimension
#we need to update the values of the coordinate
#pull in those new values as a dummy coordinate
ds.coords['time_new'] = time_new
#change the time dimension to a dummy dimension
#following code automatically creates a dummy dimension of integers of the same length
ds = ds.swap_dims({'time': 'dummy'})
#ds1.ls_forcing_left_pt
#drop the time coordinate 
ds = ds.drop('time')
#rename the dummy dimension as time
ds = ds.rename({'dummy': 'time'})
#rename the new time coordinate as time
ds = ds.rename({'time_new': 'time'})


ds.coords['time_new'] = time_new
#ds1 = ds1.swap_dims({'time': 'time_new'})
ds = ds.swap_dims({'time_rad': 'dummy'})
#ds1.ls_forcing_left_pt
#ds1.coords['time_updated'] = 
ds = ds.drop('time_rad')
#rename the dimension
ds = ds.rename({'dummy': 'time_rad'})
#rename the coordinate
ds = ds.rename({'time_new': 'time_rad'})

ds2 = xr.Dataset()
ds2 = ds.copy()

#%%
ches_dynamic_ds.close()
ches_dynamic_ds = xr.Dataset()
#combine time series
var_list = []
for var in ds1.data_vars:
    if "init" not in var:
        print(var)    
        var_list.append(var)
    


ches_dynamic_ds = xr.merge([ds1[var_list], ds2[var_list]],compat='no_conflicts')#,ds3[var_list],ds4[var_list],ds5[var_list]],compat='no_conflicts')
#ches_dynamic_ds.init_atmosphere_pt

#ches_dynamic_ds.surface_forcing_surface_pressure.time


#pick initial profiles and soil data
for var in ds1.data_vars:
    if "init" in var:
        print(var)    
        ches_dynamic_ds[var] = ds1[var]

#the attributes xu and yv are not present in the ds1 for this run, iop03
#so I'm including them from the earlier version. This won't be necessary for later runs


#%%
#copy the file we prepared now to a separate dataset
test_ds = ches_dynamic_ds.copy()

#read in the current dynamic file for the missing dimensions
#ches_dynamic_ds = xr.open_dataset('ches_90m_dynamic1_1hour')

#pull in those
test_ds['x'] = ches_dynamic_ds.x
test_ds['y'] = ches_dynamic_ds.y
test_ds['xu'] = ches_dynamic_ds.xu
test_ds['yv'] = ches_dynamic_ds.yv

res_origin = '90x90 m'

# test_ds.ls_forcing_left_w.attrs
# test_dynamic_ds.ls_forcing_left_w.attrs

#pulling in initial soil data just to be sure :)
test_ds['init_soil_m'] = ches_dynamic_ds['init_soil_m']
test_ds['init_soil_t'] = ches_dynamic_ds['init_soil_t']

## add attributes
test_ds.attrs['description'] = f'Contains dynamic data from NCEP HRRR mesoscale. HRRR output file'
test_ds.attrs['author'] = 'Sreenath Paleri (paleri@wisc.edu)'
test_ds.attrs['history'] = 'Created at ' + time.ctime(time.time())
test_ds.attrs['source']= 'netCDF4 python'
test_ds.attrs['origin_lat'] = ds1.origin_lat
test_ds.attrs['origin_lon'] = ds1.origin_lon
test_ds.attrs['z'] = np.float(0)
test_ds.attrs['x'] = np.float(0)
test_ds.attrs['y'] = np.float(0)

test_ds.attrs['rotation_angle'] = np.float(0)
#ds = xr.open_dataset('ches_90m_dynamic1_1hour')
test_ds.attrs['origin_time'] =  ds1.origin_time
#ds = xr.open_dataset('ches_90m_dynamic5_1hour')
test_ds.attrs['end_time'] =  ds2.end_time

# test_ds.ls_forcing_vg

test_ds.to_netcdf('ches_90m_dynamic_1hour_52x48_L1_07232021.nc')

#test_dynamic_ds.close()
#test_static_ds.close()
test_ds.close()
ches_dynamic_ds.close()






