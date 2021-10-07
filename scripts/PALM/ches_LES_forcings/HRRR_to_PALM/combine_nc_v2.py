# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 01:42:24 2021

@author: Sreenath
"""

#%% clear environment
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#%%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import datetime
#import seaborn as sns
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import dask
from metpy.interpolate import interpolate_1d
#%%

def rho(T, p):
    
    """
    
    Calculates air density (rho)
    
    """
    
    
    Rd = 287.0

#    Tv   = T * (1+0.61*qv) # Virtual temperature

    rho = p / (Rd * T) # Air density [kg m^-3]

    return(rho)
#%% pull out and transform all the necessary variables and write them our as separate files

day = '20190926'

source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\HRRRB\hrrr' + '\\' + day + '\data'

#move into the data folder
os.chdir(source)

#get list of all the files
file_list = os.listdir(source)


for i, file_name in enumerate(file_list):
    if file_name.endswith('.nc'):
        print(file_name)
        print(i)
        #setup the xtime and pres variables
        ds1_check = xr.open_dataset(file_name)
        hour = np.asarray([float(pd.to_datetime(ds1_check.initial_time0_hours[0].values).hour*60)])
        initial_time0_hours = ds1_check['initial_time0_hours'].values
        XTIME = xr.DataArray(hour, coords={'initial_time0_hours': initial_time0_hours},
                     dims=['initial_time0_hours'])
        ds1_check['XTIME'] = XTIME
        ds1_check['XTIME'].attrs['long_name'] = 'minutes from forcings start'
        ds1_check['XTIME'].attrs['units'] = ''
        
        lv_ISBL0 = ds1_check.lv_ISBL0.values 
        dummy = np.empty((40,110,111))
        
        for i in range(np.shape(dummy)[1]):
            for j in range(np.shape(dummy)[2]):
                dummy[:,i,j] = lv_ISBL0
        
        gridlat_0 = ds1_check.gridlat_0.values
        gridlon_0 = ds1_check.gridlon_0.values
        
        Pres = xr.DataArray(dummy, coords={'lv_ISBL0': lv_ISBL0,
                                               'gridlat_0':(['ygrid_0','xgrid_0'],gridlat_0),
                                               'gridlon_0':(['ygrid_0','xgrid_0'],gridlon_0)},
                     dims=['lv_ISBL0','ygrid_0','xgrid_0'])
        del dummy
        
        ds1_check['Pres'] = Pres
        ds1_check.Pres.attrs = ds1_check['Pres'].attrs
        
        
        #rotate wind
        ds1_check['Vearth'] = ds1_check.VGRD_P0_L100_GLC0*np.cos(ds1_check.gridrot_0) - \
                              ds1_check.UGRD_P0_L100_GLC0*np.sin(ds1_check.gridrot_0)
        ds1_check['Uearth'] = ds1_check.VGRD_P0_L100_GLC0*np.sin(ds1_check.gridrot_0) + \
                              ds1_check.UGRD_P0_L100_GLC0*np.cos(ds1_check.gridrot_0)
        
        ds1_check['Vearth_surf'] = ds1_check.VGRD_P0_L103_GLC0*np.cos(ds1_check.gridrot_0) - \
                                   ds1_check.UGRD_P0_L103_GLC0*np.sin(ds1_check.gridrot_0)
        ds1_check['Uearth_surf'] = ds1_check.VGRD_P0_L103_GLC0*np.sin(ds1_check.gridrot_0) + \
                                   ds1_check.UGRD_P0_L103_GLC0*np.cos(ds1_check.gridrot_0)
        #convert to potential temp
        P0 = 100000 #Pa
        ds1_check['Pt'] = ds1_check.TMP_P0_L100_GLC0*((P0/ds1_check.Pres)**(0.286))
        
        #convert to vertical velocity
        g = 9.80665 #m/s**2
        #per http://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml
        ds1_check['WVEL'] =  -(ds1_check['VVEL_P0_L100_GLC0']/ (rho(ds1_check.TMP_P0_L100_GLC0,ds1_check.Pres)*g))
        
        check1 = ds1_check[['lv_DBLL7_l1','lv_DBLL7_l0','lv_HTGL1','LAND_P0_L1_GLC0','TSOIL_P0_2L106_GLC0','HGT_P0_L100_GLC0','HGT_P0_L1_GLC0',
             'XTIME','TMP_P0_L1_GLC0', 'TMP_P0_L100_GLC0','PRES_P0_L1_GLC0','Pres','SOILW_P0_2L106_GLC0','SPFH_P0_L100_GLC0',
             'SPFH_P0_L103_GLC0','Uearth','Vearth','Uearth_surf','Vearth_surf','WVEL','Pt','DSWRF_P0_L1_GLC0','DLWRF_P0_L1_GLC0']]
        
        #average the values horizontally, except for soil data

        # for var in check1.data_vars:
        #     if var in list_to_avg:
        #         #print(var)
        #         check1[var] = check1[var].mean('xgrid_0').mean('ygrid_0')

        check1.to_netcdf(f'hrrr_forcings\{file_name}')
        print(f'{file_name} done')
        check1.close()
        ds1_check.close()

#%% checking the pt conversion
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\HRRRB\hrrr\20190924\data')


with xr.open_dataset('20190924_hrrr.t16z.wrfprsf00.small.nc') as ds_hrrr:
    print(ds_hrrr.keys())
    

data_hrrr     = '20190924_hrrr.t17z.wrfprsf00.small.nc'
ds_hrrr = xr.open_dataset(data_hrrr)
ds_hrrr.close()

ds_hrrr.PRES_P0_L1_GLC0
ds_hrrr.TMP_P0_L1_GLC0
ds_hrrr.HGT_P0_L1_GLC0
ds_hrrr.UGRD_P0_L100_GLC0[:,39].mean('xgrid_0').mean('ygrid_0')
ds_hrrr.UGRD_P0_L103_GLC0.mean('xgrid_0').mean('ygrid_0')

ds1_check = ds_hrrr.copy()

lv_ISBL0 = ds1_check.lv_ISBL0.values 
dummy = np.empty((40,110,111))

for i in range(np.shape(dummy)[1]):
    for j in range(np.shape(dummy)[2]):
        dummy[:,i,j] = lv_ISBL0

gridlat_0 = ds1_check.gridlat_0.values
gridlon_0 = ds1_check.gridlon_0.values

Pres = xr.DataArray(dummy, coords={'lv_ISBL0': lv_ISBL0,
                                       'gridlat_0':(['ygrid_0','xgrid_0'],gridlat_0),
                                       'gridlon_0':(['ygrid_0','xgrid_0'],gridlon_0)},
             dims=['lv_ISBL0','ygrid_0','xgrid_0'])
del dummy

ds1_check['Pres'] = Pres
ds1_check.Pres.attrs = ds1_check.lv_ISBL0.attrs



#ds1_check = ds1_check.sortby('lv_ISBL0',ascending=False).copy()

ds1_check = ds1_check.mean('xgrid_0').mean('ygrid_0')
#ds1_check.lv_ISBL0

# ds1_check.Pres
# ds1_check.Pt
# ds1_check.TMP_P0_L100_GLC0
P0 = 100000 #Pa

ds1_check['Pt'] = ds1_check.TMP_P0_L100_GLC0*((P0/ds1_check.Pres)**(0.286))

#%% checking temperature and pt plots
fig = plt.figure(figsize=(8,8))
ax = fig.subplots(1,1)
ax.plot(np.asarray(ds1_check.TMP_P0_L100_GLC0.sel(initial_time0_hours=ds1_check.initial_time0_hours[0])),np.asarray(ds1_check.HGT_P0_L100_GLC0.sel(initial_time0_hours=ds1_check.initial_time0_hours[0])),'-')
ax.set_ylabel('Geopotential height (m)',fontsize=18)
ax.set_xlabel('Temp (K)',fontsize=18)
#fig.legend(fontsize=12)
#plt.gca().invert_yaxis()
plt.ylim((0., 2000.))
plt.xlim((275., 310.))

fig = plt.figure(figsize=(8,8))
ax = fig.subplots(1,1)
ax.plot(np.asarray(ds1_check.Pt.sel(initial_time0_hours=ds1_check.initial_time0_hours[0])),np.asarray(ds1_check.HGT_P0_L100_GLC0.sel(initial_time0_hours=ds1_check.initial_time0_hours[0])),'-')
ax.set_xlabel('Pot. Temp (K)',fontsize=18)
#fig.legend(fontsize=12)
#plt.gca().invert_yaxis()
plt.ylim((0., 2000.))
plt.xlim((275., 310.))
        


#%% combine all the necessary .nc files

#combine files for day1
day = '20190924'

source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\HRRRB\hrrr' + '\\' + day + '\data\hrrr_forcings'
#source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\HRRRB\hrrr\scratch'
#move into the data folder
os.chdir(source)

#get list of all the files
file_list = os.listdir(source)

forcings_list = []
for i, file_name in enumerate(file_list):
    if (i in range(22,24)): #checking for pt, #starting at midnight CDT on day1 5,24
        print(i)
        forcings_list.append(file_list[i])
        
# #data sanity check
# ds1_check = xr.open_dataset(forcings_list[18])
# ds1_check.close()

# ds.close()
#ds1 = xr.open_dataset(file_list[18])

ds1 = xr.open_mfdataset(forcings_list,combine = 'by_coords', concat_dim="initial_time0_hours")
ds1 = ds1.sortby('lv_ISBL0',ascending=False).copy()


#combine files for day2
day = '20190925'

source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\HRRRB\hrrr' + '\\' + day + '\data\hrrr_forcings'

#move into the data folder
os.chdir(source)

#get list of all the files
file_list = os.listdir(source)
forcings_list = []
for i, file_name in enumerate(file_list):
    #if i in range(0,7):#ending at 01:00 CDT day 2
    if i in range(0,24):#for day2 
        print(i)
        forcings_list.append(file_list[i])

ds2 = xr.open_mfdataset(forcings_list,combine = 'by_coords', concat_dim="initial_time0_hours")    
ds2 = ds2.sortby('lv_ISBL0',ascending=False).copy()
    
#combine files for day3
day = '20190926'

source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\HRRRB\hrrr' + '\\' + day + '\data\hrrr_forcings'

#move into the data folder
os.chdir(source)

#get list of all the files
file_list = os.listdir(source)
forcings_list = []
for i, file_name in enumerate(file_list):
    if i in range(0,6):#ending at midnight CDT day 3
        print(i)
        forcings_list.append(file_list[i])

ds3 = xr.open_mfdataset(forcings_list,combine = 'by_coords', concat_dim="initial_time0_hours")    
ds3 = ds3.sortby('lv_ISBL0',ascending=False).copy()


check_ds = xr.concat([ds1, ds2, ds3], dim="initial_time0_hours")

#check_ds.initial_time0_hours

#check_ds = check_ds.sortby('lv_ISBL0',ascending=False).copy()

#check the plots of Temp and Pt and Q once.
#%% checking temperature and pt plots


# check_ds = check_ds.mean('xgrid_0').mean('ygrid_0')

# os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master')


# #plotting HRRR temperature profiles
# lvl_isb =  np.asarray(ds1['lv_ISBL0'])
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,len(ds1.initial_time0_hours)): 
#     ax.plot(np.asarray(ds1.TMP_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),np.asarray(ds1.HGT_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),'-',label=str(pd.to_datetime(ds1.initial_time0_hours[i].values).hour))
#     ax.set_ylabel('Geopotential height (m)',fontsize=18)
# fig.legend(fontsize=12)
# #plt.gca().invert_yaxis()
# plt.ylim((450., 2500.))
# plt.xlim((275., 310.))

# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,len(ds1.initial_time0_hours),2): 
#     ax.plot(np.asarray(ds1.SPFH_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),np.asarray(ds1.HGT_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),'-',label=str(pd.to_datetime(ds1.initial_time0_hours[i].values).hour))
#     ax.set_ylabel('Geopotential height (m)',fontsize=18)
# fig.legend(fontsize=12)
# #plt.gca().invert_yaxis()
# plt.ylim((0., 2000.))
# #plt.xlim((275., 310.))

# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,len(ds1.initial_time0_hours)): 
#     ax.plot(np.asarray(ds1.Pt.sel(initial_time0_hours=ds1.initial_time0_hours[i])),np.asarray(ds1.HGT_P0_L100_GLC0.sel(initial_time0_hours=ds1.initial_time0_hours[i])),'-',label=str(pd.to_datetime(ds1.initial_time0_hours[i].values).hour))
#     ax.set_ylabel('Geopotential height (m)',fontsize=18)
# fig.legend(fontsize=12)
# #plt.gca().invert_yaxis()
# plt.ylim((500., 2000.))
# plt.xlim((275., 310.))

# P0 = 100000 #Pa

# #check data first
# for i in range(0,len(ds1.initial_time0_hours)):
#     Pt_surface = ds1.TMP_P0_L1_GLC0[i]*((P0/ds1.PRES_P0_L1_GLC0[i])**(0.286))
#     print('Pt surface init: ', ds1.Pt[i,0].values, 'Pt Surface new : ', Pt_surface.values)

# #sub the surface values for the Pt profile.
# ds1['Pt_updated'] = ds1['Pt']
# for i in range(0,len(ds1.initial_time0_hours)):
#     Pt_surface = ds1.TMP_P0_L1_GLC0[i]*((P0/ds1.PRES_P0_L1_GLC0[i])**(0.286))
#     ds1.Pt_updated[i,0].load()
#     (ds1.Pt_updated[i,0].load()).values = Pt_surface.values


# #reading in the ches dynamic file anc checking surface pt values.

# source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master\dynamic_files'

# os.chdir(source)

# ches_dynamic = 'ches_dynamic_test'
# ches_dynamic_ds = xr.open_dataset(ches_dynamic)

# ches_dynamic_ds = ches_dynamic_ds.mean('x').mean('y').mean('yv').mean('xu')

# #sub the values to the dynamic file
# z_forcing = np.asarray(ches_dynamic_ds['z'])
# zw_forcing = np.asarray(ches_dynamic_ds['zw'])

# mean = ches_dynamic_ds[['ls_forcing_left_pt','ls_forcing_right_pt',
#                                             'ls_forcing_south_pt','ls_forcing_north_pt']].to_array(dim='new').mean('new')
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,len(mean.time)): 
#     ax.plot(np.asarray(mean.sel(time=mean.time[i])),z_forcing,'-o',label=str(mean.time[i].data/(1E9*3600)))
#     ax.set_ylim([0., 1500.])
# ax.set_ylabel('z (m)',fontsize=18)
# ax.set_xlabel('theta HRRR (K)',fontsize=18)
# fig.legend(fontsize=12)


# #check data first
# for i in range(0,len(ds1.initial_time0_hours)):
#     Pt_surface = ds1.TMP_P0_L1_GLC0[i]*((P0/ds1.PRES_P0_L1_GLC0[i])**(0.286))
#     print('Pt surface init: ', mean[i,0].values, 'Pt Surface new : ', Pt_surface.values)
    
# #sub the surface values for the Pt profile.
# for i in range(0,len(ds1.initial_time0_hours)):
#     Pt_surface = ds1.TMP_P0_L1_GLC0[i]*((P0/ds1.PRES_P0_L1_GLC0[i])**(0.286))
#     mean[i,0] = Pt_surface.values
#     print('Pt surface init: ', mean[i,0].values, 'Pt Surface new : ', Pt_surface.values)
    
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,len(mean.time),2): 
#     ax.plot(np.asarray(mean.sel(time=mean.time[i])),z_forcing,'-',label=str(mean.time[i].data/(1E9*3600)))
#     ax.set_ylim([0., 1500.])
# ax.set_ylabel('z (m)',fontsize=18)
# ax.set_xlabel('theta HRRR (K)',fontsize=18)
# fig.legend(fontsize=12)

#%% back to main code
#substitute surface temp and q values for 
#if doing part by part
#check_ds = ds1.copy()
#check_ds = ds2.copy()


ds1 = check_ds.copy()
#ds1 = ds1.mean('xgrid_0').mean('ygrid_0')

#drop the first two pressure levels and update the third one with the surface values

ds2 = ds1.where( (ds1.lv_ISBL0!=(101320)) & (ds1.lv_ISBL0!=(100000)), drop=True)

#update surface height and offset
Ht_surf = ds2.HGT_P0_L100_GLC0[:,0].values
H = np.empty_like(ds2.HGT_P0_L100_GLC0.values)

for i in tqdm(range(0,len(ds2.initial_time0_hours)),ascii=True):
    ds2.HGT_P0_L100_GLC0[i,0].values = Ht_surf[i]
    for j in range(ds2.HGT_P0_L100_GLC0.values.shape[1]):
        H[i,j] = ds2.HGT_P0_L100_GLC0[i,j].values - Ht_surf[i]


#now need to update the profiles of the geostrophic winds
#need 3d pres and temp data for it :(
#so do this at the beginning
#need to run the create dynamic script to calculate the temp arrays too
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master')


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

pres_tmp = np.empty((times.shape[0], z.shape[0],pres_hint.shape[2], pres_hint.shape[3]))
tk_tmp = np.empty((times.shape[0], z.shape[0],tk_hint.shape[2], tk_hint.shape[3]))

for l_idx, l in tqdm(enumerate(z), desc="Interpolating unstaggered vertical levels"):
    for t in range(0, times.shape[0]):
        pres_tmp[t, int(l_idx), :, :] = interplevel(pres_hint[t, :, :, :], z_wrf_hint[t,:,:,:], l).data
        tk_tmp[t, int(l_idx), :, :] = interplevel(tk_hint[t, :, :, :], z_wrf_hint[t,:,:,:], l).data


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
for t in range(0,geo_wind_v_f.shape[0],2): 
    ax.plot(geo_wind_v_f[t],z,'-',label=str(t))
    ax.set_ylabel('Geopotential height (m)',fontsize=18)
fig.legend(fontsize=12)
#plt.gca().invert_yaxis()
plt.ylim((0., 2000.))
#plt.xlim((-0.02, 0.02))

# #check previous dynamic forcings data
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,ches_dynamic_ds.time.shape[0],2): 
#     ax.plot(np.asarray(ches_dynamic_ds.ls_forcing_vg.sel(time=ches_dynamic_ds.time[i])),z_forcing,'-',label=str(i))
#     ax.set_ylim([0., 2000.])
# ax.set_ylabel('z (m)',fontsize=18)
# ax.set_xlabel('Vg',fontsize=18)
# fig.legend(fontsize=12)


#looks weird near the surface, below 250m 
#looks like there are issues while interpolating the surface pressure values
#pres_tmp[t, 0, :, :] has some values lower than pres_tmp[t, 1, :, :]
#and z_wrf_hint[t,1,:,:] is ~ 220. So, we can set the first grid value at
#z_wrf = 0 to be 0 and then interpolate to the next z_wrf level to create the final profiles.

#check the test_dynamic file for geostrophic profiles to be sure.
# test_dynamic_urban_ds
#test_dynamic_ds  miskor_dynamic

# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,test_dynamic_ds.time.shape[0],2): 
#     ax.plot(np.asarray(test_dynamic_ds.ls_forcing_vg.sel(time=test_dynamic_ds.time[i])),np.asarray(test_dynamic_ds.z),'-',label=str(i))
#     ax.set_ylim([0., 2000.])
# ax.set_ylabel('z (m)',fontsize=18)
# ax.set_xlabel('Vg',fontsize=18)
# fig.legend(fontsize=12)

#looks like surface values are non zero.
#extract the surface height value from the previous profile
#find ug and vg at that height and linearly interpolate to the z_wrf second grid point value

# #check previous dynamic forcings data
# #Ug
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,ches_dynamic_ds.time.shape[0],2): 
#     ax.plot(np.asarray(ches_dynamic_ds.ls_forcing_ug.sel(time=ches_dynamic_ds.time[i])),z_forcing,'-',label=str(i))
#     ax.set_ylim([0., 2000.])
# ax.set_ylabel('z (m)',fontsize=18)
# ax.set_xlabel('Ug',fontsize=18)
# fig.legend(fontsize=12)

# #Vg
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for i in range(1,ches_dynamic_ds.time.shape[0],2): 
#     ax.plot(np.asarray(ches_dynamic_ds.ls_forcing_vg.sel(time=ches_dynamic_ds.time[i])),z_forcing,'-',label=str(i))
#     ax.set_ylim([0., 1000.])
# ax.set_ylabel('z (m)',fontsize=18)
# ax.set_xlabel('Vg',fontsize=18)
# fig.legend(fontsize=12)


#what are the heights we need the
#looks like should linearly interpolate below a height of z[21]!
#down to the surface

#better than that, sub values from previous calculation, ches_dynamic values
#for height levels below z[21], indices from 0 to 20 - 21 values
#sub values fromz[27] to next 21 values

# geo_wind_u_f[1,21]
# geo_wind_v_f[1,21]


# ches_dynamic_ds.ls_forcing_ug[1,48]

# #does not look promising for now, but let's try substituting and see
# geo_wind_u_f_updated = np.zeros((geo_wind_u_f.shape[0],geo_wind_u_f.shape[1]))
# geo_wind_v_f_updated = np.zeros((geo_wind_v_f.shape[0],geo_wind_v_f.shape[1]))

# for t in range(geo_wind_u_f_updated.shape[0]):
#     #t=0
#     geo_wind_u_f_updated[t,:21] = ches_dynamic_ds.ls_forcing_ug[t,27:48]
#     geo_wind_u_f_updated[t,21:] = geo_wind_u_f[t,21:]

# #check
# #Uearth_geo
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for t in range(0,geo_wind_u_f_updated.shape[0],2): 
#     ax.plot(geo_wind_u_f_updated[t],z,'-',label=str(t))
#     ax.set_ylabel('Geopotential height (m)',fontsize=18)
# fig.legend(fontsize=12)
# #plt.gca().invert_yaxis()
# plt.ylim((0., 2000.))
# #plt.xlim((-0.02, 0.02))

#does not work, just linearly interpolate to surface and be done with it for now.
#del (geo_wind_u_f_updated)

#height for the second model grid point above surface
# z_wrf[0,1].mean()
# #Out[819]: 219.34665256076389

# for t in range(z_wrf.shape[0]):
#     print(z_wrf[t,1].mean())

# # find that height from the interpolated PALM grid
# #= z[18]

# #try plotting from z[18] as a check
# #Uearth_geo
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for t in range(0,geo_wind_u_f.shape[0],2): 
#     ax.plot(geo_wind_u_f[t,22:],z[22:],'-',label=str(t))
#     ax.set_ylabel('Geopotential height (m)',fontsize=18)
# fig.legend(fontsize=12)
# #plt.gca().invert_yaxis()
# plt.ylim((0., 2000.))
# #plt.xlim((-0.02, 0.02))



# vertical interpolation to PALM grid
#linearly interpolating the wind values from ~ 250 m to the surface
Ugearth_int = np.zeros((geo_wind_u_f.shape[0],geo_wind_u_f.shape[1]))

for t in range(0,Ugearth_int.shape[0]): 
    fx = np.polyfit(z[22:63], geo_wind_u_f[t,22:63],2)
    f = np.poly1d(fx)
    Ugearth_int[t,:22] = f(z[:22])   # use interpolation function returned by `interp1d`
    Ugearth_int[t,22:] = geo_wind_u_f[t,22:]   # use interpolation function returned by `interp1d`
    
#try plotting as a check
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
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master\wrf_output')


#finally works!!!!!!


#go on with the rest of the data processing
check_ds = check_ds.mean('xgrid_0').mean('ygrid_0')


#after averaging
check_ds.to_netcdf('hrrr_check2_20210425.nc')
check_ds.close()
ds1.close()
ds2.close()
# test data
#check_ds = xr.open_dataset('hrrr_check2.nc')
check_ds = xr.Dataset()
ds1 = xr.Dataset()
ds2= xr.Dataset()

#reading data back again to convert from dask arrays to xarrays
ds1 = xr.open_dataset('hrrr_check2_20210425.nc')

# #sub the surface values for the Pt profile.
# for i in range(0,len(ds1.initial_time0_hours)):
#     Pt_surface = ds1.TMP_P0_L1_GLC0[i]*((P0/ds1.PRES_P0_L1_GLC0[i])**(0.286))
#     ds1.Pt[i,0] = Pt_surface
    
# #Update the surface values for q
# for i in range(0,len(ds1.initial_time0_hours)):
#     ds1.SPFH_P0_L100_GLC0[i,0] = ds1.SPFH_P0_L103_GLC0[i]

# #Update the lowest geopotential height to 0, for surface
# for i in range(0,len(ds1.initial_time0_hours)):
#     ds1.HGT_P0_L100_GLC0[i,0] = 0

#update the profile data with first data point at the surface
#update the heights by offsetting for surface height
#then pickout the next data point at the next highest pressure level

# #check surface pressures
# for i in range(0,len(ds1.initial_time0_hours)):
#     #check difference between surface pressure and pressure level
#     print(ds1.PRES_P0_L1_GLC0[i].data - ds1.Pres[i,3].data )
#     #4th pressure level at 95000 Pa seems to be AGL always
#     #use that as second grid point

# #check the surface heights
# for i in range(0,len(ds1.initial_time0_hours)):
#     #check for the difference between surface heights and pressure level above
#     print(ds1.HGT_P0_L100_GLC0[i,3].data - ds1.HGT_P0_L1_GLC0[i].data)

# #will have to interpolate for the first 180 m at night!

# #check the surface heights between the first two grid points
# for i in range(0,len(ds1.initial_time0_hours)):
#     #check for the difference between sea level and pressure level at 100000 Pa
#     print(ds1.HGT_P0_L100_GLC0[i,1].data - ds1.HGT_P0_L100_GLC0[i,0].data)
#     #always ~ 110 m


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

#looks cool!   

#interpolate the pt and q data to the dynamic file vertical levels
#get the height data to a variable called H
#H = ds2['HGT_P0_L100_GLC0']

#pt interpolation
pt_int = np.empty((ds2.Pt.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,pt_int.shape[0]): 
    pt_int[t] = np.interp(z, H[t], ds2.Pt[t])
    
#check the interpolated data
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
qv_int = np.empty((ds2.SPFH_P0_L100_GLC0.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,qv_int.shape[0]): 
    qv_int[t] = np.interp(z, H[t], ds2.SPFH_P0_L100_GLC0[t])
    
#check the interpolated data
#pt
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

#now start updating the velocity profiles
#the lowest grid point value needs updating, it is now at 975 hpa, below surface
#try and plot the velocity profiles

#u/v
fig = plt.figure(figsize=(8,8))
ax = fig.subplots(1,1)
for i in range(1,ds1.initial_time0_hours.shape[0],2): 
    ax.plot(np.asarray(ds2.Uearth.sel(initial_time0_hours=ds2.initial_time0_hours[i])),H[i],'-',label=str(i))
    ax.set_ylabel('Geopotential height (m)',fontsize=18)
fig.legend(fontsize=12)
#plt.gca().invert_yaxis()
plt.ylim((0., 2000.))
#plt.xlim((275., 310.))

#check the surface value for velocities in ds2 with the actual surface value
for t in range(0,ds2.initial_time0_hours.shape[0]):
    print(ds2.Uearth_surf[t,0,0].data,ds2.Uearth[t,0].data)
    
#there are differences, better to create a new array with data at 10 and 80 m added
# and a separate height vector just for the velocities. For vertical velocities looks like it
#is enough to set it equal to the 97500 value. Keep it at that
#tried this but the profiles look discontinuous when interpolated
#continue with using just the Uearth and Vearth arrays because surface values
#are not very off

{
# #create a separate height array for velocities with 10 and 80 m as the first two data points
# H_vel = np.zeros((H.shape[0],H.shape[1]+1))
# for t in range(0,ds2.initial_time0_hours.shape[0]):
#     H_vel[t,1:] = H[t].copy()
#     H_vel[t,1] = 80
#     H_vel[t] = np.concatenate(([10],H_vel[t,1:]))



# #for Uearth, creating a data array 
# Uearth = np.zeros((ds2.Uearth.shape[0],ds2.Uearth.shape[1]+1))

# #t=0
# for t in range(0,ds2.initial_time0_hours.shape[0]):
#     #add in the 10m and 80m values with the Uearth values
#     Uearth[t] = np.concatenate(([ds2.Uearth_surf[t,0,0].data,ds2.Uearth_surf[t,1,0].data],ds2.Uearth.data[t,1:]))

# # vertical interpolation to PALM grid
# Uearth_int = np.zeros((H_vel.shape[0],z.shape[0]))

# #define z from the creat_dynamic code
# #interpolate to z
# for t in range(0,Uearth_int.shape[0]): 
#     Uearth_int[t] = np.interp(z, H_vel[t], Uearth[t])

# #check the interpolated data
# #Uearth
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for t in range(0,Uearth_int.shape[0],2): 
#     ax.plot(Uearth_int[t],z,'-',label=str(t))
#     ax.set_ylabel('Geopotential height (m)',fontsize=18)
# fig.legend(fontsize=12)
# #plt.gca().invert_yaxis()
# plt.ylim((0., 2000.))
# #plt.xlim((275., 310.))
# #does not look very cool, continue using just Uearth and Vearth!

####repeat the same for Vearth

# #for Vearth, creating a data array 
# Vearth = np.zeros((ds2.Vearth.shape[0],ds2.Vearth.shape[1]+1))

# #t=0
# for t in range(0,ds2.initial_time0_hours.shape[0]):
#     #add in the 10m and 80 m  value
#     Vearth[t] = np.concatenate(([ds2.Vearth_surf[t,0,0].data,ds2.Vearth_surf[t,1,0].data],ds2.Vearth.data[t,1:]))

# # vertical interpolation to PALM grid
# Vearth_int = np.zeros((H_vel.shape[0],z.shape[0]))

# #define z from the creat_dynamic code
# #interpolate to z
# for t in range(0,Uearth_int.shape[0]): 
#     Vearth_int[t] = np.interp(z, H_vel[t], Vearth[t])

# #check the interpolated data
# #Uearth
# fig = plt.figure(figsize=(8,8))
# ax = fig.subplots(1,1)
# for t in range(0,Uearth_int.shape[0],2): 
#     ax.plot(Vearth_int[t],z,'-',label=str(t))
#     ax.set_ylabel('Geopotential height (m)',fontsize=18)
# fig.legend(fontsize=12)
# #plt.gca().invert_yaxis()
# plt.ylim((0., 2000.))
# #plt.xlim((275., 310.))
# #looks cool!   
}


#interpolate Uearth and Vearth to the PALM grid

#Uearth interpolation
Uearth_int = np.empty((ds2.Uearth.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,qv_int.shape[0]): 
    Uearth_int[t] = np.interp(z, H[t], ds2.Uearth[t])

#check the interpolated data

#u
fig = plt.figure(figsize=(8,8))
ax = fig.subplots(1,1)
for t in range(0,pt_int.shape[0],2): 
    ax.plot(Uearth_int[t],z,'-',label=str(t))
    ax.set_ylabel('Geopotential height (m)',fontsize=18)
fig.legend(fontsize=12)
#plt.gca().invert_yaxis()
plt.ylim((0., 2000.))
#plt.xlim((275., 310.))
#looks cool!   

#repeat the same for Vearth

#Vearth interpolation
Vearth_int = np.empty((ds2.Vearth.shape[0],z.shape[0]))

#define z from the creat_dynamic code
#interpolate to z
for t in range(0,qv_int.shape[0]): 
    Vearth_int[t] = np.interp(z, H[t], ds2.Vearth[t])

#check the interpolated data

#v
fig = plt.figure(figsize=(8,8))
ax = fig.subplots(1,1)
for t in range(0,pt_int.shape[0],2): 
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


#update the ches_dynamic file with the updated arrays
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\realistic_IOP03\trials\static_dynamic_coarse5_20hrs_spinup\input')

ches_dynamic = 'ches_scratch2_dynamic'
ches_dynamic_ds = xr.open_dataset(ches_dynamic)
ches_dynamic_ds.close()

#pt
ches_dynamic_ds.init_atmosphere_pt.values = pt_int[0]
ches_dynamic_ds.ls_forcing_left_pt.values = pt_int
ches_dynamic_ds.ls_forcing_right_pt.values = pt_int
ches_dynamic_ds.ls_forcing_south_pt.values = pt_int
ches_dynamic_ds.ls_forcing_north_pt.values = pt_int
ches_dynamic_ds.ls_forcing_top_pt.values = pt_int[:,239]

#qv
ches_dynamic_ds.init_atmosphere_qv.values = qv_int[0]
ches_dynamic_ds.ls_forcing_left_qv.values = qv_int
ches_dynamic_ds.ls_forcing_right_qv.values = qv_int
ches_dynamic_ds.ls_forcing_south_qv.values = qv_int
ches_dynamic_ds.ls_forcing_north_qv.values = qv_int
ches_dynamic_ds.ls_forcing_top_qv.values = qv_int[:,239]

#u
ches_dynamic_ds.init_atmosphere_u.values = Uearth_int[0]
ches_dynamic_ds.ls_forcing_left_u.values = Uearth_int
ches_dynamic_ds.ls_forcing_right_u.values = Uearth_int
ches_dynamic_ds.ls_forcing_south_u.values = Uearth_int
ches_dynamic_ds.ls_forcing_north_u.values = Uearth_int
ches_dynamic_ds.ls_forcing_top_u.values = Uearth_int[:,239]

#v
ches_dynamic_ds.init_atmosphere_v.values = Vearth_int[0]
ches_dynamic_ds.ls_forcing_left_v.values = Vearth_int
ches_dynamic_ds.ls_forcing_right_v.values = Vearth_int
ches_dynamic_ds.ls_forcing_south_v.values = Vearth_int
ches_dynamic_ds.ls_forcing_north_v.values = Vearth_int
ches_dynamic_ds.ls_forcing_top_v.values = Vearth_int[:,239]

#w
ches_dynamic_ds.init_atmosphere_w.values = WVEL_int[0]
ches_dynamic_ds.ls_forcing_left_w.values = WVEL_int
ches_dynamic_ds.ls_forcing_right_w.values = WVEL_int
ches_dynamic_ds.ls_forcing_south_w.values = WVEL_int
ches_dynamic_ds.ls_forcing_north_w.values = WVEL_int
ches_dynamic_ds.ls_forcing_top_w.values = WVEL_int[:,238]

ches_dynamic_ds.surface_forcing_surface_pressure.values = ds2.PRES_P0_L1_GLC0[:,0].values

ches_dynamic_ds.ls_forcing_ug.values = Ugearth_int
ches_dynamic_ds.ls_forcing_vg.values = Vgearth_int


#updating the lod field

#pt
ches_dynamic_ds.ls_forcing_left_pt.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_right_pt.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_south_pt.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_north_pt.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_top_pt.attrs['lod'] = np.int32(1)

#qv
ches_dynamic_ds.ls_forcing_left_qv.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_right_qv.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_south_qv.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_north_qv.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_top_qv.attrs['lod'] = np.int32(1)

#u
ches_dynamic_ds.ls_forcing_left_u.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_right_u.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_south_u.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_north_u.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_top_u.attrs['lod'] = np.int32(1)

#v
ches_dynamic_ds.ls_forcing_left_v.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_right_v.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_south_v.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_north_v.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_top_v.attrs['lod'] = np.int32(1)

#w
ches_dynamic_ds.ls_forcing_left_w.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_right_w.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_south_w.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_north_w.attrs['lod'] = np.int32(1)
ches_dynamic_ds.ls_forcing_top_w.attrs['lod'] = np.int32(1)

#ches_dynamic_ds.attrs['end_time'] =  str(times[-1]) + ' UTC'

#done!

#check plots before writing

#write the file
source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master\dynamic_files'

os.chdir(source)

#first updated version with the correct profiles 
#ches_dynamic_ds.to_netcdf('ches_90m_dynamic_1hour_52x48_L1_04252021.nc'), 


#ches_dynamic_ds.to_netcdf('ches_90m_dynamic_1hour_52x48_L1_05032021.nc')
#ches_dynamic_ds.close()

#update the dynamic file vertical extent and the forcing values on top
#changing the vertical extent from nz = 240 to nz = 168

#reading in the current data
ches_dynamic = 'ches_90m_dynamic_1hour_52x48_L1_05022021.nc'
ches_dynamic_ds = xr.open_dataset(ches_dynamic)
ches_dynamic_ds.close()

#subsetting for the nz
ds =  ches_dynamic_ds.isel(z=slice(0,180),zw=slice(0,179)).copy()

#updating the time series values at the domain top
ds.ls_forcing_top_pt.values = ds.ls_forcing_left_pt.isel(z=179).values
ds.ls_forcing_top_qv.values = ds.ls_forcing_left_qv.isel(z=179).values
ds.ls_forcing_top_u.values = ds.ls_forcing_left_u.isel(z=179).values
ds.ls_forcing_top_v.values = ds.ls_forcing_left_v.isel(z=179).values
ds.ls_forcing_top_w.values = ds.ls_forcing_left_w.isel(zw=178).values


ds.to_netcdf('ches_90m_dynamic_1hour_52x48_L1_05022021.nc')
ds.close()

ds = xr.open_dataset('ches_90m_dynamic_1hour_52x48_L1_05032021.nc')
ds.close()

#subsetting only for the morning CBL from 0900 to 1500 local
#subset the time array

time_ches_morning = ds.time[0:7].copy()
#subset data

ds_morning = ds.isel(time=slice(9,16),time_rad=slice(9,16)).copy()

ds_morning = ds_morning.assign_coords(time = time_ches_morning)
ds_morning = ds_morning.assign_coords(time_rad = time_ches_morning)


ds_morning.attrs['origin_time'] = str((datetime(2019, 9, 24,14,))) + 'UTC'
ds_morning.attrs['end_time'] = str((datetime(2019, 9, 24,20,))) + 'UTC'

ds_morning.to_netcdf('ches_90m_dynamic_1hour_52x48_L1_05032021_morning.nc')


#%%
source = r'C:\Users\Sreenath\Documents\palm\Cheyenne\ches_forcings\WRF4PALM-master\dynamic_files'

os.chdir(source)

ds1 = xr.open_dataset('ches_90m_dynamic1_1hour')
ds1 = ds1.mean('x').mean('y').mean('yv').mean('xu')

###testing the dynamic file with updated pt

#%% setting up the time delay
ds = xr.open_dataset('ches_90m_dynamic2_1hour')


delay = 18000000000000+3600000000000

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time')
ds = ds.rename({'time_updated': 'time'})

ds.time

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time_rad': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time_rad')
ds = ds.rename({'time_updated': 'time_rad'})

ds.time_rad
ds2 = ds.copy()
ds2 = ds2.mean('x').mean('y').mean('yv').mean('xu')
#%%
ds = xr.open_dataset('ches_90m_dynamic3_1hour')


delay = 39600000000000+3600000000000

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time')
ds = ds.rename({'time_updated': 'time'})

ds.time

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time_rad': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time_rad')
ds = ds.rename({'time_updated': 'time_rad'})

ds.time_rad
ds3 = ds.copy()
ds3 = ds3.mean('x').mean('y').mean('yv').mean('xu')


#ds3.time
#%%
ds = xr.open_dataset('ches_90m_dynamic4_1hour')


delay = 61200000000000+3600000000000

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time')
ds = ds.rename({'time_updated': 'time'})

ds.time

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time_rad': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time_rad')
ds = ds.rename({'time_updated': 'time_rad'})

ds.time_rad

ds4 = ds.copy()
ds4 = ds4.mean('x').mean('y').mean('yv').mean('xu')

#%%
ds = xr.open_dataset('ches_90m_dynamic5_1hour')


delay = 68400000000000+3600000000000

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time')
ds = ds.rename({'time_updated': 'time'})

ds.time

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.swap_dims({'time_rad': 'time_updated'})
#ds2 = ds2.swap_dims({'time_rad': 'time_updated'})

ds.coords['time_updated'] = ds.time_rad + np.array(delay).astype('timedelta64[ns]')
ds = ds.drop('time_rad')
ds = ds.rename({'time_updated': 'time_rad'})

ds.time_rad

ds5 = ds.copy()
ds5 = ds5.mean('x').mean('y').mean('yv').mean('xu')

#%%
# list_var = ches_dynamic_ds.var()

#combine time series
var_list = []
for var in ds1.data_vars:
    if "init" not in var:
        print(var)    
        var_list.append(var)
    


ches_dynamic_ds = xr.merge([ds1[var_list], ds2[var_list],ds3[var_list],ds4[var_list],ds5[var_list]],compat='no_conflicts')
#ches_dynamic_ds.init_atmosphere_pt

#ches_dynamic_ds.surface_forcing_surface_pressure.time


#pick initial profiles and soil data
for var in ds1.data_vars:
    if "init" in var:
        print(var)    
        ches_dynamic_ds[var] = ds1[var]


#%%

#test_ds = ches_dynamic_ds.mean('x').mean('y').mean('yv').mean('xu')

test_ds = ches_dynamic_ds.copy()

ches_dynamic_ds = xr.open_dataset('ches_90m_dynamic1_1hour')


test_ds['x'] = ches_dynamic_ds.x
test_ds['y'] = ches_dynamic_ds.y
test_ds['xu'] = ches_dynamic_ds.xu
test_ds['yv'] = ches_dynamic_ds.yv

res_origin = '90x90 m'

# test_ds.ls_forcing_left_w.attrs
# test_dynamic_ds.ls_forcing_left_w.attrs


test_ds['init_soil_m'] = ches_dynamic_ds['init_soil_m']
test_ds['init_soil_t'] = ches_dynamic_ds['init_soil_t']

# test_ds['rad_lw_in'] = ches_dynamic_ds['rad_lw_in']
# test_ds['rad_sw_in'] = ches_dynamic_ds['rad_sw_in'].attrs



# output boundary conditions to PALM input
# directions: 0 left, 1 right
#             0 south, 1 north


test_ds['rad_lw_in'].attrs = {'units':'W.m^2', 'lod':np.int32(1), 'source':'WRF', 'long_name':'Downward long-wave rad. flux'}
test_ds['rad_sw_in'].attrs = {'units':'W.m^2', 'lod':np.int32(1), 'source':'WRF', 'long_name':'Downward short-wave rad. flux'}


test_ds['init_atmosphere_pt'].attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}

test_ds['ls_forcing_left_pt'].attrs={'units':'K', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}

test_ds['ls_forcing_right_pt'].attrs={'units':'K','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}

test_ds['ls_forcing_south_pt'].attrs={'units':'K','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_north_pt'].attrs={'units':'K','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_top_pt'].attrs={'units':'K','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}

test_ds['init_atmosphere_qv'].attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_left_qv'].attrs={'units':'kg/kg','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_right_qv'].attrs={'units':'kg/kg','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_south_qv'].attrs={'units':'kg/kg','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_north_qv'].attrs={'units':'kg/kg','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_top_qv'].attrs={'units':'kg/kg','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}

test_ds['init_atmosphere_u'].attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_left_u'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_right_u'].attrs={'units':'m/s','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_south_u'].attrs={'units':'m/s','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_north_u'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_top_u'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}

test_ds['init_atmosphere_v'].attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_left_v'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_right_v'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_south_v'].attrs={'units':'m/s','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_north_v'].attrs={'units':'m/s','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_top_v'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}

test_ds['init_atmosphere_w'].attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_left_w'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_right_w'].attrs={'units':'m/s','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_south_w'].attrs={'units':'m/s','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_north_w'].attrs={'units':'m/s','lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_top_w'].attrs={'units':'m/s', 'lod':np.int32(1),'source':'WRF', 'res_origin':res_origin}

test_ds['surface_forcing_surface_pressure'].attrs={'units':'Pa', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin}


test_ds['ls_forcing_ug'].attrs={'units':'m/s', 'lod':np.int32(1),'long_name':'u wind component geostrophic', 'source':'WRF', 'res_origin':res_origin}
test_ds['ls_forcing_vg'].attrs={'units':'m/s', 'lod':np.int32(1),'long_name':'v wind component geostrophic', 'source':'WRF', 'res_origin':res_origin}




## add attributes
test_ds.attrs['description'] = f'Contains dynamic data from NCEP HRRR mesoscale. HRRR output file'
test_ds.attrs['author'] = 'Sreenath Paleri (paleri@wisc.edu)'
test_ds.attrs['history'] = 'Created at ' + time.ctime(time.time())
test_ds.attrs['source']= 'netCDF4 python'
test_ds.attrs['origin_lat'] = ds.origin_lat
test_ds.attrs['origin_lon'] = ds.origin_lon
test_ds.attrs['z'] = np.float(0)
test_ds.attrs['x'] = np.float(0)
test_ds.attrs['y'] = np.float(0)

test_ds.attrs['rotation_angle'] = np.float(0)
ds = xr.open_dataset('ches_90m_dynamic1_1hour')
test_ds.attrs['origin_time'] =  ds.origin_time
ds = xr.open_dataset('ches_90m_dynamic5_1hour')
test_ds.attrs['end_time'] =  ds.end_time

# test_ds.ls_forcing_vg

test_ds.to_netcdf('ches_90m_dynamic_1hour_52x48_L1_rad_04122021.nc')

#test_dynamic_ds.ls_forcing_left_w


#test_dynamic_ds.close()
#test_static_ds.close()

test_ds.close()
ches_dynamic_ds.close()

#%%
os.chdir(r'C:\Users\Sreenath\Documents\palm\Cheyenne\realistic_IOP03\INPUT')

ches_dynamic = 'ches_90m_dynamic_1hour_52x48_L1_rad_03192021.nc'

ches_dynamic_ds = xr.open_dataset(ches_dynamic)

ches_dynamic_ds.rad_lw_in.plot()

ds.drop('time')

ches_dynamic_ds.close()
del ches_dynamic_ds
