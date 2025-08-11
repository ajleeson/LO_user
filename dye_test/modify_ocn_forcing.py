'''
Modified Jilian's tracer nudge code to modify ocean forcing files
'''


# modify nudge_coef for tracer
import xarray as xr

# old forcing files from Parker's apogee account
print('Getting original ocean forcing')
ds_ini = xr.open_dataset('../../LO_output/forcing/cas7/f2012.10.07/ocnG00/ocean_ini_OG.nc')
ds_bry = xr.open_dataset('../../LO_output/forcing/cas7/f2012.10.07/ocnG00/ocean_bry_OG.nc')
ds_clm = xr.open_dataset('../../LO_output/forcing/cas7/f2012.10.07/ocnG00/ocean_clm_OG.nc')

################################################################
# add dye to the different datasets, using a copy of temp as a reference

# Ocean initial conditions
ds_ini['dye_01'] = xr.zeros_like(ds_ini['temp'])
ds_ini['dye_01'].attrs['long_name'] = 'Passive tracer dye concentration'
ds_ini['dye_01'].attrs['units'] = 'mg m-3' 

# Ocean boundary conditions
# North
dye_time_north = xr.zeros_like(ds_bry['temp_north']).rename({'temp_time': 'dye_time'})
ds_bry['dye_01_north'] = dye_time_north
ds_bry['dye_01_north'].attrs['long_name'] = 'Passive tracer dye northern boundary conditions'
ds_bry['dye_01_north'].attrs['units'] = 'mg m-3' 
# South
dye_time_south = xr.zeros_like(ds_bry['temp_south']).rename({'temp_time': 'dye_time'})
ds_bry['dye_01_south'] = dye_time_south
ds_bry['dye_01_south'].attrs['long_name'] = 'Passive tracer dye southern boundary conditions'
ds_bry['dye_01_south'].attrs['units'] = 'mg m-3' 
# East
dye_time_east = xr.zeros_like(ds_bry['temp_east']).rename({'temp_time': 'dye_time'})
ds_bry['dye_01_east'] = dye_time_east
ds_bry['dye_01_east'].attrs['long_name'] = 'Passive tracer dye eastern boundary conditions'
ds_bry['dye_01_east'].attrs['units'] = 'mg m-3' 
# West
dye_time_west = xr.zeros_like(ds_bry['temp_west']).rename({'temp_time': 'dye_time'})
ds_bry['dye_01_west'] = dye_time_west
ds_bry['dye_01_west'].attrs['long_name'] = 'Passive tracer dye western boundary conditions'
ds_bry['dye_01_west'].attrs['units'] = 'mg m-3' 
# Add dye_time coordinate
ds_bry['dye_time'].attrs = ds_bry['temp_time'].attrs.copy()

# Ocean climatology
dye_time = xr.zeros_like(ds_clm['temp']).rename({'temp_time': 'dye_time'})
ds_clm['dye_01'] = dye_time
ds_clm['dye_01'].attrs['long_name'] = 'Passive tracer dye concentration'
ds_clm['dye_01'].attrs['units'] = 'mg m-3' 
# Add dye_time coordinate
ds_clm['dye_time'].attrs = ds_clm['temp_time'].attrs.copy()


################################################################
# Save .nc files
print('Saving new files...')
ds_ini.to_netcdf(path='../../LO_output/forcing/cas7/f2012.10.07/ocnG00/ocean_ini.nc')
ds_bry.to_netcdf(path='../../LO_output/forcing/cas7/f2012.10.07/ocnG00/ocean_bry.nc')
ds_clm.to_netcdf(path='../../LO_output/forcing/cas7/f2012.10.07/ocnG00/ocean_clm.nc')

