# createVar_MERRA_T2M_TSOIL1_RH_NuwaraEliya_biasCorrected.py
#
# The MERRA 2-m and soil temperatures (T2M, TSOIL1) for Nuwara Eliya are systematically too high.
# This script takes MERRA T2M and TSOIL1, subsets the data to Nuwara Eliya, transforms them
# based on a linear regression of T2M to observed air temperature (WMO), and saves the results as
# netcdf files (to be used by createMetInput.py).
# We also recalculate and save bias-corrected RH as it is derived from T2M, but we didn't end up
# using the bias-corrected RH because its values were unrealistic and the original RH was more appropriate.
# Note: We transform TSOIL1 based on the linreg of T2M to observed air temperature because we don't
#       have observations of soil temperature to regress TSOIL1 against.

# import required packages
# ************************************************
import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import sys

# import required functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/1_metInputs")
from createVar_MERRA_RH import createVar_MERRA_RH

## Define directories and filepaths
#  ***********************************************************************************
fDir_MERRA = '/mnt/redwood/local_drive/chanud/RawData/MERRA_1980-2021_SL/'
fDir_WMO   = '/mnt/redwood/local_drive/chanud/RawData/observations/'
fPath_in_MERRA_T2M     = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.nc'
fPath_out_MERRA_T2M    = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.NuwaraEliya.biasCorrected.nc'
fPath_in_MERRA_TSOIL1  = fDir_MERRA + 'MERRA_hourly_LND/MERRA2_tavg1_2d_lnd_Nx.19800101_20211231.SUB.TSOIL1.nc'
fPath_out_MERRA_TSOIL1 = fDir_MERRA + 'MERRA_hourly_LND/MERRA2_tavg1_2d_lnd_Nx.19800101_20211231.SUB.TSOIL1.NuwaraEliya.biasCorrected.nc'
fPath_in_WMO    = fDir_WMO + 'NuwaraEliya_WMO_2001-2020.csv'


## Import and extract MERRA data
#  ***********************************************************************************

# Load the MERRA reanalysis data as xarray datasets
T2M_xrds = xr.open_dataset(fPath_in_MERRA_T2M)
TSOIL1_xrds = xr.open_dataset(fPath_in_MERRA_TSOIL1)

# Coordinates of Nuwara Eliya
lat = 6.9497  # deg N
lon = 80.7891 # deg E

# Extract data for the specific coordinate
T2M_xrda = T2M_xrds['T2M'].sel(lat=lat, lon=lon, method='nearest') # xarray data array
T2M_df = T2M_xrda.drop(['lat', 'lon']).to_dataframe() # pandas dataframe
T2M_s = T2M_df['T2M'] - 273.15 # pandas series; convert from K to C
TSOIL1_xrda = TSOIL1_xrds['TSOIL1'].sel(lat=lat, lon=lon, method='nearest') # xarray data array
TSOIL1_df = TSOIL1_xrda.drop(['lat', 'lon']).to_dataframe() # pandas dataframe
TSOIL1_s = TSOIL1_df['TSOIL1'] - 273.15 # pandas series; convert from K to C


## Do linear regression
#  ***********************************************************************************
#  Regressing MERRA-2 avg daily temps with WMO avg daily temps.

# Calculate daily average MERRA temp
T2Mavg = T2M_s.resample('D').mean()

# Load the WMO data
# Read the CSV file into a pandas DataFrame
data_WMO = pd.read_csv(fPath_in_WMO)
data_WMO.set_index(pd.to_datetime(data_WMO['DATE']), inplace=True)
# Replace the missing value for temperatures (9999.9) with nan
data_WMO['TEMP'].replace(9999.9, np.nan, inplace=True)
# Extract the daily average temperatures (to use for the MERRA bias correction)
Tavg_WMO = (data_WMO['TEMP'] - 32) * (5/9) # convert from F to C

# Subset MERRA-2 and WMO data to just the common dates for which they both have non-missing data.
common_dates_avg = T2Mavg.index.intersection(Tavg_WMO.index)
common_dates_avg = common_dates_avg[(~T2Mavg.loc[common_dates_avg].isna()) & (~Tavg_WMO.loc[common_dates_avg].isna())]
Tavg_WMO_subset = Tavg_WMO.loc[common_dates_avg]
T2Mavg_subset = T2Mavg.loc[common_dates_avg]

# Perform linear regression
slope_Tavg, int_Tavg, _, _, _ = linregress(Tavg_WMO_subset.values, T2Mavg_subset.values)


## Do bias correction
#  ***********************************************************************************
#  Bias correct the MERRA-2 hourly data (both T2M and TSOIL1) based on the regression.

# Create bias-corrected MERRA temperatures
T2M_bc_s = (T2M_s - int_Tavg)/slope_Tavg
T2M_bc_s = T2M_bc_s + 273.15  # convert from C to K (to match original MERRA units)
TSOIL1_bc_s = (TSOIL1_s - int_Tavg)/slope_Tavg
TSOIL1_bc_s = TSOIL1_bc_s + 273.15  # convert from C to K (to match original MERRA units)

# Save bias corrected data to new NetCDF files
# First turn into xarray dataArray, then Dataset, then save that.
T2M_bc_xrda = xr.DataArray(T2M_bc_s, coords=T2M_xrda.coords, attrs=T2M_xrda.attrs)
T2M_bc_xrds = xr.Dataset({'T2M': T2M_bc_xrda})
#T2M_bc_xrds.to_netcdf(fPath_out_MERRA_T2M)
TSOIL1_bc_xrda = xr.DataArray(TSOIL1_bc_s, coords=TSOIL1_xrda.coords, attrs=TSOIL1_xrda.attrs)
TSOIL1_bc_xrds = xr.Dataset({'TSOIL1': TSOIL1_bc_xrda})
#TSOIL1_bc_xrds.to_netcdf(fPath_out_MERRA_TSOIL1)

# Close the original datasets
T2M_xrds.close()
TSOIL1_xrds.close()




## Create RH variable based on bias-corrected T2M
#  ***********************************************************************************
# Note: the RH values this yields don't look good! They're regularly wayyy above 100%. 
# We should stick to using RH calculated based on the original T2M (rather than the 
# bias-corrected T2M).
# Why does RH look weird? My guess is that there's some systematic bias in QV2M and/or PS,
# and that these would need to be corrected before using them with T2M_bc to calculate RH.

## Define additional directories and filepaths for RH calculation
fPath_in_MERRA_QV2M = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.QV2M.nc'
fPath_in_MERRA_PS   = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.PS.nc'
fPath_out_MERRA_RH  = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.' + 'RH.NuwaraEliya.biasCorrected' + '.nc'

# Import and extract MERRA data
# Load the MERRA reanalysis data as xarray datasets
QV2M_xrds = xr.open_dataset(fPath_in_MERRA_QV2M)  # 2-m specific humidity
PS_xrds   = xr.open_dataset(fPath_in_MERRA_PS)    # surface pressure

# Extract data for the specific coordinate
QV2M_xrda =  QV2M_xrds['QV2M'].sel(lat=lat, lon=lon, method='nearest') # xarray data array
PS_xrda   =  PS_xrds['PS'].sel(lat=lat, lon=lon, method='nearest')     # xarray data array

createVar_MERRA_RH(QV2M_xrda, T2M_bc_xrda, PS_xrda, fPath_out_MERRA_RH)
