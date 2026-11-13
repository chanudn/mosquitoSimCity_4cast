# createVar_MERRA_T2M_TSOIL1_RH_biasCorrected.py
#
# The MERRA 2-m and soil temperatures (T2M, TSOIL1) for Nuwara Eliya are systematically too high.
# This script takes MERRA T2M and TSOIL1, subsets the data to Nuwara Eliya, transforms them
# based on a linear regression of T2M to observed air temperature (GSOD), and saves the results as
# netcdf files (to be used by createMetInput.py).
# We also recalculate and save bias-corrected RH as it is derived from T2M, but we didn't end up
# using the bias-corrected RH because its values were unrealistic and the original RH was more appropriate.
# Note: We transform TSOIL1 based on the linreg of T2M to observed air temperature because we don't
#       have observations of soil temperature to regress TSOIL1 against.
#
# Revision history:
# 2025-09-14, CNY - Changed bias correction to match new method used for GEOS data. Now (a) the linear regression
#                   happens separately for each calendar month and (b) we convert back to hourly data by applying
#                   applying the linear regression to just the daily means and then adding on hourly anomalies
#                   (rather than applying linreg to all of the hourly data).
#                   Deleted unused script for generating RH based on bias-corrected T2M.
# 2025-10-21, CNY - Changed bias correction to match new method used for GEOS data. Now the bias correction is
#                   based on quantile mapping for each calendar month.

# import required packages
# ************************************************
import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import sys
import os


# import required functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/1_metInputs")
from createVar_MERRA_RH import createVar_MERRA_RH


## Define parameters for describing the files
#  ***********************************************************************************
varNames   = ['T2M', 'TSOIL1'] # names in main dataset (MERRA)

# City metadata
cities = {
    "Negombo": {"lat": 7.2008, "lon": 79.8737, "elev": 2., "obs_file": "Katunayake_WMO_2001-2020.csv"},
    "NuwaraEliya": {"lat": 6.9497, "lon": 80.7891, "elev": 1868., "obs_file": "NuwaraEliya_WMO_2001-2020.csv"},
    "Jaffna": {"lat": 9.6615, "lon": 80.0255, "elev": 5., "obs_file": "Jaffna_WMO_2009-2020.csv"}
}


## Define directories and filepaths
#  ***********************************************************************************
fDir_MERRA = '/mnt/redwood/local_drive/chanud/RawData/MERRA_1980-2021_SL/'
fDir_GSOD   = '/mnt/redwood/local_drive/chanud/RawData/observations/'
fPath_in_MERRA_T2M     = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.nc'
fPath_in_MERRA_TSOIL1  = fDir_MERRA + 'MERRA_hourly_LND/MERRA2_tavg1_2d_lnd_Nx.19800101_20211231.SUB.TSOIL1.nc'


## Define helper functions
#  ***********************************************************************************

# Load in daily GSOD observation data (convert temp to C, precip to total mm)
def load_obsData(filepath):

    df = pd.read_csv(filepath)
    df["DATE"] = pd.to_datetime(df["DATE"])
    df.set_index("DATE", inplace=True)
    
    for col in ["TEMP", "MAX", "MIN"]: # for daily average, max, and min temps
        df[col] = df[col].replace(9999.9, np.nan) # replace missing value for temperatures (9999.9) with nan
        df[col] = (df[col] - 32) * 5 / 9  # convert from F to C

    df["PRCP"] = df["PRCP"].replace(99.99, np.nan)  # replace missing value for precipitation (99.99) with nan
    df["PRCP"] = df["PRCP"] * 25.4                  # convert from in to mm (daily total)
    return df


# Load in main data (convert temp to C)
def extract_mainData_temp(ds, varName, lat, lon):
    df = ds.sel(lat=lat, lon=lon, method='nearest').drop_vars(['lat', 'lon']).to_dataframe()
    return df[varName] - 273.15  # convert from K to C


# Get common data between obs and main data
def get_common_data(T_obsData, T_mainData):
    common_dates = T_obsData.index.intersection(T_mainData.index)
    common_dates = common_dates[(~T_mainData.loc[common_dates].isna()) & (~T_obsData.loc[common_dates].isna())] # remove dates corresponding to nan values
    T_obsData_common, T_mainData_common = T_obsData.loc[common_dates], T_mainData.loc[common_dates]
    return T_obsData_common, T_mainData_common


# Bias correct data via quantile mapping for each calendar month
def bias_correct_mainData_byMon(var_obsData_daily, var_mainData_hourly):

    ## Do quantile mapping on daily data for each month
    var_mainData_daily = var_mainData_hourly.resample('D').mean()
    var_mainData_daily_bc = pd.Series(index=var_mainData_daily.index, dtype='float64')  # create empty series to store bias-corrected values
    for month in range(1, 12+1):
        # Extract common data for the current month
        obs_mon      = var_obsData_daily[var_obsData_daily.index.month == month]
        mainData_mon = var_mainData_daily[var_mainData_daily.index.month == month]
        obs_mon_common, mainData_mon_common = get_common_data(obs_mon, mainData_mon)

        # Skip if no common data
        if len(obs_mon_common) == 0 or len(mainData_mon_common) == 0:
            print(f'Error: no common data for month {month}!')
            continue

        # Perform quantile mapping for the current month
        quantiles = np.linspace(0, 1, 100+1)
        obs_q = np.quantile(obs_mon_common, quantiles)
        mainData_q = np.quantile(mainData_mon_common, quantiles)
        mainData_mon_bc = pd.Series(np.interp(mainData_mon, mainData_q, obs_q, left=obs_q[0], right=obs_q[-1]), index=mainData_mon.index)
        var_mainData_daily_bc.loc[mainData_mon_bc.index] = mainData_mon_bc

    ## Convert bias-corrected data back from daily to hourly

    # Convert hourly mainData to the hourly value as an anomaly relative to the daily mean
    var_mainData_daily4eachHour = var_mainData_daily.reindex(var_mainData_hourly.index, method='ffill')
    var_mainData_hourly_asDailyAnom = var_mainData_hourly - var_mainData_daily4eachHour

    # Add anomalous hourly value with bias-corrected daily data
    var_mainData_daily4eachHour_bc = var_mainData_daily_bc.reindex(var_mainData_hourly.index, method='ffill')
    var_mainData_hourly_bc = var_mainData_hourly_asDailyAnom + var_mainData_daily4eachHour_bc

    return var_mainData_hourly_bc




# DEPRECATED
# This version of bias correction take the hourly data for each calendar month, converts to daily data,
# does a linear regression against the observations, applies this linear regression to the daily data, 
# and reconstructs the hourly data as the bias-corrected daily data + the hourly anomalies.
#
#def bias_correct_mainData_linReg_byMon(var_obsData_daily, var_mainData_hourly):
#
#    # Do linear regression on daily data for each month
#    var_mainData_daily = var_mainData_hourly.resample('D').mean()
#    var_mainData_daily_bc = pd.Series(index=var_mainData_daily.index, dtype='float64')  # create empty series to store bias-corrected values
#    for month in range(1, 12+1):
#        # Extract common data for the current month
#        obs_mon      = var_obsData_daily[var_obsData_daily.index.month == month]
#        mainData_mon = var_mainData_daily[var_mainData_daily.index.month == month]
#        obs_mon_common, mainData_mon_common = get_common_data(obs_mon, mainData_mon)
#
#        # Skip if no common data
#        if len(obs_mon_common) == 0 or len(mainData_mon_common) == 0:
#            print(f'Error: no common data for month {month}!')
#            continue
#
#        # Perform linear regression for the current month
#        slope, intercept, r, *_ = linregress(mainData_mon_common, obs_mon_common)
#        mainData_mon_bc = slope*mainData_mon + intercept
#        var_mainData_daily_bc.loc[mainData_mon_bc.index] = mainData_mon_bc
#    
#    ## Convert bias-corrected data back from daily to hourly
#    
#    # Convert hourly mainData to the hourly value as an anomaly relative to the daily value
#    var_mainData_daily4eachHour = var_mainData_daily.reindex(var_mainData_hourly.index, method='ffill')
#    var_mainData_hourly_asDailyAnom = var_mainData_hourly - var_mainData_daily4eachHour
#
#    # Add anomalous hourly value with bias-corrected daily data
#    var_mainData_daily4eachHour_bc = var_mainData_daily_bc.reindex(var_mainData_hourly.index, method='ffill')
#    var_mainData_hourly_bc = var_mainData_hourly_asDailyAnom + var_mainData_daily4eachHour_bc
#
#    return var_mainData_hourly_bc




## Main script
#  ***********************************************************************************

for city, cityInfo in cities.items():
    print(f'{city}')

    # Load observation dataset (we use GSOD air temperature for bias correcting all temperatures in the main dataset)
    fPath_in_obsData = os.path.join(fDir_GSOD, cityInfo["obs_file"])
    df_obsData_daily   = load_obsData(fPath_in_obsData) # has daily Tavg, Tmin, Tmax
    var_obsData_daily = df_obsData_daily["TEMP"] # get just Tavg, which we use for all the bias correction

    for varName in varNames:
        print(f'  {varName}')

        # Set data filepaths based on variable
        if varName == 'T2M':
            fPath_in_mainData  = fPath_in_MERRA_T2M
            fPath_out_mainData = os.path.join(fDir_MERRA, f'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.{city}.biasCorrected.nc')
        elif varName == 'TSOIL1':
            fPath_in_mainData  = fPath_in_MERRA_TSOIL1 
            fPath_out_mainData = os.path.join(fDir_MERRA, f'MERRA_hourly_LND/MERRA2_tavg1_2d_lnd_Nx.19800101_20211231.SUB.TSOIL1.{city}.biasCorrected.nc')
        else:
            print(f'Error: unknown variable name {varName}')
            sys.exit()

        # Load main dataset
        ds_mainData = xr.open_dataset(fPath_in_mainData)
        var_mainData_hourly = extract_mainData_temp(ds_mainData, varName, cityInfo["lat"], cityInfo["lon"])

        # Do bias correction
        var_mainData_hourly_bc = bias_correct_mainData_byMon(var_obsData_daily, var_mainData_hourly)
        var_mainData_hourly_bc = var_mainData_hourly_bc + 273.15 # convert from C to K (to match original units)

        # Save bias corrected data to new NetCDF files
        # First turn into xarray dataArray, then Dataset, then save that.
        da_mainData = ds_mainData[varName].sel(lat=cityInfo["lat"], lon=cityInfo["lon"], method='nearest')
        var_bc_da = xr.DataArray(var_mainData_hourly_bc, coords=da_mainData.coords, attrs=da_mainData.attrs)
        var_bc_ds = xr.Dataset({varName: var_bc_da})
        var_bc_ds.to_netcdf(fPath_out_mainData)

        # Close original datasets
        ds_mainData.close()
