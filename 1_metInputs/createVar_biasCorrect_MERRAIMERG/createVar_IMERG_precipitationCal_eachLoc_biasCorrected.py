# createVar_IMERG_precipitationCal_eachLoc_biasCorrected.py
#
# This code is based on 2_biasCorrect_PRECTOT_eachLoc_GEOS.py.
#
# The IMERG precipitation values don't match the distribution of WMO observed precipitation.
# This script takes IMERG precipitationCal, subsets the data to each location in turn, 
# transforms the data  based on a fit of precipitationCal to observed precipitation (WMO), 
# and saves the results as netcdf files (to be used by createMetInput.py).

# import required packages
# ************************************************
import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress
import sys
import os

## Import constants
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants")
from CONST_FILE_DIR import (DIR_OBS, PATH_IMERG_PRECIP)


## Define parameters for describing the files
#  ***********************************************************************************
varName_obsData   = 'PRCP'    # name in WMO
varName_mainData  = 'precipitationCal' # name in IMERG

# City metadata
cities = {
    "Negombo": {"lat": 7.2008, "lon": 79.8737, "elev": 2., "obs_file": "Katunayake_WMO_2001-2020.csv"},
    "NuwaraEliya": {"lat": 6.9497, "lon": 80.7891, "elev": 1868., "obs_file": "NuwaraEliya_WMO_2001-2020.csv"},
    "Jaffna": {"lat": 9.6615, "lon": 80.0255, "elev": 5., "obs_file": "Jaffna_WMO_2009-2020.csv"}
}


## Define directories and filepaths
#  ***********************************************************************************
# Observation data used for bias correction (WMO)
fDir_obsData  = DIR_OBS

# Main data to bias correct (GEOS)
fPath_in_mainData = PATH_IMERG_PRECIP


## Define helper functions
#  ***********************************************************************************

# Load in daily WMO observation data (convert to total mm)
def load_obsData(filepath):

    df = pd.read_csv(filepath)
    df["DATE"] = pd.to_datetime(df["DATE"])
    df.set_index("DATE", inplace=True)
    
    for col in ["TEMP", "MAX", "MIN"]: # for daily average, max, and min temps
        df[col] = df[col].replace(9999.9, np.nan) # replace missing value for temperatures (9999.9) with nan
        df[col] = (df[col] - 32) * 5 / 9  # convert from F to C

    df[varName_obsData] = df[varName_obsData].replace(99.99, np.nan)  # replace missing value for precipitation (99.99) with nan
    df[varName_obsData] = df[varName_obsData] * 25.4                  # convert from in to mm (daily total)
    return df


# Load in hourly IMERG or GEOS data (convert to total mm)
def extract_mainData_prcp(ds, mainDataName, varName_mainData, lat, lon):

    if "bnds" in ds[varName_mainData].dims:
        ds = ds[varName_mainData].isel(bnds=0)  # pick the first bound
    else:
        ds = ds[varName_mainData]
    df = ds.sel(lat=lat, lon=lon, method='nearest').drop_vars(['lat', 'lon']).to_dataframe()
    if mainDataName == 'GEOS':
        varData = df[varName_mainData] * 3600  # convert from kg/m^2/s to total mm 
    else: # if IMERG
        varData = df[varName_mainData]         # already in total mm
    return varData


# Get subsets of observation and main data that correspond to the same times
def get_common_data(var_obsData, var_mainData):

    common_dates = var_obsData.index.intersection(var_mainData.index)
    common_dates = common_dates[(~var_obsData.loc[common_dates].isna()) & (~var_mainData.loc[common_dates].isna())] # remove dates corresponding to nan values
    var_obsData_common, var_mainData_common = var_obsData.loc[common_dates], var_mainData.loc[common_dates]
    return var_obsData_common, var_mainData_common


# Bias correct data via quantile mapping for each calendar month
def bias_correct_mainData_byMon(var_obsData_daily, var_mainData_hourly):

    ## Do quantile mapping on daily data
    var_mainData_daily = var_mainData_hourly.resample('D').sum()
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

    # Convert hourly mainData to the hourly value as a fraction of the daily total
    var_mainData_daily4eachHour = var_mainData_daily.reindex(var_mainData_hourly.index, method='ffill')
    var_mainData_hourly_asDailyFrac = var_mainData_hourly / var_mainData_daily4eachHour
    var_mainData_hourly_asDailyFrac[var_mainData_daily4eachHour == 0] = 0 # fix any divide by zeros that we just did

    # Multiply fractional hourly value with bias-corrected daily data
    var_mainData_daily4eachHour_bc = var_mainData_daily_bc.reindex(var_mainData_hourly.index, method='ffill')
    var_mainData_hourly_bc = var_mainData_hourly_asDailyFrac * var_mainData_daily4eachHour_bc

    return var_mainData_hourly_bc


if True:
    for city, cityInfo in cities.items():
        print(f'   {city}')

        # Load observation dataset
        fPath_in_obsData = os.path.join(fDir_obsData, cityInfo["obs_file"])
        df_obsData_daily   = load_obsData(fPath_in_obsData)  # in total mm
        var_obsData_daily = df_obsData_daily[varName_obsData]

        print(f'  Working on {varName_mainData}...')

        # Load main dataset
        ds_mainData = xr.open_dataset(fPath_in_mainData)
        var_mainData_hourly = extract_mainData_prcp(ds_mainData, 'GEOS', varName_mainData, cityInfo["lat"], cityInfo["lon"]) # function converts from kg/m^2/s to total mm

        # Define path for output data
        fDir_mainData = os.path.dirname(fPath_in_mainData)
        fName_mainData = os.path.basename(fPath_in_mainData)
        fName_mainData_root, fName_mainData_ext = os.path.splitext(fName_mainData)
        fPath_out_mainData = os.path.join(fDir_mainData, f'{fName_mainData_root}_{city}_biasCorrected{fName_mainData_ext}')

        # Do bias correction
        var_mainData_hourly_bc = bias_correct_mainData_byMon(var_obsData_daily, var_mainData_hourly)
        # No units conversion needed since we've been working in IMERG units (mm) the whole time

        # Save bias corrected data to new NetCDF files
        # First turn into xarray dataArray, then Dataset, then save that.
        da_mainData = ds_mainData[varName_mainData].sel(lat=cityInfo["lat"], lon=cityInfo["lon"], method='nearest')
        var_bc_da = xr.DataArray(var_mainData_hourly_bc, coords=da_mainData.coords, attrs=da_mainData.attrs)
        var_bc_ds = xr.Dataset({varName_mainData: var_bc_da})
        var_bc_ds.to_netcdf(fPath_out_mainData)

        # Close original datasets
        ds_mainData.close()

