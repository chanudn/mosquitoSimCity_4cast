# createMetInput.py
#
# Ported to python from original ncl code (createMetInput.ncl).

# import required packages
# ************************************************
import xarray as xr
import pandas as pd
import sys
import re
import os

# import constants
# ************************************************
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants')
from CONST_FILE_DIR import (DIR_GEOS, FILE_GEOS_PREFIX,
                            PATH_MERRA_T, PATH_MERRA_TSOIL, PATH_MERRA_RH, PATH_MERRA_LC, PATH_MERRA_MC, PATH_MERRA_HC, 
                            PATH_IMERG_PRECIP, PATH_IMERG_PRECIP_NEGOMBO_BIASCORRECTED, PATH_IMERG_PRECIP_JAFFNA_BIASCORRECTED, PATH_IMERG_PRECIP_NUWARAELIYA_BIASCORRECTED, 
                            PATH_MERRA_T_NUWARAELIYA_BIASCORRECTED, PATH_MERRA_TSOIL_NUWARAELIYA_BIASCORRECTED)


# ********************************************************************************************************
#   Main functions (called by other files)
# ********************************************************************************************************


# create_met_input()
# ************************************************
# Create a WHATCH'EM meteorological input text file based on hourly meteorological data.
# Note: the format for the output filename is hardcoded here because it MUST match the 
# filename format in the WHATCH'EM code (run_container.pl). Careful if/when changing it!
#
# Arguments:
#   dataSrc    - string; 'MERRA_IMERG', 'GEOS'
#   loc        - string; 'Negombo', 'NuwaraEliya', 'Jaffna'
#   t_start    - int; must be in YYYYMMDDHH format
#   t_end      - int; must be in YYYYMMDDHH format
#   fpath_out  - string; filepath for the created file
#   lead_time  - int; indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#   ens_num    - int; indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
def create_met_input(dataSrc, loc, t_start, t_end, fdir_out, lead_time=None, ens_num=None):

    # Check args
    if dataSrc not in ['MERRA_IMERG', 'GEOS']:
        raise ValueError(f'Invalid argument: {dataSrc}. Please provide one of: MERRA_IMERG, GEOS')
    if loc not in ['Negombo', 'NuwaraEliya', 'Jaffna']:
        raise ValueError(f'Invalid argument: {loc}. Please provide one of: Negombo, NuwaraEliya, Jaffna')
    if not re.compile(r'^\d{10}$').match(str(t_start)):
        raise ValueError(f'Invalid argument: {t_start}. Please provide a number in YYYYMMDDHH format.')
    if not re.compile(r'^\d{10}$').match(str(t_end)):
        raise ValueError(f'Invalid argument: {t_start}. Please provide a number in YYYYMMDDHH format.')

    # Validate GEOS-specific args, then set string used in data source-specific directories/filenames
    if dataSrc == 'GEOS':
        # Error if GEOS-specific args not provided when they're expected
        if lead_time is None or ens_num is None:
            raise ValueError('For dataSrc="GEOS", both lead_time and ens must be provided.')
        dataSrc_extended = f'{dataSrc}_lead{lead_time}_ens{ens_num}'
    else:
        # Warn if GEOS-specific args were provided but won't be used
        if lead_time is not None or ens_num is not None:
            print('Warning: provided lead_time and ens are ignored when dataSrc is not "GEOS".')
        dataSrc_extended = dataSrc

    df = extract_data_to_df(dataSrc, loc, t_start, t_end, lead_time, ens_num) # subset data to desired location/timespan and put in one dataframe
    df = format_df_for_txtfile(df)                                            # format dataframe before saving to text file

    # Save dataframe to text file
    name_srcLocTime = f'{dataSrc_extended}-{loc}-{t_start}_{t_end}'
    fpath_out = os.path.join(fdir_out, f'container_input_{name_srcLocTime}.txt')
    df.to_csv(fpath_out, sep='\t', index=False)






# ********************************************************************************************************
#   Helper functions (only used within this file)
# ********************************************************************************************************


# extract_data_to_df()
# ************************************************
# Subset data to desired location/timespan and put it all in one dataframe.
# Importantly, also converts units to be WHATCH'EM-compatible.
#
# Arguments:
#   dataSrc    - string; 'MERRA_IMERG', 'GEOS'
#   loc        - string; 'Negombo', 'NuwaraEliya', 'Jaffna'
#   t_start     - int; must be in YYYYMMDDHH format
#   t_end      - int; must be in YYYYMMDDHH format
#   lead_time  - int; indicates GEOS lead time (optional, only use if dataSrc=GEOS)
#   ens_num    - int; indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
def extract_data_to_df(dataSrc, loc, t_start, t_end, lead_time=None, ens_num=None):

    if dataSrc not in ['MERRA_IMERG', 'GEOS']:
        raise ValueError(f'Invalid data source: {dataSrc}')
    if dataSrc == 'GEOS' and (lead_time is None or ens_num is None):
        raise ValueError("For GEOS data source, both 'lead_time' and 'ens_num' must be provided.")

    # Define data variable names and filepaths
    # MERRA/IMERG:
    #   T..................MERRA-2 T2M [K]
    #   RH.................MERRA-2 RH, calculated from MERRA-2 QV2M, T2M, PS [%]
    #   LC / MC / HC.......MERRA-2 CLDLOW / CLDMID / CLDHGH [fraction]
    #   TSOIL (50mm).......MERRA-2 TSOIL1 (0-0.0988m) [K]
    #   PRECIP.............IMERG precipitationCal [mm/hr]
    # GEOS:
    #   T..................GEOS T2M [K]
    #   RH.................GEOS RH, calculated from GEOS QV2M, T2M, PS [%]
    #   LC / MC / HC.......GEOS CLDTT (we use CLDTT for all three cloud levels) [fraction]
    #   TSOIL (50mm).......GEOS TS [K]
    #   PRECIP.............GEOS PRECTOT [kg/m^2/s]
    fname_GEOS_base = f'{FILE_GEOS_PREFIX}_lead{lead_time}_ens{ens_num}'
    data_dict = {
        'T':      {'MERRA_IMERG': {'varName': 'T2M', 
                                   'fpath'  : {'default': PATH_MERRA_T, 
                                               'NuwaraEliya_bc': PATH_MERRA_T_NUWARAELIYA_BIASCORRECTED}},
                   'GEOS':        {'varName': 'T2M', 
                                   'fpath'  : {'default': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_T2M.nc'),
                                               'NuwaraEliya_bc': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_T2M_NuwaraEliya_biasCorrected.nc')}}},
        'TSOIL':  {'MERRA_IMERG': {'varName': 'TSOIL1', 
                                   'fpath'  : {'default': PATH_MERRA_TSOIL, 
                                               'NuwaraEliya_bc': PATH_MERRA_TSOIL_NUWARAELIYA_BIASCORRECTED}},
                   'GEOS':        {'varName': 'TS', 
                                   'fpath'  : {'default': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_TS.nc'),
                                               'NuwaraEliya_bc': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_TS_NuwaraEliya_biasCorrected.nc')}}},
        'RH':     {'MERRA_IMERG': {'varName': 'RH', 
                                   'fpath'  : {'default': PATH_MERRA_RH}},
                   'GEOS':        {'varName': 'RH', 
                                   'fpath'  : {'default': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_RH.nc')}}},
        'LC':     {'MERRA_IMERG': {'varName': 'CLDLOW', 
                                   'fpath'  : {'default': PATH_MERRA_LC}},
                   'GEOS':        {'varName': 'CLDLO', 
                                   'fpath'  : {'Negombo_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDLO_Negombo.nc'),
                                               'Jaffna_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDLO_Jaffna.nc'),
                                               'NuwaraEliya_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDLO_NuwaraEliya.nc')}}},
        'MC':     {'MERRA_IMERG': {'varName': 'CLDMID', 
                                   'fpath'  : {'default': PATH_MERRA_MC}},
                   'GEOS':        {'varName': 'CLDMD', 
                                   'fpath'  : {'Negombo_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDMD_Negombo.nc'),
                                               'Jaffna_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDMD_Jaffna.nc'),
                                               'NuwaraEliya_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDMD_NuwaraEliya.nc')}}},
        'HC':     {'MERRA_IMERG': {'varName': 'CLDHGH', 
                                   'fpath'  : {'default': PATH_MERRA_HC}},
                   'GEOS':        {'varName': 'CLDHI', 
                                   'fpath'  : {'Negombo_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDHI_Negombo.nc'),
                                               'Jaffna_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDHI_Jaffna.nc'),
                                               'NuwaraEliya_est': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_CLDHI_NuwaraEliya.nc')}}},
        'PRECIP': {'MERRA_IMERG': {'varName': 'precipitationCal', 
                                   'fpath'  : {'default': PATH_IMERG_PRECIP,
                                               'Negombo_bc': PATH_IMERG_PRECIP_NEGOMBO_BIASCORRECTED,
                                               'Jaffna_bc': PATH_IMERG_PRECIP_JAFFNA_BIASCORRECTED,
                                               'NuwaraEliya_bc': PATH_IMERG_PRECIP_NUWARAELIYA_BIASCORRECTED}},
                   'GEOS':        {'varName': 'PRECTOT', 
                                   'fpath'  : {'default': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_PRECTOT.nc'),
                                               'Negombo_bc': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_PRECTOT_Negombo_biasCorrected.nc'),
                                               'Jaffna_bc': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_PRECTOT_Jaffna_biasCorrected.nc'),
                                               'NuwaraEliya_bc': os.path.join(DIR_GEOS, f'{fname_GEOS_base}_PRECTOT_NuwaraEliya_biasCorrected.nc')}}}
    }

    # Load data
    if loc == 'Negombo': # for Negombo, load in regular temperatures
        file_T      = xr.open_dataset(data_dict['T'][dataSrc]['fpath']['default'])
        file_TSOIL  = xr.open_dataset(data_dict['TSOIL'][dataSrc]['fpath']['default'])
    elif loc == 'Jaffna': # for Jaffna, load in regular temperatures
        file_T      = xr.open_dataset(data_dict['T'][dataSrc]['fpath']['default'])
        file_TSOIL  = xr.open_dataset(data_dict['TSOIL'][dataSrc]['fpath']['default'])
    elif loc == 'NuwaraEliya': # for Nuwara Eliya, load in bias-corrected temperatures
        file_T      = xr.open_dataset(data_dict['T'][dataSrc]['fpath'][f'{loc}_bc'])
        file_TSOIL  = xr.open_dataset(data_dict['TSOIL'][dataSrc]['fpath'][f'{loc}_bc'])
    file_PRECIP = xr.open_dataset(data_dict['PRECIP'][dataSrc]['fpath'][f'{loc}_bc']) # bias-corrected precip
    file_RH    = xr.open_dataset(data_dict['RH'][dataSrc]['fpath']['default']) # regular RH
    if dataSrc == 'GEOS': # estimated cloud cover
        file_LC    = xr.open_dataset(data_dict['LC'][dataSrc]['fpath'][f'{loc}_est'])
        file_MC    = xr.open_dataset(data_dict['MC'][dataSrc]['fpath'][f'{loc}_est'])
        file_HC    = xr.open_dataset(data_dict['HC'][dataSrc]['fpath'][f'{loc}_est'])
    elif dataSrc == 'MERRA_IMERG': # regular cloud cover
        file_LC    = xr.open_dataset(data_dict['LC'][dataSrc]['fpath']['default'])
        file_MC    = xr.open_dataset(data_dict['MC'][dataSrc]['fpath']['default'])
        file_HC    = xr.open_dataset(data_dict['HC'][dataSrc]['fpath']['default'])
    
    # Subset data based on location
    site_lat, site_lon = get_site_coordinates(loc)
    
    # the bias-corrected data is already subset by location, so we don't need to do anything to those
    if loc == 'Negombo':
        site_T      = file_T.sel(lat=site_lat, lon=site_lon, method='nearest')
        site_TSOIL  = file_TSOIL.sel(lat=site_lat, lon=site_lon, method='nearest')
    elif loc == 'Jaffna':
        site_T      = file_T.sel(lat=site_lat, lon=site_lon, method='nearest')
        site_TSOIL  = file_TSOIL.sel(lat=site_lat, lon=site_lon, method='nearest')
    elif loc == 'NuwaraEliya':
        site_T      = file_T
        site_TSOIL  = file_TSOIL
    site_PRECIP = file_PRECIP
    site_RH     = file_RH.sel(lat=site_lat, lon=site_lon, method='nearest')
    if dataSrc == 'GEOS': # already subset by location
        site_LC     = file_LC
        site_MC     = file_MC
        site_HC     = file_HC
    elif dataSrc == 'MERRA_IMERG': # needs to be subset by location
        site_LC     = file_LC.sel(lat=site_lat, lon=site_lon, method='nearest')
        site_MC     = file_MC.sel(lat=site_lat, lon=site_lon, method='nearest')
        site_HC     = file_HC.sel(lat=site_lat, lon=site_lon, method='nearest')


    # Subset data based on time span
    # Need to manually add an hour to the end time since the last index is not included
    site_time = slice(pd.to_datetime(t_start, format='%Y%m%d%H'), 
                      pd.to_datetime(t_end, format='%Y%m%d%H') + pd.Timedelta(hours=1))
    site_T      = site_T.sel(time=site_time)
    site_RH     = site_RH.sel(time=site_time)
    site_LC     = site_LC.sel(time=site_time)
    site_MC     = site_MC.sel(time=site_time)
    site_HC     = site_HC.sel(time=site_time)
    site_TSOIL  = site_TSOIL.sel(time=site_time)
    site_PRECIP = site_PRECIP.sel(time=site_time)

    # Convert units
    site_T     = site_T - 273.15         # convert from K to deg C
    site_TSOIL = site_TSOIL - 273.15     # convert from K to deg C
    if dataSrc == 'GEOS':
        site_PRECIP = site_PRECIP * 3600 # convert from kg/m^2/s to mm/hr

    # Convert data to pandas DataFrame
    df = pd.DataFrame({
        'YYYY': site_T['time'].dt.year.values,
        'MM': site_T['time'].dt.month.values,
        'DD': site_T['time'].dt.day.values,
        'HH': site_T['time'].dt.hour.values,
        'DOY': site_T['time'].dt.dayofyear.values,
        'T': site_T[data_dict['T'][dataSrc]['varName']].values,                 # [C] (converted from K)
        'RH': site_RH[data_dict['RH'][dataSrc]['varName']].values,              # [%]
        'PRECIP': site_PRECIP[data_dict['PRECIP'][dataSrc]['varName']].values,  # [mm/hr] (if GEOS, converted from kg/m^2/s)
        'LC': site_LC[data_dict['LC'][dataSrc]['varName']].values,              # [fraction]
        'MC': site_MC[data_dict['MC'][dataSrc]['varName']].values,              # [fraction]
        'HC': site_HC[data_dict['HC'][dataSrc]['varName']].values,              # [fraction]
        'TSOIL': site_TSOIL[data_dict['TSOIL'][dataSrc]['varName']].values      # [C] (converted from K)
    })

    return df





# get_site_coordinates()
# ************************************************
#
# Get location coordinates based on location name.
#
# Arguments:
#   loc - string; 'Negombo', 'NuwaraEliya', 'Jaffna'
#                 (can handle 'Colombo', 'Kandy', 'Trincomalee', but those are currently unused in the overall code)
def get_site_coordinates(loc):
    coordinates = {    # [deg N, deg E]
        'Negombo':     (7.2008, 79.8737),
        'Colombo':     (6.9271, 79.8612),
        'Kandy':       (7.2906, 80.6337),
        'NuwaraEliya': (6.9497, 80.7891),
        'Trincomalee': (8.5874, 81.2152),
        'Jaffna':      (9.6615, 80.0255)
    }
    return coordinates[loc]




# format_df_for_txtfile()
# ************************************************
# Format dataframe before saving to text file.
# Date columns are padded with leading zeroes (e.g., the first month is '01' rather than '1').
# Other columns are rounded to either two or four decimal places.
#
# Arguments:
#   df_init  - input dataframe of hourly data
def format_df_for_txtfile(df_init):

    df = df_init.copy(deep=True)

    cols_zfill2 = ['MM', 'DD', 'HH']
    df[cols_zfill2] = df[cols_zfill2].astype(str).apply(lambda x: x.str.zfill(2))

    cols_zfill3 = ['DOY']
    df[cols_zfill3] = df[cols_zfill3].astype(str).apply(lambda x: x.str.zfill(3))

    cols_zfill4 = ['YYYY']
    df[cols_zfill4] = df[cols_zfill4].astype(str).apply(lambda x: x.str.zfill(4))

    cols_round2 = ['T', 'RH', 'PRECIP', 'TSOIL']
    df[cols_round2] = df[cols_round2].round(2).apply(lambda x: x.apply(lambda y: '{:.2f}'.format(y)))

    cols_round4 = ['LC', 'MC', 'HC']
    df[cols_round4] = df[cols_round4].round(4).apply(lambda x: x.apply(lambda y: '{:.4f}'.format(y)))

    return df




# calc_climatology()
# ************************************************
#
# ***********************NOT USED IN A LONG TIME; PROBABLY NEEDS FIXING!****************************************************
#
# Calculate climatology (averaging hourly data across multiple years).
# Sets the YYYY column to a missing value of '9999'.
#
# Arguments:
#   df_init  - input dataframe of hourly data
def calc_climatology(df_init):

    df = df_init.copy(deep=True)

    # Compute climatology by grouping by MMDDHH
    df['MMDDHH'] = df['MM'].astype(str).str.zfill(2) + df['DD'].astype(str).str.zfill(2) + df['HH'].astype(str).str.zfill(2)
    df = df.drop(['YYYY', 'MM', 'DD', 'HH', 'DOY'], axis=1)  # don't want these for climatology calc
    clim_df_B = df.groupby('MMDDHH').mean().reset_index()    # reset_index ensures MMDDHH remains a column of the df

    # Recreate date/time columns
    clim_df_A = pd.DataFrame({
        'YYYY': '9999', # climatology doesn't have a year; use a placeholder value with the same number of digits
        'MM': clim_df_B['MMDDHH'].str[:2],
        'DD': clim_df_B['MMDDHH'].str[2:4],
        'HH': clim_df_B['MMDDHH'].str[4:]
    })

    # Generate DOY column
    sampleLeapYr = 2000 # temp year value (a leap year) for calculating DOY
    date_str = str(sampleLeapYr) + '-' + clim_df_A['MM'].astype(str) + '-' + clim_df_A['DD'].astype(str)
    clim_df_A['DOY'] = pd.to_datetime(date_str, format='%Y-%m-%d').dt.dayofyear

    # Format climatology data into one dataframe
    clim_df = pd.concat([clim_df_A, clim_df_B], axis=1)
    clim_df = clim_df.drop('MMDDHH', axis=1) # drop unwanted columns

    return clim_df




