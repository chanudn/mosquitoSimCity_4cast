#*******************************************************
# fns_analysis.py
#*******************************************************
#
# The functions in this file process the data from the
# modeling pipeline. In some cases the resulting data
# is saved to file. In other cases they're returned by
# the function (e.g., as a precursor to plotting the data).
#
#
# Revision history
# 2024-10-19, CNY: file created
# 2025-05-26, CNY: added fileDir argument for custom save location
# 2025-05-31, CNY: added stats calculations to compare GEOS vs MERRA/IMERG
#
# Functions:
#   calc_tercile_crosstab_splitByMon()
#   calc_R_RMSE_bias()
#   compute_rps()
#   compute_rpss()
#   calc_table_comparisonStats()
#   get_table_comparisonStats()



# import required packages
# ************************************************
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.stats import pearsonr
from sklearn.metrics import root_mean_squared_error
import calendar


# import constants and custom functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants") # directory of constants file
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/4_analysis") # directory of helper functions file
from fns_helper import categorizeSeries_LowMedHigh_byMon




def calc_tercile_crosstab_splitByMon(df_in, varX_colName, varY_colName, fileDir, fileName):

   df = df_in.copy(deep=True)  # create a temporary copy of the df so that we don't alter the original
   
   # Categorize data into Low, Medium, High, on a monthly basis (and then recombine all monthly results into one series)
   varX_terciles = categorizeSeries_LowMedHigh_byMon(df[varX_colName])    
   varY_terciles = categorizeSeries_LowMedHigh_byMon(df[varY_colName])

   # If categorization failed for either variable, we're not going to do the calculation.
   if varX_terciles.isnull().all():
      print(f'Error: categorization failed for variable {varX_colName}. Skipping calculation.')
      return
   if varY_terciles.isnull().all():
      print(f'Error: categorization failed for variable {varY_colName}. Skipping calculation.')
      return

   # Add tercile columns to the DataFrame
   varXtercile_colName = varX_colName + '_tercile'
   varYtercile_colName = varY_colName + '_tercile'
   df[varXtercile_colName] = varX_terciles
   df[varYtercile_colName] = varY_terciles

   # Step 2: Create a crosstab to calculate the counts of varY terciles within each varX tercile
   crosstab = pd.crosstab(df[varXtercile_colName], df[varYtercile_colName])

   # Step 3: Normalize the crosstab to get percentages
   crosstab_normalized = crosstab.div(crosstab.sum(axis=1), axis=0) * 100

   # Save to csv files
   os.makedirs(fileDir, exist_ok=True) # make save directory if it doesn't exist
   crosstab.to_csv(os.path.join(fileDir, f'{fileName}.csv'))
   crosstab_normalized.to_csv(os.path.join(fileDir, f'{fileName}_norm.csv'))



# Calculate R (and its associated p-value), RMSE, and bias between two datasets
def calc_R_RMSE_bias(main_vals, ref_vals):

    # Check for matching number of data points and non-NaNs and then compute R, RMSE, and bias
    if len(main_vals) == len(ref_vals) and len(main_vals) > 1:
        r, p_r = pearsonr(ref_vals, main_vals)
        rmse   = root_mean_squared_error(ref_vals, main_vals)
        bias   = np.mean(main_vals - ref_vals) # mean bias of GEOS ensemble mean relative to MERRA/IMERG
    else:
        r    = np.nan
        p_r  = np.nan
        rmse = np.nan
        bias = np.nan

    return r, p_r, rmse, bias



# Compute RPS for forecast-observation pair
#   forecast_probs   - probabilities for low, med, high terciles
#   observed_tercile - tercile of observation (0=low, 1=med, or 2=high)
# This function is only used within compute_rpss().
def compute_rps(forecast_probs, observed_tercile):

    obs_cdf = np.zeros(3)
    obs_cdf[observed_tercile:] = 1

    forecast_cdf = np.cumsum(forecast_probs)

    rps = np.sum((forecast_cdf - obs_cdf) ** 2)
    return rps



# Compute RPSS for forecast-observation pairs using defined bins (e.g., terciles)
#   observed_vals  - observations, one per year (in this case 20 years)
#   ens_vals_by_yr - list of N lists, each with k ens members (in this case N=20 years, k=5 ensemble members)
def compute_rpss(observed_vals, ens_vals_by_yr):

    observed_vals = np.array(observed_vals)
    n_years = len(observed_vals)

    # Step 1: Define tercile bins separately for observed values and for ensemble values
    tercile_bins_obs = np.percentile(observed_vals, [33.33, 66.67])
    ens_vals_by_yr_flattened = [v for sublist in ens_vals_by_yr for v in sublist]
    tercile_bins_ens = np.percentile(ens_vals_by_yr_flattened,  [33.33, 66.67])

    # Map observations to tercile category: 0 = low, 1 = mid, 2 = high
    obs_terciles = np.digitize(observed_vals, tercile_bins_obs)

    rps_all = []       # to store the RPS results for each year (forecasts based on ensemble members)
    rps_clim_all = []  # to store the RPSclim results for each year (forecasts based on climatology)
    for year in range(n_years):
        ens_vals = ens_vals_by_yr[year]
        obs_tercile = obs_terciles[year]

        # Count how many ensemble members fall into each tercile
        ens_terciles = np.digitize(ens_vals, tercile_bins_ens)
        probs = np.array([np.sum(ens_terciles == 0), np.sum(ens_terciles == 1), np.sum(ens_terciles == 2)]) / len(ens_vals)

        rps = compute_rps(probs, obs_tercile)
        rps_all.append(rps)

        # Climatology forecast (uniform 1/3 probability for each tercile)
        clim_probs = np.array([1/3, 1/3, 1/3])
        rps_clim = compute_rps(clim_probs, obs_tercile)
        rps_clim_all.append(rps_clim)
    
    mean_rps = np.mean(rps_all)
    mean_rps_clim = np.mean(rps_clim_all)

    # Final RPSS
    if mean_rps_clim == 0:
        return np.nan
    else:
        rpss = 1 - (mean_rps / mean_rps_clim)
        return rpss



# Calculate and compile statistics (R, RMSE, bias, RPSS) comparing GEOS and MERRA/IEMRG data for one variable
def calc_table_comparisonStats(dfs_GEOS_ens_in, df_GEOS_ensMean_in, df_MERRAIMERG_in, var_colName, fileDir, fileName):

   ## Copy all input dfs so we don't modify them
   dfs_GEOS_ens    = [df.copy() for df in dfs_GEOS_ens_in]
   df_GEOS_ensMean = df_GEOS_ensMean_in.copy()
   df_MERRAIMERG   = df_MERRAIMERG_in.copy()

   # Add month and year columns to support calculations
   for df in dfs_GEOS_ens:
       df['Month'] = df.index.month
       df['Year']  = df.index.year
   df_GEOS_ensMean['Month'] = df_GEOS_ensMean.index.month
   df_GEOS_ensMean['Year'] = df_GEOS_ensMean.index.year
   df_MERRAIMERG['Month'] = df_MERRAIMERG.index.month
   df_MERRAIMERG['Year'] = df_MERRAIMERG.index.year

   # Prepare a structure to store all stats
   stats_by_month = {'RPSS': [], 'R': [], 'p_R': [], 'RMSE': [], 'Bias': []}
   for month in range(1, 12+1):
      ## Get data for the month and variable in question
      GEOS_ensMean_vals_oneMon = df_GEOS_ensMean[df_GEOS_ensMean['Month'] == month][var_colName].values # GEOS ensemble mean
      MERRAIMERG_vals_oneMon = df_MERRAIMERG[df_MERRAIMERG['Month'] == month][var_colName].values       # MERRA/IMERG
      
      # Gather GEOS ensemble members for RPSS computation
      ens_vals_byYr_oneMon = [] # a list of lists (a list of 5 ensemble members for each year)         
      for i in range(len(MERRAIMERG_vals_oneMon)): # for each year
         ens_vals_i = [] # a list of 5 ensemble members for the given year
         for df_ens in dfs_GEOS_ens: # for each ensemble member
            ens_vals_month = df_ens[df_ens['Month'] == month][var_colName].tolist()
            if len(ens_vals_month) > i:
               ens_vals_i.append(ens_vals_month[i]) # append the value corresponding to the ith year
         ens_vals_byYr_oneMon.append(ens_vals_i)

      ## Calculate R (and its associated p-value), RMSE, and bias between GEOS ensemble mean and MERRA/IMERG data
      r, p_r, rmse, bias = calc_R_RMSE_bias(GEOS_ensMean_vals_oneMon, MERRAIMERG_vals_oneMon)

      ## Compute RPSS between GEOS ensemble members and MERRA/IMERG data
      if all(len(ens) == len(dfs_GEOS_ens) for ens in ens_vals_byYr_oneMon): # check if all years have same number of ensemble values
         rpss = compute_rpss(MERRAIMERG_vals_oneMon, ens_vals_byYr_oneMon)
      else:
         rpss = np.nan

      # Store stats
      stats_by_month['RPSS'].append(rpss)
      stats_by_month['R'].append(r)
      stats_by_month['p_R'].append(p_r)
      stats_by_month['RMSE'].append(rmse)
      stats_by_month['Bias'].append(bias)

   # Create dataframe with stats as rows, months as columns
   month_names = calendar.month_abbr[1:]  # ['Jan', 'Feb', ..., 'Dec']
   comparison_table = pd.DataFrame(stats_by_month, index=month_names).T
   
   # Save to csv files
   os.makedirs(fileDir, exist_ok=True) # make save directory if it doesn't exist
   comparison_table.to_csv(os.path.join(fileDir, fileName))



# Retrieves previously calculated statistics from saved csv file
# Note: tspan must be list of two dates ['YYYY-MM-DD', 'YYYY-MM-DD']
def get_table_comparisonStats(dataSrc_main, dataSrc_ref, lead_time, varName, loc, dengueDataYrs_noOutbreaks_opt, survCurve, doAnom_opt):

   if doAnom_opt:
      fileStr_doAnom = 'Anom_'
   else:
      fileStr_doAnom = ''
   dataSrcs_extended = f'{fileStr_doAnom}{dataSrc_main}_lead{lead_time}_VS_{dataSrc_ref}'

   if survCurve == 'Eisen':
       fileStr_survCurve = '_eisenQuadFit'
   else: # if survCurve is 'Focks'
       fileStr_survCurve = ''
   if dengueDataYrs_noOutbreaks_opt:
      fileStr_dengueDataYrs_noOutbreaks = f'_2007-2020_noOutbreaks{fileStr_survCurve}'
   else:
      fileStr_dengueDataYrs_noOutbreaks = f'_2001-2020{fileStr_survCurve}'
       

   fileDir  = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, 'comparisonStats_GEOS_VS_MERRAIMERG')
   fileName = f'comparisonStats_{dataSrcs_extended}_{varName}_{loc}{fileStr_dengueDataYrs_noOutbreaks}'
   filePath = os.path.join(fileDir, f'{fileName}.csv')

   if os.path.exists(filePath):
      df_comparisonStats = pd.read_csv(filePath, index_col=0)
   else:
      print(f'Error: file not found at {filePath}')
      df_comparisonStats = None  # or raise an error

   return df_comparisonStats

