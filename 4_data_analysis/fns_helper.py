#*******************************************************
# fn_helper.py
#*******************************************************
# Revision history
# 2023-11-29, CNY: Edited read_WHATCHEM_output() so that metadata is stored as a dictionary (rather than just a list of strings)
# 2024-03-04, CNY: Edited read_WHATCHEM_output() to account for new/renamed column names (PRCP, PRCPEFF)
# 2024-05-14, CNY: Added convert_hrly2dlyMeanMinMax()
# 2024-05-16, CNY: Added fillMissingVals_linReg()
# 2024-10-19, CNY: Added a whole bunch of functions while reorganizing code
# 2025-04-14, CNY: Fixed the date parsing in read_WHATCHEM_output(), which had stopped working (perhaps due to conda env updates).
# 2025-04-17, CNY: Updated fns to support Eisen survival curves (combine_dfs_overInitDays(), get_df_vectorPop_dly_merged()).
# 2025-05-15, CNY: Updated fns to support GEOS-based WHATCH'EM output (WHATCHEM_specs_to_fpath(), combine_dfs_overInitDays(), get_df_vectorPop_dly_merged()).
#                  Streamlined get_tbounds_yyyymmddhh().
# 2025-05-22, CNY: Added get_df_vectorPop_dly_merged_ensMean() for getting df of ensemble mean for GEOS-based model runs
#                  
#
# Helper functions:
#   get_tbounds_yyyymmddhh()
#   WHATCHEM_specs_to_fpath()
#   read_WHATCHEM_output()
#   convert_hrly2dlyMeanMinMax()
#   fillMissingVals_linReg()
#   MORE!



# import required packages
# ************************************************
import pandas as pd
import numpy as np
import re
import sys
from sklearn.linear_model import LinearRegression
import calendar
import os


# import constants
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants")
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA_BASE




# Convert a year and month into start and end times in the format
# we're using for the WHATCH'EM files (yyyymmddhh).
def get_tbounds_yyyymmddhh(year, month):

    days_in_month = calendar.monthrange(year, month)[1]
    t_start = year * 1000000 + month * 10000 + 100
    t_end = year * 1000000 + month * 10000 + days_in_month * 100 + 23
    return t_start, t_end




# Read in a WHATCHEM output file (container_output*.txt) and return
# it as a dataframe. The dataframe has a metadata attribute in which
# the header data lives (as a dictionary).
def read_WHATCHEM_output(filepath_in):

    with open(filepath_in) as f:
        fLines = f.read().splitlines()
    fLines = [line.strip() for line in fLines] # remove whitespace from each line
    fLines = [i for i in fLines if i]          # remove lines that are just empty strings

    # Get metadata from header and format it into a dictionary
    metadata = fLines[0:4]
    prefixesToRemove = ["Location: ", "Container Specs: ", "Physical Specs: "] # remove these substrings so they don't mess up
                                                                               # splitting of strings based on colons
    metadata_new = []
    for line in metadata:
        for prefix in prefixesToRemove:
            line = line.replace(prefix, "")
        metadata_new.append(line)

    metadata_dict = {}
    # Iterate through each string in the list
    for line in metadata_new:
        # Split each string based on colons and double spaces
        pairs = re.split(":|  ", line)
        pairs = [str.strip() for str in pairs] # get rid of extra white space
        pairs = [str for str in pairs if str]  # remove empty strings
        for i in range(0, len(pairs), 2):
            key = pairs[i]
            value = pairs[i + 1] if i + 1 < len(pairs) else None
            metadata_dict[key] = value

    # Create dataframe with metadata stored in attrs
    colNames = ["Year", "Month", "Day", "Hour", "SWDOWN (W)", "SWSIDE (W)", \
                "LWDOWN (W)", "LWUP (W)", "LWIN (W)", "LWOUT (W)", "SHUP (W)", \
                "SHOUT (W)", "LHUP (W)", "GHW (W)", "GHC (W)", "COND (W)", \
                "BAL_W (W)", "BAL_C (W)", "PRCP (mm)", "PRCPEFF (mm)", "EVAP (mm)", "WH (mm)", \
                "TG (C)", "TA (C)", "TCON (C)", "TW (C)", "VPD (kPa)"]
    df_out = pd.read_table(filepath_in, sep=r"\s+", skiprows=7, names=colNames, na_values="-9999.00", engine="python")
    df_out["Date"] = pd.to_datetime(df_out[["Year", "Month", "Day", "Hour"]]) # combine columns into new date column
    df_out = df_out.drop(columns=["Year", "Month", "Day", "Hour"])            # remove unneeded columns
    df_out = df_out.set_index("Date")
    df_out.attrs["metadata"] = metadata_dict

    return df_out




# Get daily mean/min/max values from hourly data.
# If there are any nan values in the hourly data, the corresponding min/max values will be nan.
def convert_hrly2dlyMeanMinMax(var, minVarName, maxVarName):

    var_mean = var.resample("D").apply(lambda x: np.nan if any(np.isnan(x)) else np.nanmean(x))
    var_min  = var.resample("D").apply(lambda x: np.nan if any(np.isnan(x)) else np.nanmin(x))
    var_min  = var_min.rename(minVarName)
    var_max  = var.resample("D").apply(lambda x: np.nan if any(np.isnan(x)) else np.nanmax(x))
    var_max  = var_max.rename(maxVarName)
    return var_mean, var_min, var_max




# Estimate and fill missing values in a target variable. The estimates are based on a linear 
# regression of the existing values in the target variable against a reference variable.
def fillMissingVals_linReg(var_target, var_ref):

    var_target_cleaned = var_target.dropna() # remove rows that have nan values
    var_target_nanFill = var_target.copy(deep=True)
    numMissingVals = len(var_target) - len(var_target_cleaned) # number of nan values in the target variable

    # Do the linear regression after checking that (a) there are at least some missing values and
    # (b) not /all/ the values are missing.
    if numMissingVals > 0 and numMissingVals < len(var_target):
        dates_missingData = var_ref.index.difference(var_target_cleaned.index) # get dates that have nan values
        var_ref_cleaned   = var_ref[~var_ref.index.isin(dates_missingData)] # remove rows of reference var corresponding to nan values in target var
        regression_model = LinearRegression()
        regression_model.fit(var_ref_cleaned.values.reshape(-1, 1), var_target_cleaned.values.reshape(-1, 1))
        var_target_pred = regression_model.predict(var_ref.loc[dates_missingData].values.reshape(-1, 1)) # predictions for days with nan values
        var_target_nanFill.loc[dates_missingData] = var_target_pred.flatten() # substitute in the predicted values for the nan values

    return var_target_nanFill




# strConvert_paren2hyphen()
# Converts a string of format '*(Z)' to '*-Z', removing extra white spaces.
# This is for converting population column names (e.g., 'pop_A(E)')
# to a format without parentheses (e.g., 'pop_A-E') for using as a filename.
def strConvert_paren2hyphen(input_str):
  
  # Find the positions of '(' and ')'
  open_bracket_index = input_str.find('(')
  close_bracket_index = input_str.find(')')
  
  # Check if '(' and ')' is found and they're in the expected order of '(' before ')'
  if open_bracket_index != -1 and close_bracket_index != -1 and open_bracket_index < close_bracket_index:
      # Extract the parts before '(', after ')', and between them
      substr_beforeBrackets  = input_str[:open_bracket_index]
      substr_betweenBrackets = input_str[open_bracket_index+1:close_bracket_index]
      substr_afterBrackets   = input_str[close_bracket_index+1:] # including this just in case, but is empty for pop col names
      
      # Combine the parts with '-', stripping leading/trailing white space is each substring
      output_str = substr_beforeBrackets.strip() + '-' + substr_betweenBrackets.strip() + substr_afterBrackets.strip()
      return output_str
  else:
      return input_str




# strConvert_parenUnitsRemove()
# Converts a string of format '* (Z)' to '*', removing extra white spaces.
# This is for converting column names with units (e.g., 'PRCP (mm)')
# to a format without units (e.g., 'PRCP') for using as a filename.
def strConvert_parenUnitsRemove(input_str):  
  
  if re.search(r'\ \(.*\)', input_str): # check if ' (*)' is found, where '*' is any number of characters
      open_bracket_index = input_str.find('(') # find the position of '('
      substr_beforeBrackets  = input_str[:open_bracket_index] # extract the part before '('
      output_str = substr_beforeBrackets.strip() # strip leading/trailing white space
      return output_str
  else:
      return input_str




# Read in a WHATCH'EM specifications and return the filepath to the corresponding WHATCHEM
# output file (container_output_*.txt). Also, return a subdirectory name and a string based
# on just the WHATCH'EM specifications (used for creating names of new files that are based 
# on the WHATCH'EM output file).
#
# Arguments:
#   dataSrc    - "MERRA_IMERG" or "GEOS"; source(s) of climate data
#   loc        - "Negombo", "Jaffna", or "NuwaraEliya"; location
#   year       - year (int)
#   month      - month (int)
#   intv       - "intv00", "intvF0", "intv0E", or "intvFE" (0 is no intv, F is fill, E is empty)
#   lead_time  - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#   ens_num    - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
def WHATCHEM_specs_to_fpath(dataSrc, loc, year, month, intv, lead_time=None, ens_num=None):

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

    tStrt, tEnd = get_tbounds_yyyymmddhh(year, month)
    tspan       = str(tStrt) + "_" + str(tEnd)  # "yyyymmddhh_yyyymmddhh"
    container = "Bucket"
    refl      = "Gray"
    shade     = "HalfShade"
    contSpecs = f'{container}-{refl}-{shade}'
    
    clouds    = "Observed" # future note: if clouds not included in data, need state specified cloud state (e.g., clouds='PartlyCloudy')

    subDir   = f'{dataSrc_extended}-{loc}-{tspan}'
    dir_data = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_BASE, dataSrc, subDir)

    fileStr_allSpecs = f'{dataSrc_extended}-{loc}-{tspan}-{intv}-{contSpecs}'
    if clouds != "Observed": # if clouds are specified (rather than from data), filename includes the 'clouds' string
        fileStr_allSpecs = fileStr_allSpecs + "-" + clouds
    filepath_data = os.path.join(dir_data, f'container_output_{fileStr_allSpecs}.txt')

    return filepath_data, subDir, fileStr_allSpecs


# Combine dfs that differ in immature introduction day into a single df.
# The combined df has population values that are the sums of the individual dfs.
#
# Arguments:
#   dataSrc    - "MERRA_IMERG" or "GEOS"; source(s) of climate data
#   loc        - "Negombo", "Jaffna", or "NuwaraEliya"; location
#   year       - year (int)
#   month      - month (int)
#   intv       - "intv00", "intvF0", "intv0E", or "intvFE" (0 is no intv, F is fill, E is empty)
#   initDays   - list of immature introduction days (int)
#   initEggs   - number of initial eggs (int)
#   initLarvae - number of initial larvae (int)
#   initPupae  - number of initial pupae (int)
#   survCurve  - "Focks" or "Eisen"; which temp-based survival curves to use
#   lead_time  - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#   ens_num    - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
def combine_dfs_overInitDays(dataSrc, loc, year, month, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time=None, ens_num=None):
    
    # Validate GEOS-specific args
    if dataSrc == 'GEOS':
        # Error if GEOS-specific args not provided when they're expected
        if lead_time is None or ens_num is None:
            raise ValueError('For dataSrc="GEOS", both lead_time and ens must be provided.')
    else:
        # Warn if GEOS-specific args were provided but won't be used
        if lead_time is not None or ens_num is not None:
            print('Warning: provided lead_time and ens are ignored when dataSrc is not "GEOS".')

    _, subDir, fileStr_allSpecs = WHATCHEM_specs_to_fpath(dataSrc, loc, year, month, intv, lead_time, ens_num)
    if survCurve == "Eisen":
        survStr = "_eisenQuadFit"
    else: # if survCurve is "Focks"
        survStr = ""
    df_list = []
    pop_cols = ['pop_E', 'pop_L', 'pop_P', 'pop_A', 'pop_L(E)', 'pop_L(L)', 
                'pop_P(E)', 'pop_P(L)', 'pop_P(P)', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)']
    for initDay in initDays:
        popStr = f'day{initDay}_{initEggs}e{initLarvae}l{initPupae}p'
        fname  = f'df_vectorPop_dly_{fileStr_allSpecs}-{popStr}{survStr}.txt'
        fpath  = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_BASE, dataSrc, subDir, fname)
        df_vectorPop_dly = pd.read_csv(fpath, sep='\t', na_values="NaN", index_col=0, parse_dates=True, comment='#')
        df_list.append(df_vectorPop_dly[pop_cols])                                                      # get just the population columns
    df_vectorPop_dly_pop_sumOverInitDays = sum(df_list)                                                 # sum over population columns
    df_vectorPop_dly_notPop = df_vectorPop_dly.drop(pop_cols, axis=1)                                   # get the non-population columns
    df_allInitDays = pd.concat([df_vectorPop_dly_pop_sumOverInitDays, df_vectorPop_dly_notPop], axis=1) # combine the summed pop cols with the non-pop cols and append
    return df_allInitDays




# Get dataframe that combines all df_vectorPop_dly data for one location across the specified parameters.
# Arguments:
#   dataSrc    - "MERRA_IMERG"; source(s) of climate data
#   loc        - "Negombo", "Jaffna", or "NuwaraEliya"; location
#   years      - list of years (int)
#   months     - list of months (int)
#   initDays   - list of immature introduction days (int)
#   initEggs   - number of initial eggs (int)
#   initLarvae - number of initial larvae (int)
#   initPupae  - number of initial pupae (int)
#   survCurve  - "Focks" or "Eisen"; which temp-based survival curves to use
#   lead_time  - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#   ens_num    - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
def get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time=None, ens_num=None):
    
    # Validate GEOS-specific args
    if dataSrc == 'GEOS':
        # Error if GEOS-specific args not provided when they're expected
        if lead_time is None or ens_num is None:
            raise ValueError('For dataSrc="GEOS", both lead_time and ens must be provided.')
    else:
        # Warn if GEOS-specific args were provided but won't be used
        if lead_time is not None or ens_num is not None:
            print('Warning: provided lead_time and ens are ignored when dataSrc is not "GEOS".')
    
    df_list = []
    for year in years:
        for month in months:
            df_allInitDays = combine_dfs_overInitDays(dataSrc, loc, year, month, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
            df_list.append(df_allInitDays)
    df_vectorPop_dly_merged = pd.concat(df_list)
    df_vectorPop_dly_merged = df_vectorPop_dly_merged.sort_index() # just in case it isn't automatically sorted in datetime order
    return df_vectorPop_dly_merged




# Get dataframe that combines all df_vectorPop_dly data for one location across the specified parameters and averaged over the specified ensemble members.
# Arguments:
#   dataSrc    - "GEOS"; source(s) of climate data
#   loc        - "Negombo", "Jaffna", or "NuwaraEliya"; location
#   years      - list of years (int)
#   months     - list of months (int)
#   initDays   - list of immature introduction days (int)
#   initEggs   - number of initial eggs (int)
#   initLarvae - number of initial larvae (int)
#   initPupae  - number of initial pupae (int)
#   survCurve  - "Focks" or "Eisen"; which temp-based survival curves to use
#   lead_time  - indicates GEOS lead time in months
#   ens_nums   - indicates GEOS ensemble numbers to average over
def get_df_vectorPop_dly_merged_ensMean(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_nums):
    
    # Check that dataSrc is 'GEOS'
    if dataSrc != 'GEOS':
        raise ValueError(f"Invalid dataSrc '{dataSrc}'. Only 'GEOS' is supported.")

    dfs = []
    for ens_num in ens_nums:
        df_vectorPop_dly_merged = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
        dfs.append(df_vectorPop_dly_merged)

    df_vectorPop_dly_merged_ensMean = pd.concat(dfs).groupby(level=0).mean()
    return df_vectorPop_dly_merged_ensMean




# Resample df from daily to monthly, using different aggregation functions for each column.
# Monthly values are obtained as follows:
#   population - last daily value
#   prcp       - sum of daily values
#   all others - mean of daily values
def convert_df_vectorPop_dly2mon(df_dly):

    # aggregation functions for resampling from daily to monthly
    agg_fns = {'pop_E':'last',          'pop_L':'last',              'pop_P':'last',       'pop_A':'last',         'pop_L(E)':'last',      'pop_L(L)':'last',
               'pop_P(E)':'last',       'pop_P(L)':'last',           'pop_P(P)':'last',    'pop_A(E)':'last',      'pop_A(L)':'last',      'pop_A(P)':'last',
               'dev_E':'mean',          'dev_L':'mean',              'dev_P':'mean',       'surv_E':'mean',        'surv_L':'mean',        'surv_P':'mean',
               'surv_temp_E':'mean',    'surv_temp_L':'mean',        'surv_temp_P':'mean', 'surv_desicc_E':'mean', 'surv_desicc_L':'mean', 'surv_desicc_P':'mean', 
               'WH (mm)':'mean',        'WH_del (mm)':'mean',        'TW (C)':'mean',      'TW_min (C)':'mean',    'TW_max (C)':'mean', 
               'TW_nanFill (C)':'mean', 'TW_min_nanFill (C)':'mean', 'TW_max_nanFill (C)':'mean',
               'PRCP (mm)':'sum',       'TA (C)':'mean',             'TA_min (C)':'mean',  'TA_max (C)':'mean',    'VPD (kPa)':'mean'
              }
    df_mon = df_dly.resample('MS').agg(agg_fns) # # resample to monthly, with time index set to beginning of month
    return df_mon




# Creates a df that contained the passed series as well as that same series lagged by
# the specified lag times.
# This is intended for creating a df of lagged monthly dengue case data.
#
# Arguments:
#   series      - series of dengue data
#   lagTimes    - lag time; how much time to shift the data back by (months)
#   baseColName - the prefix to use for all lagged column names
#   startDate   - beginning date of interest for series
#   endDate     - ending date of interest for series
def create_df_withLags(series, lagTimes, baseColName, startDate, endDate):

    df_lagged = pd.DataFrame()
    df_lagged[baseColName + '-lag0'] = series[startDate:endDate] # start with original, non-lagged series

    for lagTime in lagTimes:
        if lagTime != 0: # we don't do anything if there's no lag, since that's already included in the df
            colName = baseColName + '-lag' + str(lagTime)
            series_lagged = series.copy(deep=True)
            series_lagged.index = series.index - pd.DateOffset(months=lagTime)
            series_lagged = series_lagged[startDate:endDate]
            df_lagged[colName] = series_lagged
    return df_lagged




# Categorizes the input series into two bins (Low, High), 
# returning a corresponding series of these labels.
# If the categorization fails, the returned series is all NaN values.
# Arguments:
#   series - series of data; must have datetime index
def categorizeSeries_LowHigh(series):
    category_labels = ['Low', 'High']
    try:
        series_categories = pd.qcut(series, q=len(category_labels), labels=category_labels)
        series_categories = pd.Categorical(series_categories, categories=category_labels, ordered=True)
    except ValueError:
        series_categories = pd.Series(index=series.index, dtype='category')
    
    series_categories = pd.Series(series_categories, index=series.index) # ensure same indices as series
    return series_categories




# Categorizes the input series into three bins (Low, Medium, High), 
# returning a corresponding series of these labels.
# If the categorization fails, the returned series is all NaN values.
# Arguments:
#   series - series of data; must have datetime index
def categorizeSeries_LowMedHigh(series):
    category_labels = ['Low', 'Medium', 'High']
    try:
        series_categories = pd.qcut(series, q=len(category_labels), labels=category_labels)
        series_categories = pd.Categorical(series_categories, categories=category_labels, ordered=True)
    except ValueError:
        series_categories = pd.Series(index=series.index, dtype='category')
    
    series_categories = pd.Series(series_categories, index=series.index) # ensure same indices as series
    return series_categories




# Categorizes the input series into three bins (Zero, Medium, High), 
# returning a corresponding series of these labels.
# If the categorization fails, the returned series is all NaN values.
# This is intended for the adult population variables, which sometimes
# have so many zero values that categorizeSeries_lowMedHigh() doesn't work.
# Arguments:
#   series - series of data; must have datetime index
def categorizeSeries_ZeroMedHigh(series):
    category_labels = ['Zero', 'Medium', 'High']

    series_categories = pd.Series(index=series.index, dtype='object')

    freqVal_mask = series == 0
    series_categories[freqVal_mask] = 'Zero'
    
    # Now we're all done with categorization if all the values are zero.
    # Otherwise we need to categorize the nonzero values, which we do now.
    series_otherVals = series[~freqVal_mask]
    if len(series_otherVals) == 1: # if there's only one nonzero value, set its category to High
        series_categories[~freqVal_mask] = 'High'
    elif len(series_otherVals) > 1: # if there's more than one nonzero value, split into Medium and High based on median (which is always Medium)
        median = series_otherVals.median()
        try:
            series_categories_otherVals = pd.cut(series_otherVals, bins=[0, median, np.inf], right=True, labels=['Medium', 'High'])
            series_categories[~freqVal_mask] = series_categories_otherVals
        except ValueError:
            return pd.Series(index=series.index, dtype='category')
    
    series_categories = pd.Categorical(series_categories, categories=category_labels, ordered=True)
    series_categories = pd.Series(series_categories, index=series.index) # ensure same indices as series

    return series_categories




# Categorizes the input series into three bins (Low, Medium, 0.99), 
# returning a corresponding series of these labels.
# If the categorization fails, the returned series is all NaN values.
# This is intended for the survival factor variables, which sometimes
# have so many 0.99 values that categorizeSeries_lowMedHigh() doesn't work.
# Arguments:
#   series - series of data; must have datetime index
def categorizeSeries_LowMed99(series):
    category_labels = ['Low', 'Medium', '0.99']

    series_categories = pd.Series(index=series.index, dtype='object')

    freqVal_mask = series == 0.99
    series_categories[freqVal_mask] = '0.99'
    
    # Now we're all done with categorization if all the values are 0.99.
    # Otherwise we need to categorize the non-0.99 values, which we do now.
    series_otherVals = series[~freqVal_mask]
    if len(series_otherVals) == 1: # if there's only one non-0.99 value, set its category to Low
        series_categories[~freqVal_mask] = 'Low'
    elif len(series_otherVals) > 1: # if there's more than one non-0.99 value, split into Low and Medium based on median (which is always Medium)
        median = series_otherVals.median()
        try:
            series_categories_otherVals = pd.cut(series_otherVals, bins=[0, median, 0.99], right=True, labels=['Low', 'Medium'])
            series_categories[~freqVal_mask] = series_categories_otherVals
        except ValueError:
            return pd.Series(index=series.index, dtype='category')
    
    series_categories = pd.Categorical(series_categories, categories=category_labels, ordered=True)
    series_categories = pd.Series(series_categories, index=series.index) # ensure same indices as series

    return series_categories




# Categorizes the input series into two bins (Low, High) on a 
# per calendar month basis, returning a corresponding series of these labels.
#
# For each calendar month we try to categorize using categorizeSeries_LowHigh().
# Arguments:
#   series - series of data; must have datetime index
def categorizeSeries_LowHigh_byMon(series):
    months = sorted(series.index.month.unique()) # get months based on series datetime index
    category_labels = ['Low', 'High']

    series_categories_list = []
    for month in months:
        series_oneMon = series[series.index.month == month]
        series_oneMon_categories = categorizeSeries_LowHigh(series_oneMon)
        if series_oneMon_categories.isnull().all():
            return pd.Series(index=series.index, dtype='category')
        series_categories_list.append(series_oneMon_categories)

    series_categories = pd.concat(series_categories_list).sort_index()
    series_categories = pd.Categorical(series_categories, categories=category_labels, ordered=True)
    series_categories = pd.Series(series_categories, index=series.index) # ensure same indices as series

    return series_categories




# Categorizes the input series into three bins (Low, Medium, High) on a 
# per calendar month basis, returning a corresponding series of these labels.
#
# For each calendar month we first try to categorize using categorizeSeries_LowMedHigh().
# If this fails we assume that's because either (a) it's a pop variable and there are too many zeros,
# or (b) it's a surv variable and there are too many 0.99 values.
# In case (a) we instead use categorizeSeries_ZeroMedHigh() and rename the 'Zero' category to 'Low'.
# In case (b) we instead use categorizeSeries_LowMed99() and rename the '0.99' category to 'High'.
# Arguments:
#   series - series of data; must have datetime index
def categorizeSeries_LowMedHigh_byMon(series):
    months = sorted(series.index.month.unique()) # get months based on series datetime index
    category_labels = ['Low', 'Medium', 'High']

    series_categories_list = []
    for month in months:
        series_oneMon = series[series.index.month == month]
        series_oneMon_categories = categorizeSeries_LowMedHigh(series_oneMon)
        if series_oneMon_categories.isnull().all():
            if series.name in ['pop_A', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)']:
                series_oneMon_categories = categorizeSeries_ZeroMedHigh(series_oneMon)
                series_oneMon_categories = series_oneMon_categories.cat.rename_categories({'Zero': 'Low'})
            elif series.name in  ['surv_E', 'surv_L', 'surv_P']:
                series_oneMon_categories = categorizeSeries_LowMed99(series_oneMon)
                series_oneMon_categories = series_oneMon_categories.cat.rename_categories({'0.99': 'High'})
            if series_oneMon_categories.isnull().all():
                return pd.Series(index=series.index, dtype='category')
        series_categories_list.append(series_oneMon_categories)

    series_categories = pd.concat(series_categories_list).sort_index()
    series_categories = pd.Categorical(series_categories, categories=category_labels, ordered=True)
    series_categories = pd.Series(series_categories, index=series.index) # ensure same indices as series

    return series_categories




# Loads crosstab data from file as a dataframe.
def load_crosstab(fPath):
    crosstab = pd.read_csv(fPath, index_col=0, header=0)
    return crosstab


