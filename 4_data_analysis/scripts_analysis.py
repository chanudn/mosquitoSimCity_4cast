# scripts_analysis.py
#
# Assorted blocks of code for analyzing data from WHATCH'EM and the vector survival simulations.
#
# Revision history
#    2025-04-17, CNY: updated scripts to account for Focks vs Eisen survival curves
#    2025-05-26, CNY: updated scripts to account for GEOS
#    2025-05-31, CNY: added stats calculation to compare GEOS vs MERRA/IMERG
#



# import required packages
# ************************************************
import sys
import numpy as np
import pandas as pd
import os
from datetime import datetime



# import constants and functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants") # directory of constants file
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/4_analysis")  # directory of helper fns file
from fns_analysis import calc_tercile_crosstab_splitByMon
from fns_helper import (strConvert_paren2hyphen, strConvert_parenUnitsRemove, get_df_vectorPop_dly_merged, get_df_vectorPop_dly_merged_ensMean,
                        convert_df_vectorPop_dly2mon, create_df_withLags, load_crosstab)
from fns_analysis import calc_table_comparisonStats





# calculate and save tercile crosstabs of monthly values of variables against dengue incidence
# (each resulting df has the crosstab values for ONE variable for a given lag+outbreakYr scenario)
# Note: NEED TO RUN THIS BEFORE RUNNING NEXT BLOCK, SINCE THAT DEPENDS ON THE CROSSTABS PRODUCED HERE!
if False:

    # parameters related to modeling pipeline
    dataSrcs   = ['GEOS']#['MERRA_IMERG', 'GEOS']
    lead_time  = 0             # only relevant for GEOS
    ens_nums   = range(1, 5+1) # only relevant for GEOS
    locs       = ['Negombo', 'NuwaraEliya', 'Jaffna']
    years      = list(range(2001, 2020+1))
    months     = list(range(1,12+1))
    intv       = 'intv00'
    initDays   = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']
    var_colNames = ['TA_min (C)', 'TA (C)' ,'TA_max (C)', 'PRCP (mm)', 'TW_min_nanFill (C)', 'TW_nanFill (C)', 'TW_max_nanFill (C)', 'WH (mm)',
                   'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'       # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False]  # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)

    # load dengue data
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    file_dengue  = os.path.join(fileDir_base, 'dengueInc_monthly_2007_2022.csv')
    df_den_monthly = pd.read_csv(file_dengue, index_col="Date", parse_dates=True)

    for dataSrc in dataSrcs:
        print(dataSrc)
        for loc in locs:
            print(loc)
            if loc == 'Negombo':
                den_loc = 'Gampaha'
            elif loc in ['NuwaraEliya', 'Jaffna']:
                den_loc = loc
            else:
                print('Error: unknown loc.')
                exit

            # subset dengue data to location of interest and create df with lagged data (which also subsets to time of interest)
            series_den_monthly_loc = df_den_monthly[den_loc]
            df_den_monthly_loc_wLags = create_df_withLags(series_den_monthly_loc, lagTimes, den_colName_base, startDate, endDate)

            # load model pipeline data at location of interest
            # also: convert from daily to monthly, subset to time of interest
            if dataSrc == 'MERRA_IMERG':
                df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
                dataSrc_extended = f'{dataSrc}'
            elif dataSrc == 'GEOS':
                df_vectorPop_dly = get_df_vectorPop_dly_merged_ensMean(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_nums)
                dataSrc_extended = f'{dataSrc}_lead{lead_time}_ensMean'
            df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
            df_vectorPop_monthly = df_vectorPop_monthly[startDate:endDate]

            # combine model pipeline data and dengue data into one df
            df_vectorPop_den_monthly = df_vectorPop_monthly.join(df_den_monthly_loc_wLags)

            if survCurve == 'Eisen':
                fileStr_survCurve = '_eisenQuadFit'
            else: # if survCurve is 'Focks'
                fileStr_survCurve = ''

            for noOutbreaks in noOutbreaks_opts:
                print(f' noOutbreaks: {noOutbreaks}')
                if noOutbreaks:
                    fileStr_noOutbreaks = '_noOutbreaks'
                    exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
                else:
                    fileStr_noOutbreaks = ''
                    exclYrs = []

                # remove data in outbreak years we want to exclude
                df_vectorPop_den_monthly_yrFiltered = df_vectorPop_den_monthly[~df_vectorPop_den_monthly.index.year.isin(exclYrs)]
                
                for i_lag, lagTime in enumerate(lagTimes):
                    print(f'   Lag (mo): {lagTime}')
                    for i_var, var_colName in enumerate(var_colNames):
                        var_pltName = var_pltNames[i_var]
                        print(f'  Var: {var_pltName}')

                        den_colName = f'{den_colName_base}-lag{lagTime}' # this must match col name format set in create_df_withLags()
                        den_pltName = den_colName # just reusing the column name here for simplicity
                        
                        fileName = f'crosstab_tercile_splitbyMon_{dataSrc_extended}_{var_pltName}_{den_pltName}_{loc}{fileStr_noOutbreaks}{fileStr_survCurve}'
                        fileDir = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc, 'crosstabs')
                        calc_tercile_crosstab_splitByMon(df_vectorPop_den_monthly_yrFiltered, var_colName, den_colName, fileDir, fileName)



# calculate and save tercile crosstab element values (e.g. High-High) of monthly values of variables against dengue incidence
# (each resulting df has the crosstab element values for ALL the variables for a given lag+outbreakYr scenario)
# Note: NEED TO RUN PRIOR BLOCK BEFORE RUNNING THIS ONE, SINCE THIS DEPENDS ON THE CROSSTABS GENERATED BY THE PRIOR BLOCK!
if False:

    # parameters related to modeling pipeline
    dataSrcs   = ['GEOS']#['MERRA_IMERG', 'GEOS']
    lead_time  = 0             # only relevant for GEOS
    ens_nums   = range(1, 5+1) # only relevant for GEOS
    locs       = ['Negombo', 'NuwaraEliya', 'Jaffna']
    years      = list(range(2001, 2020+1))
    months     = list(range(1,12+1))
    intv       = 'intv00'
    initDays   = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']
    var_colNames = ['TA_min (C)', 'TA (C)', 'TA_max (C)', 'PRCP (mm)', 'TW_min_nanFill (C)', 'TW_nanFill (C)', 'TW_max_nanFill (C)', 'WH (mm)',
                    'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'       # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False]  # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)

    # parameters related to the crosstab element of interest
    # format is 'terc1-terc2' where terc1 is tercile of the variable, terc2 is tercile of dengue
    crosstab_elts = ['HI-HI', 'MED-HI', 'LO-HI', 'HI-LO', 'MED-LO', 'LO-LO']
    map_tercAbbr_to_tercName = {'HI': 'High', 'MED': 'Medium', 'LO': 'Low'} # map to names of rows/columns in crosstab csv file

    # load dengue data
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    file_dengue  = os.path.join(fileDir_base, 'dengueInc_monthly_2007_2022.csv')
    df_den_monthly = pd.read_csv(file_dengue, index_col="Date", parse_dates=True)

    for dataSrc in dataSrcs:
        print(dataSrc)
        for loc in locs:
            print(loc)
            if loc == 'Negombo':
                den_loc = 'Gampaha'
            elif loc in ['NuwaraEliya', 'Jaffna']:
                den_loc = loc
            else:
                print('Error: unknown loc.')
                exit

            # subset dengue data to location of interest and create df with lagged data (which also subsets to time of interest)
            series_den_monthly_loc = df_den_monthly[den_loc]
            df_den_monthly_loc_wLags = create_df_withLags(series_den_monthly_loc, lagTimes, den_colName_base, startDate, endDate)

            # load model pipeline data at location of interest
            # also: convert from daily to monthly, subset to time of interest
            if dataSrc == 'MERRA_IMERG':
                df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
                dataSrc_extended = f'{dataSrc}'
            elif dataSrc == 'GEOS':
                df_vectorPop_dly = get_df_vectorPop_dly_merged_ensMean(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_nums)
                dataSrc_extended = f'{dataSrc}_lead{lead_time}_ensMean'
            df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)
            df_vectorPop_monthly = df_vectorPop_monthly[startDate:endDate]

            # combine model pipeline data and dengue data into one df
            df_vectorPop_den_monthly = df_vectorPop_monthly.join(df_den_monthly_loc_wLags)

            if survCurve == 'Eisen':
                fileStr_survCurve = '_eisenQuadFit'
            else: # if survCurve is 'Focks'
                fileStr_survCurve = ''

            for noOutbreaks in noOutbreaks_opts:
                print(' noOutbreaks: ' + str(noOutbreaks))
                if noOutbreaks:
                    fileStr_noOutbreaks = '_noOutbreaks'
                    exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
                else:
                    fileStr_noOutbreaks = ''
                    exclYrs = []

                # remove data in outbreak years we want to exclude
                df_vectorPop_den_monthly_yrFiltered = df_vectorPop_den_monthly[~df_vectorPop_den_monthly.index.year.isin(exclYrs)]
                
                for i_lag, lagTime in enumerate(lagTimes):
                    print(f'   Lag (mo): {lagTime}')
                    for crosstab_elt in crosstab_elts:
                        print(f'    Crosstab element: {crosstab_elt}')

                        tercileStr_modelVar, tercileStr_den = crosstab_elt.split('-')
                        dir_crosstabs = f'crosstabs_{tercileStr_modelVar}_{tercileStr_den}'
                        tercileName_modelVar = map_tercAbbr_to_tercName.get(tercileStr_modelVar) # get corresponding row name in crosstabs csv file
                        tercileName_den      = map_tercAbbr_to_tercName.get(tercileStr_den)      # get corresponding col name in crosstabs csv file
                        colName_4_df = f'{tercileName_modelVar}-{tercileName_den} Value' # (e.g., 'High-High Value')

                        crosstab_elt_vals = {}      # dict to store tercile-tercile values for each variable
                        crosstab_elt_vals_norm = {} # dict to store tercile-tercile values for each variable (normalized to a %)

                        for i_var, var_colName in enumerate(var_colNames):
                            var_pltName = var_pltNames[i_var]
                            print(f'  Var: {var_pltName}')

                            den_colName = f'{den_colName_base}-lag{lagTime}' # this must match col name format set in create_df_withLags()
                            den_pltName = den_colName # just reusing the column name here for simplicity
                            
                            fileName = f'crosstab_tercile_splitbyMon_{dataSrc_extended}_{var_pltName}_{den_pltName}_{loc}{fileStr_noOutbreaks}{fileStr_survCurve}'
                            filePath = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc, 'crosstabs', f'{fileName}.csv')
                            filePath_norm = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc, 'crosstabs', f'{fileName}_norm.csv')
                            crosstab = load_crosstab(filePath)
                            crosstab_norm = load_crosstab(filePath_norm)

                            crosstab_elt_vals[var_colName]      = crosstab.loc[tercileName_modelVar, tercileName_den] if tercileName_modelVar in crosstab.index and tercileName_den in crosstab.columns else np.nan
                            crosstab_elt_vals_norm[var_colName] = crosstab_norm.loc[tercileName_modelVar, tercileName_den] if tercileName_modelVar in crosstab_norm.index and tercileName_den in crosstab_norm.columns else np.nan
                    
                        # Convert the dictionary to a DataFrame
                        df_crosstab_elt = pd.DataFrame(list(crosstab_elt_vals.items()), columns=['Variable', colName_4_df])
                        df_crosstab_elt_norm = pd.DataFrame(list(crosstab_elt_vals_norm.items()), columns=['Variable', colName_4_df])
            
                        # Save to csv files
                        fileName = f'crosstab_tercile_splitbyMon_{dataSrc_extended}_{tercileStr_modelVar}-{tercileStr_den}_{den_pltName}_{loc}{fileStr_noOutbreaks}{fileStr_survCurve}'
                        fileDir  = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc, dir_crosstabs)
                        os.makedirs(fileDir, exist_ok=True) # make save directory if it doesn't exist
                        filePath      = os.path.join(fileDir, f'{fileName}.csv')
                        filePath_norm = os.path.join(fileDir, f'{fileName}_norm.csv')
                        df_crosstab_elt.to_csv(filePath, index=False)
                        df_crosstab_elt_norm.to_csv(filePath_norm, index=False)




# calculate and save statistics comparing GEOS to MERRA/IMERG
if False:
    # parameters related to modeling pipeline
    dataSrc_main  = 'GEOS'
    dataSrc_ref   = 'MERRA_IMERG'
    lead_times    = [0]           # only relevant for GEOS
    ens_nums      = range(1, 5+1) # only relevant for GEOS
    locs          = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years         = list(range(2001, 2020+1))
    months        = list(range(1, 12+1))
    intv          = 'intv00'
    initDays      = list(range(1, 7+1))
    initEggs      = 1000
    initLarvae    = 1000
    initPupae     = 1000
    survCurve     = 'Focks' # ['Focks', 'Eisen']
    var_colNames  = ['PRCP (mm)', 'WH (mm)', 'TA_max (C)', 'TA_min (C)', 'TA (C)', 'TW_max (C)', 'TW_min (C)', 'TW (C)',
                    'dev_E', 'dev_L', 'dev_P', 'surv_E', 'surv_L', 'surv_P', 'pop_A(E)_frac', 'pop_A(L)_frac', 'pop_A(P)_frac', 'VPD (kPa)']
    var_pltNames  = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]
    doAnom_opts   = [True, False] # whether to compute statistics for data as anomalies (relative to monthly climatology) or as is

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True, False] # whether to consider (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]


    for doAnom in doAnom_opts:
        print(f'doAnom: {doAnom}')
        if doAnom:
            fileStr_doAnom = 'Anom_'
        else:
            fileStr_doAnom = ''
        for loc in locs:
            print(f' {loc}')
            for lead_time in lead_times:
                print(f'  Lead time: {lead_time}')

                # load model pipeline data at location of interest
                # also: convert from daily to monthly
                dataSrcs_extended = f'{fileStr_doAnom}{dataSrc_main}_lead{lead_time}_VS_{dataSrc_ref}'

                dfs_vectorPop_monthly_GEOS_ens = []
                for ens_num in ens_nums:
                    df_vectorPop_dly_GEOS_ens     = get_df_vectorPop_dly_merged(dataSrc_main, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
                    df_vectorPop_monthly_GEOS_ens = convert_df_vectorPop_dly2mon(df_vectorPop_dly_GEOS_ens)
                    dfs_vectorPop_monthly_GEOS_ens.append(df_vectorPop_monthly_GEOS_ens)
                df_vectorPop_monthly_GEOS_ensMean = pd.concat(dfs_vectorPop_monthly_GEOS_ens).groupby(level=0).mean()

                df_vectorPop_dly_MERRAIMERG     = get_df_vectorPop_dly_merged(dataSrc_ref, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
                df_vectorPop_monthly_MERRAIMERG = convert_df_vectorPop_dly2mon(df_vectorPop_dly_MERRAIMERG)

                if survCurve == 'Eisen':
                    fileStr_survCurve = '_eisenQuadFit'
                else: # if survCurve is 'Focks'
                    fileStr_survCurve = ''

                for dengueDataYrs_noOutbreaks in dengueDataYrs_noOutbreaks_opts:
                    print(' dengueDataYrs_noOutbreaks: ' + str(dengueDataYrs_noOutbreaks))
                    if dengueDataYrs_noOutbreaks:
                        fileStr_dengueDataYrs_noOutbreaks = '_' + str(dengueDataYrs[0]) + '-' + str(dengueDataYrs[-1]) + '_noOutbreaks' + fileStr_survCurve
                        exclYrs = list(set([year for year in years if year not in dengueDataYrs] + outbreakYrs))  # years to exclude from plot
                    else:
                        fileStr_dengueDataYrs_noOutbreaks = '_' + str(years[0]) + '-' + str(years[-1]) + fileStr_survCurve
                        exclYrs = []
                    
                    # remove data in outbreak years we want to exclude
                    df_vectorPop_monthly_GEOS_ensMean_yrFiltered = df_vectorPop_monthly_GEOS_ensMean[~df_vectorPop_monthly_GEOS_ensMean.index.year.isin(exclYrs)]
                    df_vectorPop_monthly_MERRAIMERG_yrFiltered = df_vectorPop_monthly_MERRAIMERG[~df_vectorPop_monthly_MERRAIMERG.index.year.isin(exclYrs)]
                    dfs_vectorPop_monthly_GEOS_ens_yrFiltered = [df[~df.index.year.isin(exclYrs)] for df in dfs_vectorPop_monthly_GEOS_ens]
                    
                    # manually create the popFrac columns and add them to vectorPop df
                    popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
                    df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac = df_vectorPop_monthly_GEOS_ensMean_yrFiltered.copy(deep=True)
                    df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_monthly_GEOS_ensMean_yrFiltered['pop_A(E)'] / popTotal_E
                    df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_monthly_GEOS_ensMean_yrFiltered['pop_A(L)'] / popTotal_E
                    df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_monthly_GEOS_ensMean_yrFiltered['pop_A(P)'] / popTotal_E
                    df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac = df_vectorPop_monthly_MERRAIMERG_yrFiltered.copy(deep=True)
                    df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered['pop_A(E)'] / popTotal_E
                    df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered['pop_A(L)'] / popTotal_E
                    df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered['pop_A(P)'] / popTotal_E
                    dfs_vectorPop_monthly_GEOS_ens_yrFiltered_withFrac = []
                    for df in dfs_vectorPop_monthly_GEOS_ens_yrFiltered:
                        df_withFrac = df.copy(deep=True)
                        df_withFrac['pop_A(E)_frac'] = df_withFrac['pop_A(E)'] / popTotal_E
                        df_withFrac['pop_A(L)_frac'] = df_withFrac['pop_A(L)'] / popTotal_E
                        df_withFrac['pop_A(P)_frac'] = df_withFrac['pop_A(P)'] / popTotal_E
                        dfs_vectorPop_monthly_GEOS_ens_yrFiltered_withFrac.append(df_withFrac)

                    # If we want to compute stats for anomalies data, then convert our dfs to datasets of anomalies (relative to monthly climatology)
                    if doAnom:
                        ## Calculate monthly climatologies for GEOS ensemble mean and MERRA/IMERG (used for computing anomalies)
                        df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac['Month'] = df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac.index.month  # needed for grouping by month
                        df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac['Month'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac.index.month  # needed for grouping by month
                        df_GEOS_ensMean_monClim = df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac.groupby('Month').mean()
                        df_MERRAIMERG_monClim = df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac.groupby('Month').mean()

                        def subtract_monClim(df, df_monClim):
                            df = df.copy()
                            df['Month'] = df.index.month
                            df_anom = df.apply(lambda row: row.drop('Month') - df_monClim.loc[row['Month']], axis=1)
                            return df_anom
                        ## Calculate anomalies of all data relative to monthly climatologies
                        # GEOS ensemble anomalies
                        dfs_GEOS_ens_anom = []
                        for df in dfs_vectorPop_monthly_GEOS_ens_yrFiltered_withFrac:
                            df_anom = subtract_monClim(df, df_GEOS_ensMean_monClim)
                            dfs_GEOS_ens_anom.append(df_anom)
                        dfs_vectorPop_monthly_GEOS_ens_yrFiltered_withFrac = dfs_GEOS_ens_anom # overwrite with anomalies dataset

                        # GEOS ensemble mean anomalies
                        df_GEOS_ensMean_anom = subtract_monClim(df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac, df_GEOS_ensMean_monClim)
                        df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac = df_GEOS_ensMean_anom  # overwrite with anomalies dataset

                        # MERRA/IMERG anomalies
                        df_MERRAIMERG_anom = subtract_monClim(df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac, df_MERRAIMERG_monClim)
                        df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac = df_MERRAIMERG_anom  # overwrite with anomalies dataset

                    for i_var, var_colName in enumerate(var_colNames):
                        print(var_colName)
                        var_pltName = var_pltNames[i_var]
                        fileDir  = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, 'comparisonStats_GEOS_VS_MERRAIMERG')
                        os.makedirs(fileDir, exist_ok=True)
                        fileName = f'comparisonStats_{dataSrcs_extended}_{var_pltName}_{loc}{fileStr_dengueDataYrs_noOutbreaks}.csv'
                        comparisonStatsTable = calc_table_comparisonStats(dfs_vectorPop_monthly_GEOS_ens_yrFiltered_withFrac, df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac, 
                                                                            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac, var_colName, fileDir, fileName)


