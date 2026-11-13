#************************************************************************
# scripts_plot4paper.py
#************************************************************************
# Revision history
#    2026-11-12, CNY: initial revision; spun off from fns_plot4paper.py
#


# import required packages
# ************************************************
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import os
import sys
import pickle


# import constants and custom functions
# ************************************************
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants') # directory of constants file
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/4_data_analysis') # directory of helper fns file
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/5_plotting') # directory of plotting fns file

from CONST_FILE_DIR import (DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, PATH_WMO_NEGOMBO, PATH_WMO_NUWARAELIYA, PATH_WMO_JAFFNA,
                            PATH_IMERG_PRECIP, PATH_MERRA_T, PATH_MERRA_TSOIL, DIR_GEOS, FILE_GEOS_PREFIX)
from fns_helper import (get_df_vectorPop_dly_merged, get_df_vectorPop_dly_merged_ensMean, convert_df_vectorPop_dly2mon,
                        strConvert_paren2hyphen, strConvert_parenUnitsRemove,  create_df_withLags)
from fns_plot4paper import (plot_tseries_denInc, plot_tseries_climCategories, plot_tseries_climCategories_allVars, 
                            plot_bar_climCategories_adultPop, plot_violinplt_allVars, plot_violinplt_allVars_allLocsOneMon, 
                            plot_violinplt_FocksVsEisen, plot_violinplt_compareDatasets, plot_violinplt_compareDatasets_allLocsOneMon, 
                            plot_scatter_compareMet_tseries, plot_scatter_compareMet_RPSSvsR, plot_scatter_compareMet_obsVsEns, 
                            plot_pdf_compareMet_9panel, plot_bar_varX_denInc, plot_scatter_crosstabs, plot_scatter_crosstabs_FocksVsEisen, 
                            plot_scatter_crosstabs_compareDatasets, load_obsData, extract_mainData_temp, extract_mainData_prcp, extract_city_data)



# make timeseries plots of dengue data
if False:

    # Read in dengue data.
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    yearStart_file = 2007
    yearEnd_file = 2022
    timespanStr_file = str(yearStart_file) + '_' + str(yearEnd_file)
    file_denInc_monthly  = fileDir_base + 'dengueInc_monthly_' + timespanStr_file + '.csv'
    df_denInc_monthly_allData  = pd.read_csv(file_denInc_monthly, index_col='Date', parse_dates=True)
    
    # Subset to time range of interest.
    yearStart_data = 2007
    yearEnd_data = 2020
    outbreakYears = [2017, 2019]
    timespanStr_data = str(yearStart_data) + '_' + str(yearEnd_data)
    df_denInc_monthly = df_denInc_monthly_allData[str(yearStart_data)+'-01-01' : str(yearEnd_data)+'-12-31'] # subset to years of interest

    fName = 'tseries_' + 'denInc_monthly_' + timespanStr_data
    plot_tseries_denInc(df_denInc_monthly, outbreakYears, fName)




# make timeseries plots of ONE VARIABLE of the model output by climate category (lo/hi temp x lo/hi prcp)
# for each variable individually
if False:

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Jaffna']#['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']
    varTemp_colName = 'TA (C)'
    varTemp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varTemp_colName))
    varPrcp_colName = 'PRCP (mm)'
    varPrcp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varPrcp_colName))
    var_colNames = ['TA (C)', 'TA_min (C)', 'TA_max (C)', 'PRCP (mm)', 'TW_nanFill (C)', 'TW_min_nanFill (C)', 'TW_max_nanFill (C)', 
                    'WH (mm)', 'dev_E', 'dev_L', 'dev_P', 'surv_E', 'surv_L', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    noOutbreaks_opts = [False]#[True, False] # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]


    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)

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
            df_vectorPop_dly_yrFiltered = df_vectorPop_dly[~df_vectorPop_dly.index.year.isin(exclYrs)]
            
            for i_var, var_colName in enumerate(var_colNames):
                var_pltName = var_pltNames[i_var]
                print('  Var: ' + var_pltName)

                for month in [1]:#months:
                    fileName = 'tseries_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + var_pltName + '_' + loc + fileStr_noOutbreaks + '_mon' + str(month).zfill(2) + fileStr_survCurve
                    plot_tseries_climCategories(df_vectorPop_dly_yrFiltered, month, var_colName, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName)




# make BIG (ALL VARIABLES) timeseries plot of model output by climate category (lo/hi temp x lo/hi prcp)
if False:
    
    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Negombo', 'NuwaraEliya']#['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']
    varTemp_colName = 'TA (C)'
    varTemp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varTemp_colName))
    varPrcp_colName = 'PRCP (mm)'
    varPrcp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varPrcp_colName))

    # parameters related to dengue data
    noOutbreaks_opts = [False]#[True, False] # whether to exclude data from outbreak years
    outbreakYrs = [2017, 2019]


    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)

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
            df_vectorPop_dly_yrFiltered = df_vectorPop_dly[~df_vectorPop_dly.index.year.isin(exclYrs)]

            # manually create the popFrac columns and add them to vectorPop df
            popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
            df_vectorPop_dly_yrFiltered_withFrac = df_vectorPop_dly_yrFiltered.copy(deep=True)
            df_vectorPop_dly_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_dly_yrFiltered['pop_A(E)'] / popTotal_E
            df_vectorPop_dly_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_dly_yrFiltered['pop_A(L)'] / popTotal_E
            df_vectorPop_dly_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_dly_yrFiltered['pop_A(P)'] / popTotal_E
        
            for month in [7]:#months:
                fileName = 'tseries_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + 'allVars' + '_' + loc + fileStr_noOutbreaks + '_mon' + str(month).zfill(2) + fileStr_survCurve
                plot_tseries_climCategories_allVars(df_vectorPop_dly_yrFiltered_withFrac, loc, month, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName)




# make tercile plots for climate categories vs adult pop
if False:
    
    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']
    varTemp_colName = 'TA (C)'
    varTemp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varTemp_colName))
    varPrcp_colName = 'PRCP (mm)'
    varPrcp_pltName = strConvert_paren2hyphen(strConvert_parenUnitsRemove(varPrcp_colName))
    adultPop_colNames = ['pop_A(E)', 'pop_A(L)', 'pop_A(P)']
    adultPop_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in adultPop_colNames]

    # parameters related to dengue data
    noOutbreaks_opts = [False]#[True, False] # whether to exclude data from outbreak years
    outbreakYrs = [2017, 2019]

    for loc in locs:
        print(loc)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)

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
            df_vectorPop_dly_yrFiltered = df_vectorPop_dly[~df_vectorPop_dly.index.year.isin(exclYrs)]

            # make tercile bar plots with all adult pop variables (E, L, P)
            fileName = 'barTerciles_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + 'popA' + '_' + loc + fileStr_noOutbreaks + fileStr_survCurve
            #plot_bar_climCategories_adultPops(df_vectorPop_dly_yrFiltered, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, fileName)
            
            # make separate tercile bar plots for each adult pop variable (E, L, P)
            for i, adultPop_colName in enumerate(adultPop_colNames):
                adultPop_pltName = adultPop_pltNames[i]
                fileName = 'barTerciles_climCat_' + varTemp_pltName + '-' + varPrcp_pltName + '_' + adultPop_pltName + '_' + loc + fileStr_noOutbreaks + fileStr_survCurve
                plot_bar_climCategories_adultPop(df_vectorPop_dly_yrFiltered, varTemp_colName, varTemp_pltName, varPrcp_colName, varPrcp_pltName, adultPop_colName, fileName)




# make violin plots for all months of ONE location
if False:
    # parameters related to modeling pipeline
    dataSrcs   = ['MERRA_IMERG'] # ['MERRA_IMERG', 'GEOS']
    lead_times = [0]           # only relevant for GEOS
    ens_nums   = range(1, 5+1) # only relevant for GEOS
    locs       = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years      = list(range(2001, 2020+1))
    months     = list(range(1, 12+1))
    intv       = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True]#[True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]


    # load dengue incidence data
    file_dengue  = '/mnt/redwood/local_drive/chanud/RawData/epiReports/' + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly_allLocs = pd.read_csv(file_dengue, index_col='Date', parse_dates=True)
    df_den_monthly_allLocs = df_den_monthly_allLocs[df_den_monthly_allLocs.index.year.isin(years)] # cut extraneous years

    for dataSrc in dataSrcs:
        print(dataSrc)
        for loc in locs:
            print(loc)
            
            # subset dengue data to location of interest
            if loc == 'Negombo':
                den_loc = 'Gampaha'
            elif loc in ['NuwaraEliya', 'Jaffna']:
                den_loc = loc
            else:
                raise SystemExit("Error: unknown loc.")
            series_den_monthly = df_den_monthly_allLocs[den_loc]
            series_den_monthly.name = 'dengueInc'
            df_den_monthly = series_den_monthly.to_frame()

            # load model pipeline data at location of interest
            # also: convert from daily to monthly

                
            if dataSrc == 'GEOS':
                for lead_time in lead_times:
                    print(f'Lead time: {lead_time}')
                    for ens_num in ens_nums:
                        print(f'ens{ens_num}')
                        dataSrc_extended = f'{dataSrc}_lead{lead_time}_ens{ens_num}'

                        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
                        df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)

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
                            df_vectorPop_monthly_yrFiltered = df_vectorPop_monthly[~df_vectorPop_monthly.index.year.isin(exclYrs)]
                            df_den_monthly_yrFiltered = df_den_monthly[~df_den_monthly.index.year.isin(exclYrs)]
                            
                            # manually create the popFrac columns and add them to vectorPop df
                            popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
                            df_vectorPop_monthly_yrFiltered_withFrac = df_vectorPop_monthly_yrFiltered.copy(deep=True)
                            df_vectorPop_monthly_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(E)'] / popTotal_E
                            df_vectorPop_monthly_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(L)'] / popTotal_E
                            df_vectorPop_monthly_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(P)'] / popTotal_E

                            fName = 'violinplot_' + dataSrc_extended + '_' + 'TA-TW' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                            #plot_violinplt_TA_TW(df_vectorPop_monthly_yrFiltered, loc, fName)
                            fName = 'violinplot_' + dataSrc_extended + '_' + 'prcp-WH' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                            #plot_violinplt_prcp_WH(df_vectorPop_monthly_yrFiltered, loc, fName)
                            fName = 'violinplot_' + dataSrc_extended + '_' + 'dev' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                            #plot_violinplt_dev(df_vectorPop_monthly_yrFiltered, loc, fName)
                            fName = 'violinplot_' + dataSrc_extended + '_' + 'surv' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                            #plot_violinplt_surv(df_vectorPop_monthly_yrFiltered, loc, fName)
                            fName = 'violinplot_' + dataSrc_extended + '_' + 'dev-surv' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                            #plot_violinplt_dev_surv(df_vectorPop_monthly_yrFiltered, loc, fName)
                            fName = 'violinplot_' + dataSrc_extended + '_' + 'popFrac' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                            #plot_violinplt_popFrac(df_vectorPop_monthly_yrFiltered_withFrac, loc, fName)

                            if dengueDataYrs_noOutbreaks: # only do this if we're restricting years to no outbreak dengue data years (2007-2020 minus 2017 and 2019)

                                # create new df combining vectorPop and denInc dfs
                                df_vectorPop_denInc_monthly_yrFiltered = df_vectorPop_monthly_yrFiltered_withFrac.copy(deep=True)
                                df_vectorPop_denInc_monthly_yrFiltered['dengueInc'] = df_den_monthly_yrFiltered['dengueInc']

                                fName = 'violinplot_'  + dataSrc_extended + '_' + 'denInc' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                                #plot_violinplt_denInc(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)
                                fName = 'violinplot_'  + dataSrc_extended + '_' + 'popFrac-denInc' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                                #plot_violinplt_popFrac_denInc(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)
                                fName = 'violinplot_'  + dataSrc_extended + '_' + 'allVars' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                                plot_violinplt_allVars(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)
                                
            else: # if MERRA_IMERG
                dataSrc_extended = f'{dataSrc}'

                df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
                df_vectorPop_monthly = convert_df_vectorPop_dly2mon(df_vectorPop_dly)

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
                    df_vectorPop_monthly_yrFiltered = df_vectorPop_monthly[~df_vectorPop_monthly.index.year.isin(exclYrs)]
                    df_den_monthly_yrFiltered = df_den_monthly[~df_den_monthly.index.year.isin(exclYrs)]
                    
                    # manually create the popFrac columns and add them to vectorPop df
                    popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
                    df_vectorPop_monthly_yrFiltered_withFrac = df_vectorPop_monthly_yrFiltered.copy(deep=True)
                    df_vectorPop_monthly_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(E)'] / popTotal_E
                    df_vectorPop_monthly_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(L)'] / popTotal_E
                    df_vectorPop_monthly_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered['pop_A(P)'] / popTotal_E

                    fName = 'violinplot_'  + dataSrc_extended + '_' + 'TA-TW' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                    #plot_violinplt_TA_TW(df_vectorPop_monthly_yrFiltered, loc, fName)
                    fName = 'violinplot_'  + dataSrc_extended + '_' + 'prcp-WH' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                    #plot_violinplt_prcp_WH(df_vectorPop_monthly_yrFiltered, loc, fName)
                    fName = 'violinplot_'  + dataSrc_extended + '_' + 'dev' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                    #plot_violinplt_dev(df_vectorPop_monthly_yrFiltered, loc, fName)
                    fName = 'violinplot_'  + dataSrc_extended + '_' + 'surv' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                    #plot_violinplt_surv(df_vectorPop_monthly_yrFiltered, loc, fName)
                    fName = 'violinplot_'  + dataSrc_extended + '_' + 'dev-surv' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                    #plot_violinplt_dev_surv(df_vectorPop_monthly_yrFiltered, loc, fName)
                    fName = 'violinplot_'  + dataSrc_extended + '_' + 'popFrac' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                    #plot_violinplt_popFrac(df_vectorPop_monthly_yrFiltered_withFrac, loc, fName)

                    if dengueDataYrs_noOutbreaks: # only do this if we're restricting years to no outbreak dengue data years (2007-2020 minus 2017 and 2019)

                        # create new df combining vectorPop and denInc dfs
                        df_vectorPop_denInc_monthly_yrFiltered = df_vectorPop_monthly_yrFiltered_withFrac.copy(deep=True)
                        df_vectorPop_denInc_monthly_yrFiltered['dengueInc'] = df_den_monthly_yrFiltered['dengueInc']

                        fName = 'violinplot_'  + dataSrc_extended + '_' + 'denInc' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                        #plot_violinplt_denInc(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)
                        fName = 'violinplot_'  + dataSrc_extended + '_' + 'popFrac-denInc' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                        #plot_violinplt_popFrac_denInc(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)
                        fName = 'violinplot_'  + dataSrc_extended + '_' + 'allVars' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                        plot_violinplt_allVars(df_vectorPop_denInc_monthly_yrFiltered, loc, fName)




# make violin plots for one month of ALL locations
if False:
    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1)) # note: in the below code we only make the plot for June and December
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True]#[True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]


    # load dengue incidence data
    file_dengue  = '/mnt/redwood/local_drive/chanud/RawData/epiReports/' + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly_allLocs = pd.read_csv(file_dengue, index_col='Date', parse_dates=True)
    df_den_monthly_allLocs = df_den_monthly_allLocs[df_den_monthly_allLocs.index.year.isin(years)] # cut extraneous years

    series_den_monthly_Jaffna = df_den_monthly_allLocs['Jaffna']
    series_den_monthly_Jaffna.name = 'dengueInc'
    df_den_monthly_Jaffna = series_den_monthly_Jaffna.to_frame()
    series_den_monthly_Negombo = df_den_monthly_allLocs['Gampaha']
    series_den_monthly_Negombo.name = 'dengueInc'
    df_den_monthly_Negombo = series_den_monthly_Negombo.to_frame()
    series_den_monthly_NuwaraEliya = df_den_monthly_allLocs['NuwaraEliya']
    series_den_monthly_NuwaraEliya.name = 'dengueInc'
    df_den_monthly_NuwaraEliya = series_den_monthly_NuwaraEliya.to_frame()

    # load model pipeline data at location of interest
    # also: convert from daily to monthly
    df_vectorPop_dly_Jaffna = get_df_vectorPop_dly_merged(dataSrc, 'Jaffna', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
    df_vectorPop_monthly_Jaffna = convert_df_vectorPop_dly2mon(df_vectorPop_dly_Jaffna)
    df_vectorPop_dly_Negombo = get_df_vectorPop_dly_merged(dataSrc, 'Negombo', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
    df_vectorPop_monthly_Negombo = convert_df_vectorPop_dly2mon(df_vectorPop_dly_Negombo)
    df_vectorPop_dly_NuwaraEliya = get_df_vectorPop_dly_merged(dataSrc, 'NuwaraEliya', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
    df_vectorPop_monthly_NuwaraEliya = convert_df_vectorPop_dly2mon(df_vectorPop_dly_NuwaraEliya)

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
        df_vectorPop_monthly_yrFiltered_Jaffna = df_vectorPop_monthly_Jaffna[~df_vectorPop_monthly_Jaffna.index.year.isin(exclYrs)]
        df_den_monthly_yrFiltered_Jaffna = df_den_monthly_Jaffna[~df_den_monthly_Jaffna.index.year.isin(exclYrs)]
        df_vectorPop_monthly_yrFiltered_Negombo = df_vectorPop_monthly_Negombo[~df_vectorPop_monthly_Negombo.index.year.isin(exclYrs)]
        df_den_monthly_yrFiltered_Negombo = df_den_monthly_Negombo[~df_den_monthly_Negombo.index.year.isin(exclYrs)]
        df_vectorPop_monthly_yrFiltered_NuwaraEliya = df_vectorPop_monthly_NuwaraEliya[~df_vectorPop_monthly_NuwaraEliya.index.year.isin(exclYrs)]
        df_den_monthly_yrFiltered_NuwaraEliya = df_den_monthly_NuwaraEliya[~df_den_monthly_NuwaraEliya.index.year.isin(exclYrs)]

        # manually create the popFrac columns and add them to vectorPop df
        popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
        df_vectorPop_monthly_yrFiltered_withFrac_Jaffna = df_vectorPop_monthly_yrFiltered_Jaffna.copy(deep=True)
        df_vectorPop_monthly_yrFiltered_withFrac_Jaffna['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered_Jaffna['pop_A(E)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_Jaffna['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered_Jaffna['pop_A(L)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_Jaffna['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered_Jaffna['pop_A(P)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_Negombo = df_vectorPop_monthly_yrFiltered_Negombo.copy(deep=True)
        df_vectorPop_monthly_yrFiltered_withFrac_Negombo['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered_Negombo['pop_A(E)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_Negombo['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered_Negombo['pop_A(L)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_Negombo['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered_Negombo['pop_A(P)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_NuwaraEliya = df_vectorPop_monthly_yrFiltered_NuwaraEliya.copy(deep=True)
        df_vectorPop_monthly_yrFiltered_withFrac_NuwaraEliya['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered_NuwaraEliya['pop_A(E)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_NuwaraEliya['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered_NuwaraEliya['pop_A(L)'] / popTotal_E
        df_vectorPop_monthly_yrFiltered_withFrac_NuwaraEliya['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered_NuwaraEliya['pop_A(P)'] / popTotal_E

        if dengueDataYrs_noOutbreaks: # only do this if we're restricting years to no outbreak dengue data years (2007-2020 minus 2017 and 2019)

            # create new df combining vectorPop and denInc dfs
            df_vectorPop_denInc_monthly_yrFiltered_Jaffna = df_vectorPop_monthly_yrFiltered_withFrac_Jaffna.copy(deep=True)
            df_vectorPop_denInc_monthly_yrFiltered_Jaffna['dengueInc'] = df_den_monthly_yrFiltered_Jaffna['dengueInc']
            df_vectorPop_denInc_monthly_yrFiltered_Negombo = df_vectorPop_monthly_yrFiltered_withFrac_Negombo.copy(deep=True)
            df_vectorPop_denInc_monthly_yrFiltered_Negombo['dengueInc'] = df_den_monthly_yrFiltered_Negombo['dengueInc']
            df_vectorPop_denInc_monthly_yrFiltered_NuwaraEliya = df_vectorPop_monthly_yrFiltered_withFrac_NuwaraEliya.copy(deep=True)
            df_vectorPop_denInc_monthly_yrFiltered_NuwaraEliya['dengueInc'] = df_den_monthly_yrFiltered_NuwaraEliya['dengueInc']

            for month in [6, 12]: # we're only doing plots for June and December
                fileStr_month = '_mon' + str(month).zfill(2)
                fName = 'violinplot_' + 'allVars' + '_' + 'allLocs' + fileStr_month + fileStr_dengueDataYrs_noOutbreaks
                plot_violinplt_allVars_allLocsOneMon(df_vectorPop_denInc_monthly_yrFiltered_Jaffna, df_vectorPop_denInc_monthly_yrFiltered_Negombo, 
                                                     df_vectorPop_denInc_monthly_yrFiltered_NuwaraEliya, month, fName)




# make violin plots for all months of ONE location comparing Focks vs Eisen
if False:
    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True]#[True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]


    # load dengue incidence data
    file_dengue  = '/mnt/redwood/local_drive/chanud/RawData/epiReports/' + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly_allLocs = pd.read_csv(file_dengue, index_col='Date', parse_dates=True)
    df_den_monthly_allLocs = df_den_monthly_allLocs[df_den_monthly_allLocs.index.year.isin(years)] # cut extraneous years

    for loc in locs:
        print(loc)
        
        # subset dengue data to location of interest
        if loc == 'Negombo':
            den_loc = 'Gampaha'
        elif loc in ['NuwaraEliya', 'Jaffna']:
            den_loc = loc
        else:
            raise SystemExit("Error: unknown loc.")
        series_den_monthly = df_den_monthly_allLocs[den_loc]
        series_den_monthly.name = 'dengueInc'
        df_den_monthly = series_den_monthly.to_frame()

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        df_vectorPop_dly_Eisen = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, 'Eisen')
        df_vectorPop_monthly_Eisen = convert_df_vectorPop_dly2mon(df_vectorPop_dly_Eisen)
        df_vectorPop_dly_Focks = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, 'Focks')
        df_vectorPop_monthly_Focks = convert_df_vectorPop_dly2mon(df_vectorPop_dly_Focks)

        for dengueDataYrs_noOutbreaks in dengueDataYrs_noOutbreaks_opts:
            print(' dengueDataYrs_noOutbreaks: ' + str(dengueDataYrs_noOutbreaks))
            if dengueDataYrs_noOutbreaks:
                fileStr_dengueDataYrs_noOutbreaks = '_' + str(dengueDataYrs[0]) + '-' + str(dengueDataYrs[-1]) + '_noOutbreaks' + '_FocksVsEisen'
                exclYrs = list(set([year for year in years if year not in dengueDataYrs] + outbreakYrs))  # years to exclude from plot
            else:
                fileStr_dengueDataYrs_noOutbreaks = '_' + str(years[0]) + '-' + str(years[-1]) + '_FocksVsEisen'
                exclYrs = []

            # remove data in outbreak years we want to exclude
            df_vectorPop_monthly_yrFiltered_Eisen = df_vectorPop_monthly_Eisen[~df_vectorPop_monthly_Eisen.index.year.isin(exclYrs)]
            df_vectorPop_monthly_yrFiltered_Focks = df_vectorPop_monthly_Focks[~df_vectorPop_monthly_Focks.index.year.isin(exclYrs)]
            df_den_monthly_yrFiltered = df_den_monthly[~df_den_monthly.index.year.isin(exclYrs)]

            # manually create the popFrac columns and add them to vectorPop df
            popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
            df_vectorPop_monthly_yrFiltered_withFrac_Eisen = df_vectorPop_monthly_yrFiltered_Eisen.copy(deep=True)
            df_vectorPop_monthly_yrFiltered_withFrac_Eisen['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered_Eisen['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_yrFiltered_withFrac_Eisen['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered_Eisen['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_yrFiltered_withFrac_Eisen['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered_Eisen['pop_A(P)'] / popTotal_E
            df_vectorPop_monthly_yrFiltered_withFrac_Focks = df_vectorPop_monthly_yrFiltered_Focks.copy(deep=True)
            df_vectorPop_monthly_yrFiltered_withFrac_Focks['pop_A(E)_frac'] = df_vectorPop_monthly_yrFiltered_Focks['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_yrFiltered_withFrac_Focks['pop_A(L)_frac'] = df_vectorPop_monthly_yrFiltered_Focks['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_yrFiltered_withFrac_Focks['pop_A(P)_frac'] = df_vectorPop_monthly_yrFiltered_Focks['pop_A(P)'] / popTotal_E

            if dengueDataYrs_noOutbreaks: # only do this if we're restricting years to no outbreak dengue data years (2007-2020 minus 2017 and 2019)

                # create new df combining vectorPop and denInc dfs
                df_vectorPop_denInc_monthly_yrFiltered_Eisen = df_vectorPop_monthly_yrFiltered_withFrac_Eisen.copy(deep=True)
                df_vectorPop_denInc_monthly_yrFiltered_Eisen['dengueInc'] = df_den_monthly_yrFiltered['dengueInc']
                df_vectorPop_denInc_monthly_yrFiltered_Focks = df_vectorPop_monthly_yrFiltered_withFrac_Focks.copy(deep=True)
                df_vectorPop_denInc_monthly_yrFiltered_Focks['dengueInc'] = df_den_monthly_yrFiltered['dengueInc']

                fName = 'violinplot_' + 'allVars' + '_' + loc + fileStr_dengueDataYrs_noOutbreaks
                plot_violinplt_FocksVsEisen(df_vectorPop_denInc_monthly_yrFiltered_Eisen, df_vectorPop_denInc_monthly_yrFiltered_Focks, loc, fName)




# make scatter plots of crosstab tercile-tercile values
if False:

    # parameters related to modeling pipeline
    dataSrcs   = ['MERRA_IMERG', 'GEOS']
    lead_time  = 0             # only relevant for GEOS
    ens_nums   = range(1, 5+1) # only relevant for GEOS
    locs       = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years      = list(range(2001, 2020+1))
    months     = list(range(1, 12+1))
    intv       = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'
    noOutbreaks_opts = [True, False]  # including vs not including outbreak years
    outbreakYrs = [2017, 2019]
    
    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)
    
    # parameters related to the crosstab element of interest
    # format is 'terc1-terc2' where terc1 is tercile of the variable, terc2 is tercile of dengue
    crosstab_elts = ['HI-HI', 'MED-HI', 'LO-HI', 'HI-LO', 'MED-LO', 'LO-LO']
    map_tercAbbr_to_tercName = {'HI': 'High', 'MED': 'Medium', 'LO': 'Low'} # map to names of rows/columns in crosstab csv file

    if survCurve == 'Eisen':
        fileStr_survCurve = '_eisenQuadFit'
    else: # if survCurve is 'Focks'
        fileStr_survCurve = ''

    for dataSrc in dataSrcs:
        print(dataSrc)
        if dataSrc == 'MERRA_IMERG':
            dataSrc_extended = f'{dataSrc}'
        elif dataSrc == 'GEOS':
            dataSrc_extended = f'{dataSrc}_lead{lead_time}_ensMean'
        for crosstab_elt in crosstab_elts:
            print(f'Crosstab element: {crosstab_elt}')

            tercileStr_modelVar, tercileStr_den = crosstab_elt.split('-')
            dir_crosstabs = f'crosstabs_{tercileStr_modelVar}_{tercileStr_den}'
            tercileName_modelVar = map_tercAbbr_to_tercName.get(tercileStr_modelVar) # get corresponding row name in crosstabs csv file
            tercileName_den      = map_tercAbbr_to_tercName.get(tercileStr_den)      # get corresponding col name in crosstabs csv file
            colName_4_df = f'{tercileName_modelVar}-{tercileName_den} Value'         # (e.g., 'High-High Value')

            for loc in locs:
                print(loc)
                
                df_list = []
                df_name_list = []
                for i_lag, lagTime in enumerate(lagTimes):
                    print(f' Lag (mo): {lagTime}')
                    for noOutbreaks in noOutbreaks_opts:
                        print(f'   noOutbreaks: {noOutbreaks}')
                        if noOutbreaks:
                            fileStr_noOutbreaks = '_noOutbreaks'
                            exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
                        else:
                            fileStr_noOutbreaks = ''
                            exclYrs = []
                        
                        den_colName = f'{den_colName_base}-lag{lagTime}' # this must match col name format set in create_df_withLags()
                        den_pltName = den_colName # just reusing the column name here for simplicity
                        
                        # read in crosstab elt values
                        fileName_in = f'crosstab_tercile_splitbyMon_{dataSrc_extended}_{tercileStr_modelVar}-{tercileStr_den}_{den_pltName}_{loc}{fileStr_noOutbreaks}{fileStr_survCurve}'
                        filePath_norm = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc, dir_crosstabs, f'{fileName_in}_norm.csv')
                        df_crosstab_elt_norm_oneScenario = pd.read_csv(filePath_norm, index_col=0, header=0)

                        # rename based on lag+outbreakYr scenario and add to list
                        df_name = f'lag{lagTime}{fileStr_noOutbreaks}'
                        df_crosstab_elt_norm_oneScenario.rename(columns={colName_4_df: df_name}, inplace=True)
                        df_list.append(df_crosstab_elt_norm_oneScenario)

                df_crosstab_elt_norm = pd.concat(df_list, axis=1)

                fileName_out = f'crosstab_tercile_splitbyMon_{dataSrc_extended}_{tercileStr_modelVar}-{tercileStr_den}_{den_colName_base}_{loc}{fileStr_survCurve}_norm'
                plot_scatter_crosstabs(df_crosstab_elt_norm, loc, tercileName_modelVar, tercileName_den, fileName_out)




# make scatter plots of crosstab tercile-tercile values comparing Focks vs Eisen
if False:

    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs    = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years   = list(range(2001, 2020+1))
    months  = list(range(1, 12+1))
    intv    = 'intv00'
    initDays   = list(range(1, 7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'
    noOutbreaks_opts = [True, False]  # including vs not including outbreak years
    outbreakYrs = [2017, 2019]
    
    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)
    
    # parameters related to the crosstab element of interest
    # format is 'terc1-terc2' where terc1 is tercile of the variable, terc2 is tercile of dengue
    crosstab_elts = ['HI-HI', 'MED-HI', 'LO-HI', 'HI-LO', 'MED-LO', 'LO-LO']
    map_tercAbbr_to_tercName = {'HI': 'High', 'MED': 'Medium', 'LO': 'Low'} # map to names of rows/columns in crosstab csv file

    for crosstab_elt in crosstab_elts:
        print('Crosstab element: ' + crosstab_elt)

        tercileStr_modelVar, tercileStr_den = crosstab_elt.split('-')
        dir_crosstabs = 'crosstabs_' + tercileStr_modelVar + '_' + tercileStr_den + '/'
        tercileName_modelVar = map_tercAbbr_to_tercName.get(tercileStr_modelVar) # get corresponding row name in crosstabs csv file
        tercileName_den      = map_tercAbbr_to_tercName.get(tercileStr_den)      # get corresponding col name in crosstabs csv file
        colName_4_df = tercileName_modelVar + '-' + tercileName_den + ' Value' # (e.g., 'High-High Value')

        for loc in locs:
            print(loc)
            
            df_list_Eisen = []
            df_name_list_Eisen = []
            df_list_Focks = []
            df_name_list_Focks = []
            for i_lag, lagTime in enumerate(lagTimes):
                print(' Lag (mo): ' + str(lagTime))
                    
                for noOutbreaks in noOutbreaks_opts:
                    print('   noOutbreaks: ' + str(noOutbreaks))
                    if noOutbreaks:
                        fileStr_noOutbreaks = '_noOutbreaks'
                        exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
                    else:
                        fileStr_noOutbreaks = ''
                        exclYrs = []
                    
                    den_colName = den_colName_base + '-lag' + str(lagTime) # this must match col name format set in create_df_withLags()
                    den_pltName = den_colName # just reusing the column name here for simplicity
                        
                    # read in crosstab hi-hi values
                    fileName_in_Eisen = 'crosstab_tercile_splitbyMon_' + tercileStr_modelVar + '-' + tercileStr_den + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks + '_eisenQuadFit'
                    filePath_norm_Eisen = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc, dir_crosstabs, f'{fileName_in_Eisen}_norm.csv')
                    df_crosstab_elt_norm_oneScenario_Eisen = pd.read_csv(filePath_norm_Eisen, index_col=0, header=0)
                    fileName_in_Focks = 'crosstab_tercile_splitbyMon_' + tercileStr_modelVar + '-' + tercileStr_den + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks
                    filePath_norm_Focks = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc, dir_crosstabs, f'{fileName_in_Focks}_norm.csv')
                    df_crosstab_elt_norm_oneScenario_Focks = pd.read_csv(filePath_norm_Focks, index_col=0, header=0)

                    # rename based on lag+outbreakYr scenario and add to list
                    df_name_Eisen = 'lag' + str(lagTime) + fileStr_noOutbreaks
                    df_crosstab_elt_norm_oneScenario_Eisen.rename(columns={colName_4_df: df_name_Eisen}, inplace=True)
                    df_list_Eisen.append(df_crosstab_elt_norm_oneScenario_Eisen)
                    df_name_Focks = 'lag' + str(lagTime) + fileStr_noOutbreaks
                    df_crosstab_elt_norm_oneScenario_Focks.rename(columns={colName_4_df: df_name_Focks}, inplace=True)
                    df_list_Focks.append(df_crosstab_elt_norm_oneScenario_Focks)

            df_crosstab_elt_norm_Eisen = pd.concat(df_list_Eisen, axis=1)
            df_crosstab_elt_norm_Focks = pd.concat(df_list_Focks, axis=1)

            fileName_out = 'crosstab_tercile_splitbyMon_' + tercileStr_modelVar + '-' + tercileStr_den + '_' + den_colName_base + '_' + loc + '_FocksVsEisen' + '_norm'
            plot_scatter_crosstabs_FocksVsEisen(df_crosstab_elt_norm_Eisen, df_crosstab_elt_norm_Focks, loc, tercileName_modelVar, tercileName_den, fileName_out)




# make tercile-tercile plots of model vars against dengue incidence
if False:
    
    # parameters related to modeling pipeline
    dataSrc = 'MERRA_IMERG'
    locs = ['Negombo', 'NuwaraEliya', 'Jaffna']
    years = list(range(2001, 2020+1))
    months = list(range(1,12+1))
    intv = 'intv00'
    initDays = list(range(1,7+1))
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    survCurve  = 'Focks' # ['Focks', 'Eisen']
    var_colNames = ['TA_min (C)', 'TA (C)', 'TA_max (C)', 'PRCP (mm)', 'TW_min_nanFill (C)', 'TW_nanFill (C)', 'TW_max_nanFill (C)', 'WH (mm)',
                   'dev_E', 'surv_E', 'dev_L', 'surv_L', 'dev_P', 'surv_P', 'pop_A(E)', 'pop_A(L)', 'pop_A(P)'] # subset of vars to actually plot
    var_pltNames = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]

    # parameters related to dengue data
    lagTimes = [0, 1, 2]             # in months
    den_colName_base = 'denInc'      # we will rename columns to this string + the lag time (original name is the district)
    noOutbreaks_opts = [True, False] # whether to include data from outbreak years
    outbreakYrs = [2017, 2019]

    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)

    # load dengue data
    fileDir_base = '/mnt/redwood/local_drive/chanud/RawData/epiReports/'
    file_dengue  = fileDir_base + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly = pd.read_csv(file_dengue, index_col="Date", parse_dates=True)

    for loc in locs:
        print(loc)
        if loc == 'Negombo':
            den_loc = 'Gampaha'
        elif loc in ['NuwaraEliya', 'Jaffna']:
            den_loc = loc
        else:
            raise SystemExit("Error: unknown loc.")

        # subset dengue data to location of interest and create df with lagged data (which also subsets to time of interest)
        series_den_monthly_loc = df_den_monthly[den_loc]
        df_den_monthly_loc_wLags = create_df_withLags(series_den_monthly_loc, lagTimes, den_colName_base, startDate, endDate)

        # load model pipeline data at location of interest
        # also: convert from daily to monthly, subset to time of interest
        df_vectorPop_dly = get_df_vectorPop_dly_merged(dataSrc, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
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

            for i_var, var_colName in enumerate(var_colNames):
                var_pltName = var_pltNames[i_var]
                print('  Var: ' + var_pltName)

                for i_lag, lagTime in enumerate(lagTimes):
                    print("   Lag (mo): " + str(lagTime))

                    den_colName = den_colName_base + '-lag' + str(lagTime) # this must match col name format set in create_df_withLags()
                    den_pltName = den_colName # just reusing the column name here for simplicity

                    fileName = 'barTerciles_' + var_pltName + '_' + den_pltName + '_' + loc + fileStr_noOutbreaks + fileStr_survCurve
                    plot_bar_varX_denInc(df_vectorPop_den_monthly_yrFiltered, var_colName, den_colName, fileName)







############################################
### Brand new plots for paper 3
############################################



## Scatter/line plots - comparing meteorology between GEOS and MERRA/IMERG
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
    var_colNames  = ['TA_max (C)', 'TA_min (C)', 'VPD (kPa)']# ['PRCP (mm)', 'TA_max (C)', 'TA_min (C)', 'TA (C)', 'VPD (kPa)']
    var_pltNames  = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]
    
    for loc in locs:
        print(loc)
        
        for i_var, var_colName in enumerate(var_colNames):
            print(f' Variable: {var_colName}')
            var_pltName = var_pltNames[i_var]
            for lead_time in lead_times:
                print(f'  Lead time: {lead_time}')

                # load model pipeline data at location of interest
                # also: convert from daily to monthly
                dataSrcs_extended = f'{dataSrc_main}_lead{lead_time}_VS_{dataSrc_ref}'

                dfs_vectorPop_monthly_GEOS_ens = []
                for ens_num in ens_nums:
                    df_vectorPop_dly_GEOS_ens     = get_df_vectorPop_dly_merged(dataSrc_main, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
                    df_vectorPop_monthly_GEOS_ens = convert_df_vectorPop_dly2mon(df_vectorPop_dly_GEOS_ens)
                    dfs_vectorPop_monthly_GEOS_ens.append(df_vectorPop_monthly_GEOS_ens)
                df_vectorPop_monthly_GEOS_ensMean = pd.concat(dfs_vectorPop_monthly_GEOS_ens).groupby(level=0).mean()

                df_vectorPop_dly_MERRAIMERG     = get_df_vectorPop_dly_merged(dataSrc_ref, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
                df_vectorPop_monthly_MERRAIMERG = convert_df_vectorPop_dly2mon(df_vectorPop_dly_MERRAIMERG)

                fName = f'compareMet_{dataSrcs_extended}_{var_pltName}_{loc}'
                plot_scatter_compareMet_tseries(dfs_vectorPop_monthly_GEOS_ens, df_vectorPop_monthly_GEOS_ensMean, df_vectorPop_monthly_MERRAIMERG, var_colName, loc, years, False, fName)

                fName = f'compareMetAnom_{dataSrcs_extended}_{var_pltName}_{loc}'
                plot_scatter_compareMet_tseries(dfs_vectorPop_monthly_GEOS_ens, df_vectorPop_monthly_GEOS_ensMean, df_vectorPop_monthly_MERRAIMERG, var_colName, loc, years, True, fName)




## Scatter plots - RPSS vs R, comparing meteorology between GEOS and MERRA/IMERG
if False:
    # parameters related to modeling pipeline
    dataSrc_main  = 'GEOS'
    dataSrc_ref   = 'MERRA_IMERG'
    lead_times    = [0]           # only relevant for GEOS
    locs          = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years         = list(range(2001, 2020+1))
    months        = list(range(1, 12+1))
    survCurves    = ['Focks'] # ['Focks', 'Eisen']
    var_colNames  = ['PRCP (mm)', 'TA_max (C)', 'TA_min (C)', 'TA (C)', 'VPD (kPa)']
                    #['PRCP (mm)', 'WH (mm)', 'TA_max (C)', 'TA_min (C)', 'TA (C)', 'TW_max (C)', 'TW_min (C)', 'TW (C)',
                    #'dev_E', 'dev_L', 'dev_P', 'surv_E', 'surv_L', 'surv_P', 'pop_A(E)_frac', 'pop_A(L)_frac', 'pop_A(P)_frac', 'VPD (kPa)']
    var_pltNames  = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]
    doAnom_opts   = [True]#[True, False] # whether to get statistics for data as anomalies (relative to monthly climatology) or as is

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [False]#[True, False] # whether to consider (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]

    # Define style per location
    loc_styles = {
        'Jaffna':       {'marker': 'D', 'color': 'purple'},    # purple diamond
        'Negombo':      {'marker': 'o', 'color': '#DAA520'},   # gold circle
        'NuwaraEliya':  {'marker': '^', 'color': 'blue'}       # blue triangle
    }

    # Gather data
    for var_pltName in var_pltNames:
        print(f'{var_pltName}')
        for lead_time in lead_times:
            print(f' Lead time: {lead_time}')
            for survCurve in survCurves:
                print(f'  Survival curve: {survCurve}')
                if survCurve == 'Eisen':
                    fileStr_survCurve = '_eisenQuadFit'
                else: # if survCurve is 'Focks'
                    fileStr_survCurve = ''
                for dengueDataYrs_noOutbreaks in dengueDataYrs_noOutbreaks_opts:
                    print(f'   dengueDataYrs_noOutbreaks: {dengueDataYrs_noOutbreaks}')
                    if dengueDataYrs_noOutbreaks:
                        fileStr_dengueDataYrs_noOutbreaks = '_' + str(dengueDataYrs[0]) + '-' + str(dengueDataYrs[-1]) + '_noOutbreaks' + fileStr_survCurve
                    else:
                        fileStr_dengueDataYrs_noOutbreaks = '_' + str(years[0]) + '-' + str(years[-1]) + fileStr_survCurve
                    for doAnom in doAnom_opts:
                        print(f'    doAnom: {doAnom}')
                        if doAnom:
                            fileStr_doAnom = 'Anom_'
                        else:
                            fileStr_doAnom = ''
                        dataSrcs_extended = f'{fileStr_doAnom}{dataSrc_main}_lead{lead_time}_VS_{dataSrc_ref}'
                        doColorBySzn = False
                        fileName = f'compareMet_RPSSvsR_{dataSrcs_extended}_{var_pltName}{fileStr_dengueDataYrs_noOutbreaks}'
                        plot_scatter_compareMet_RPSSvsR(dataSrc_main, dataSrc_ref, lead_time, var_pltName, locs, loc_styles, dengueDataYrs_noOutbreaks, survCurve, doAnom, doColorBySzn, fileName)
                        doColorBySzn = True
                        fileName = f'compareMet_RPSSvsR_{dataSrcs_extended}_{var_pltName}{fileStr_dengueDataYrs_noOutbreaks}_colorBySzn'
                        plot_scatter_compareMet_RPSSvsR(dataSrc_main, dataSrc_ref, lead_time, var_pltName, locs, loc_styles, dengueDataYrs_noOutbreaks, survCurve, doAnom, doColorBySzn, fileName)



## Scatter plots - observations vs ensemble, comparing meteorology between GEOS and MERRA/IMERG and checking RPSS/R
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
    var_colNames  = ['PRCP (mm)', 'TA (C)']
                    #['PRCP (mm)', 'WH (mm)', 'TA_max (C)', 'TA_min (C)', 'TA (C)', 'TW_max (C)', 'TW_min (C)', 'TW (C)',
                    # 'dev_E', 'dev_L', 'dev_P', 'surv_E', 'surv_L', 'surv_P', 'pop_A(E)_frac', 'pop_A(L)_frac', 'pop_A(P)_frac', 'VPD (kPa)']
    var_pltNames  = [strConvert_paren2hyphen(strConvert_parenUnitsRemove(x)) for x in var_colNames]
    doAnom_opts   = [True]#[True, False] # whether to compute statistics for data as anomalies (relative to monthly climatology) or as is

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [False]#[True, False] # whether to consider (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
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
                        fileName = f'scatter_obsVsEns_{dataSrcs_extended}_{var_pltName}_{loc}{fileStr_dengueDataYrs_noOutbreaks}'
                        plot_scatter_compareMet_obsVsEns(dfs_vectorPop_monthly_GEOS_ens_yrFiltered_withFrac, df_vectorPop_monthly_GEOS_ensMean_yrFiltered_withFrac, 
                                                         df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac, var_colName, fileName)




## Violin plots - comparing GEOS ensemble mean and MERRA/IMERG

# make violin plots for all months of ONE location (GEOS ens mean in color, MERRA/IMERG in gray)
if False:
    # parameters related to modeling pipeline
    dataSrc_main = 'GEOS'
    dataSrc_ref  = 'MERRA_IMERG'
    lead_times   = [0]           # only relevant for GEOS
    ens_nums     = range(1, 5+1) # only relevant for GEOS
    locs         = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years        = list(range(2001, 2020+1))
    months       = list(range(1, 12+1))
    intv         = 'intv00'
    initDays     = list(range(1, 7+1))
    initEggs     = 1000
    initLarvae   = 1000
    initPupae    = 1000
    survCurve    = 'Focks' # ['Focks', 'Eisen']

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True]#[True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]

    # load dengue incidence data
    file_dengue  = '/mnt/redwood/local_drive/chanud/RawData/epiReports/' + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly_allLocs = pd.read_csv(file_dengue, index_col='Date', parse_dates=True)
    df_den_monthly_allLocs = df_den_monthly_allLocs[df_den_monthly_allLocs.index.year.isin(years)] # cut extraneous years

    for loc in locs:
        print(loc)
        
        # subset dengue data to location of interest
        if loc == 'Negombo':
            den_loc = 'Gampaha'
        elif loc in ['NuwaraEliya', 'Jaffna']:
            den_loc = loc
        else:
            raise SystemExit("Error: unknown loc.")
        series_den_monthly = df_den_monthly_allLocs[den_loc]
        series_den_monthly.name = 'dengueInc'
        df_den_monthly = series_den_monthly.to_frame()

        for lead_time in lead_times:
            print(f'Lead time: {lead_time}')

            # load model pipeline data at location of interest
            # also: convert from daily to monthly
            dataSrcs_extended = f'{dataSrc_main}_lead{lead_time}_ensMean_VS_{dataSrc_ref}'

            df_vectorPop_dly_GEOS = get_df_vectorPop_dly_merged_ensMean(dataSrc_main, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_nums)
            df_vectorPop_monthly_GEOS = convert_df_vectorPop_dly2mon(df_vectorPop_dly_GEOS)

            df_vectorPop_dly_MERRAIMERG = get_df_vectorPop_dly_merged(dataSrc_ref, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
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
                df_vectorPop_monthly_GEOS_yrFiltered = df_vectorPop_monthly_GEOS[~df_vectorPop_monthly_GEOS.index.year.isin(exclYrs)]
                df_vectorPop_monthly_MERRAIMERG_yrFiltered = df_vectorPop_monthly_MERRAIMERG[~df_vectorPop_monthly_MERRAIMERG.index.year.isin(exclYrs)]
                df_den_monthly_yrFiltered = df_den_monthly[~df_den_monthly.index.year.isin(exclYrs)]
                
                # manually create the popFrac columns and add them to vectorPop df
                popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae
                df_vectorPop_monthly_GEOS_yrFiltered_withFrac = df_vectorPop_monthly_GEOS_yrFiltered.copy(deep=True)
                df_vectorPop_monthly_GEOS_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered['pop_A(E)'] / popTotal_E
                df_vectorPop_monthly_GEOS_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered['pop_A(L)'] / popTotal_E
                df_vectorPop_monthly_GEOS_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered['pop_A(P)'] / popTotal_E
                df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac = df_vectorPop_monthly_MERRAIMERG_yrFiltered.copy(deep=True)
                df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac['pop_A(E)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered['pop_A(E)'] / popTotal_E
                df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac['pop_A(L)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered['pop_A(L)'] / popTotal_E
                df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac['pop_A(P)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered['pop_A(P)'] / popTotal_E

                if dengueDataYrs_noOutbreaks: # only do this if we're restricting years to no outbreak dengue data years (2007-2020 minus 2017 and 2019)

                    # create new df combining vectorPop and denInc dfs
                    df_vectorPop_denInc_monthly_GEOS_yrFiltered = df_vectorPop_monthly_GEOS_yrFiltered_withFrac.copy(deep=True)
                    df_vectorPop_denInc_monthly_GEOS_yrFiltered['dengueInc'] = df_den_monthly_yrFiltered['dengueInc']
                    df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered = df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac.copy(deep=True)
                    df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered['dengueInc'] = df_den_monthly_yrFiltered['dengueInc']

                    fName = f'violinplot_{dataSrcs_extended}_allVars_{loc}{fileStr_dengueDataYrs_noOutbreaks}'
                    plot_violinplt_compareDatasets(df_vectorPop_denInc_monthly_GEOS_yrFiltered, dataSrc_main, 
                                                   df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered, dataSrc_ref,
                                                   loc, survCurve, fName)




# make violin plots for one month of ALL locations (GEOS ens mean in color, MERRA/IMERG in gray)
if False:
    # parameters related to modeling pipeline
    dataSrc_main = 'GEOS'
    dataSrc_ref  = 'MERRA_IMERG'
    lead_times   = [0]           # only relevant for GEOS
    ens_nums     = range(1, 5+1) # only relevant for GEOS
    locs         = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years        = list(range(2001, 2020+1))
    months       = list(range(1, 12+1))
    intv         = 'intv00'
    initDays     = list(range(1, 7+1))
    initEggs     = 1000
    initLarvae   = 1000
    initPupae    = 1000
    survCurve    = 'Focks' # ['Focks', 'Eisen']

    # parameters related to dengue data
    dengueDataYrs_noOutbreaks_opts = [True]#[True, False] # whether to plot (dengue data years w/o outbreak years (2007-20 w/o 2017, 2019)) or (all years)
    dengueDataYrs = list(range(2007, 2020+1))
    outbreakYrs = [2017, 2019]


    # load dengue incidence data for each location
    file_dengue  = '/mnt/redwood/local_drive/chanud/RawData/epiReports/' + 'dengueInc_monthly_2007_2022.csv'
    df_den_monthly_allLocs = pd.read_csv(file_dengue, index_col='Date', parse_dates=True)
    df_den_monthly_allLocs = df_den_monthly_allLocs[df_den_monthly_allLocs.index.year.isin(years)] # cut extraneous years

    series_den_monthly_Jaffna = df_den_monthly_allLocs['Jaffna']
    series_den_monthly_Jaffna.name = 'dengueInc'
    df_den_monthly_Jaffna = series_den_monthly_Jaffna.to_frame()
    series_den_monthly_Negombo = df_den_monthly_allLocs['Gampaha']
    series_den_monthly_Negombo.name = 'dengueInc'
    df_den_monthly_Negombo = series_den_monthly_Negombo.to_frame()
    series_den_monthly_NuwaraEliya = df_den_monthly_allLocs['NuwaraEliya']
    series_den_monthly_NuwaraEliya.name = 'dengueInc'
    df_den_monthly_NuwaraEliya = series_den_monthly_NuwaraEliya.to_frame()

    for lead_time in lead_times:
        print(f'Lead time: {lead_time}')

        # load model pipeline data
        # also: convert from daily to monthly
        dataSrcs_extended = f'{dataSrc_main}_lead{lead_time}_ensMean_VS_{dataSrc_ref}'

        df_vectorPop_dly_GEOS_Jaffna = get_df_vectorPop_dly_merged_ensMean(dataSrc_main, 'Jaffna', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_nums)
        df_vectorPop_monthly_GEOS_Jaffna = convert_df_vectorPop_dly2mon(df_vectorPop_dly_GEOS_Jaffna)
        df_vectorPop_dly_GEOS_Negombo = get_df_vectorPop_dly_merged_ensMean(dataSrc_main, 'Negombo', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_nums)
        df_vectorPop_monthly_GEOS_Negombo = convert_df_vectorPop_dly2mon(df_vectorPop_dly_GEOS_Negombo)
        df_vectorPop_dly_GEOS_NuwaraEliya = get_df_vectorPop_dly_merged_ensMean(dataSrc_main, 'NuwaraEliya', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_nums)
        df_vectorPop_monthly_GEOS_NuwaraEliya = convert_df_vectorPop_dly2mon(df_vectorPop_dly_GEOS_NuwaraEliya)

        df_vectorPop_dly_MERRAIMERG_Jaffna = get_df_vectorPop_dly_merged(dataSrc_ref, 'Jaffna', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
        df_vectorPop_monthly_MERRAIMERG_Jaffna = convert_df_vectorPop_dly2mon(df_vectorPop_dly_MERRAIMERG_Jaffna)
        df_vectorPop_dly_MERRAIMERG_Negombo = get_df_vectorPop_dly_merged(dataSrc_ref, 'Negombo', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
        df_vectorPop_monthly_MERRAIMERG_Negombo = convert_df_vectorPop_dly2mon(df_vectorPop_dly_MERRAIMERG_Negombo)
        df_vectorPop_dly_MERRAIMERG_NuwaraEliya = get_df_vectorPop_dly_merged(dataSrc_ref, 'NuwaraEliya', years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
        df_vectorPop_monthly_MERRAIMERG_NuwaraEliya = convert_df_vectorPop_dly2mon(df_vectorPop_dly_MERRAIMERG_NuwaraEliya)

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
            df_vectorPop_monthly_GEOS_yrFiltered_Jaffna = df_vectorPop_monthly_GEOS_Jaffna[~df_vectorPop_monthly_GEOS_Jaffna.index.year.isin(exclYrs)]
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_Jaffna = df_vectorPop_monthly_MERRAIMERG_Jaffna[~df_vectorPop_monthly_MERRAIMERG_Jaffna.index.year.isin(exclYrs)]
            df_den_monthly_yrFiltered_Jaffna = df_den_monthly_Jaffna[~df_den_monthly_Jaffna.index.year.isin(exclYrs)]

            df_vectorPop_monthly_GEOS_yrFiltered_Negombo = df_vectorPop_monthly_GEOS_Negombo[~df_vectorPop_monthly_GEOS_Negombo.index.year.isin(exclYrs)]
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_Negombo = df_vectorPop_monthly_MERRAIMERG_Negombo[~df_vectorPop_monthly_MERRAIMERG_Negombo.index.year.isin(exclYrs)]
            df_den_monthly_yrFiltered_Negombo = df_den_monthly_Negombo[~df_den_monthly_Negombo.index.year.isin(exclYrs)]

            df_vectorPop_monthly_GEOS_yrFiltered_NuwaraEliya = df_vectorPop_monthly_GEOS_NuwaraEliya[~df_vectorPop_monthly_GEOS_NuwaraEliya.index.year.isin(exclYrs)]
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_NuwaraEliya = df_vectorPop_monthly_MERRAIMERG_NuwaraEliya[~df_vectorPop_monthly_MERRAIMERG_NuwaraEliya.index.year.isin(exclYrs)]
            df_den_monthly_yrFiltered_NuwaraEliya = df_den_monthly_NuwaraEliya[~df_den_monthly_NuwaraEliya.index.year.isin(exclYrs)]


            # manually create the popFrac columns and add them to vectorPop df
            popTotal_E = len(initDays) * initEggs # this code only works for eggs AND larvae and pupae as long as initEggs = initLarvae = initPupae

            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Jaffna = df_vectorPop_monthly_GEOS_yrFiltered_Jaffna.copy(deep=True)
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Jaffna['pop_A(E)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_Jaffna['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Jaffna['pop_A(L)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_Jaffna['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Jaffna['pop_A(P)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_Jaffna['pop_A(P)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Jaffna = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Jaffna.copy(deep=True)
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Jaffna['pop_A(E)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Jaffna['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Jaffna['pop_A(L)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Jaffna['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Jaffna['pop_A(P)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Jaffna['pop_A(P)'] / popTotal_E

            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Negombo = df_vectorPop_monthly_GEOS_yrFiltered_Negombo.copy(deep=True)
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Negombo['pop_A(E)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_Negombo['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Negombo['pop_A(L)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_Negombo['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Negombo['pop_A(P)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_Negombo['pop_A(P)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Negombo = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Negombo.copy(deep=True)
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Negombo['pop_A(E)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Negombo['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Negombo['pop_A(L)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Negombo['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Negombo['pop_A(P)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_Negombo['pop_A(P)'] / popTotal_E

            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_NuwaraEliya = df_vectorPop_monthly_GEOS_yrFiltered_NuwaraEliya.copy(deep=True)
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_NuwaraEliya['pop_A(E)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_NuwaraEliya['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_NuwaraEliya['pop_A(L)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_NuwaraEliya['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_GEOS_yrFiltered_withFrac_NuwaraEliya['pop_A(P)_frac'] = df_vectorPop_monthly_GEOS_yrFiltered_NuwaraEliya['pop_A(P)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_NuwaraEliya = df_vectorPop_monthly_MERRAIMERG_yrFiltered_NuwaraEliya.copy(deep=True)
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_NuwaraEliya['pop_A(E)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_NuwaraEliya['pop_A(E)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_NuwaraEliya['pop_A(L)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_NuwaraEliya['pop_A(L)'] / popTotal_E
            df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_NuwaraEliya['pop_A(P)_frac'] = df_vectorPop_monthly_MERRAIMERG_yrFiltered_NuwaraEliya['pop_A(P)'] / popTotal_E

            if dengueDataYrs_noOutbreaks: # only do this if we're restricting years to no outbreak dengue data years (2007-2020 minus 2017 and 2019)

                # create new df combining vectorPop and denInc dfs
                df_vectorPop_denInc_monthly_GEOS_yrFiltered_Jaffna = df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Jaffna.copy(deep=True)
                df_vectorPop_denInc_monthly_GEOS_yrFiltered_Jaffna['dengueInc'] = df_den_monthly_yrFiltered_Jaffna['dengueInc']
                df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_Jaffna = df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Jaffna.copy(deep=True)
                df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_Jaffna['dengueInc'] = df_den_monthly_yrFiltered_Jaffna['dengueInc']

                df_vectorPop_denInc_monthly_GEOS_yrFiltered_Negombo = df_vectorPop_monthly_GEOS_yrFiltered_withFrac_Negombo.copy(deep=True)
                df_vectorPop_denInc_monthly_GEOS_yrFiltered_Negombo['dengueInc'] = df_den_monthly_yrFiltered_Negombo['dengueInc']
                df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_Negombo = df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_Negombo.copy(deep=True)
                df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_Negombo['dengueInc'] = df_den_monthly_yrFiltered_Negombo['dengueInc']

                df_vectorPop_denInc_monthly_GEOS_yrFiltered_NuwaraEliya = df_vectorPop_monthly_GEOS_yrFiltered_withFrac_NuwaraEliya.copy(deep=True)
                df_vectorPop_denInc_monthly_GEOS_yrFiltered_NuwaraEliya['dengueInc'] = df_den_monthly_yrFiltered_NuwaraEliya['dengueInc']
                df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_NuwaraEliya = df_vectorPop_monthly_MERRAIMERG_yrFiltered_withFrac_NuwaraEliya.copy(deep=True)
                df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_NuwaraEliya['dengueInc'] = df_den_monthly_yrFiltered_NuwaraEliya['dengueInc']

                for month in [6, 12]:
                    print(f'  Month: {month}')
                    fileStr_month = 'mon' + str(month).zfill(2)
                    fileName = f'violinplot_{dataSrcs_extended}_allVars_allLocs_{fileStr_month}{fileStr_dengueDataYrs_noOutbreaks}'
                    plot_violinplt_compareDatasets_allLocsOneMon(df_vectorPop_denInc_monthly_GEOS_yrFiltered_Jaffna, df_vectorPop_denInc_monthly_GEOS_yrFiltered_Negombo,
                                                                 df_vectorPop_denInc_monthly_GEOS_yrFiltered_NuwaraEliya, dataSrc_main,
                                                                 df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_Jaffna, df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_Negombo, 
                                                                 df_vectorPop_denInc_monthly_MERRAIMERG_yrFiltered_NuwaraEliya, dataSrc_ref, month, survCurve, fileName)




## Crosstab scatter plots - comparing GEOS ensemble mean and MERRA/IMERG

# make scatter plots of crosstab tercile-tercile values (GEOS ens mean in color, MERRA/IMERG in gray)
if False:

    # parameters related to modeling pipeline
    dataSrc_main = 'GEOS'
    dataSrc_ref  = 'MERRA_IMERG'
    lead_time    = 0             # only relevant for GEOS
    ens_nums     = range(1, 5+1) # only relevant for GEOS
    locs         = ['Negombo', 'Jaffna', 'NuwaraEliya']
    years        = list(range(2001, 2020+1))
    months       = list(range(1, 12+1))
    intv         = 'intv00'
    initDays     = list(range(1, 7+1))
    initEggs     = 1000
    initLarvae   = 1000
    initPupae    = 1000

    # parameters related to dengue data
    lagTimes = [0, 1, 2]              # in months
    den_colName_base = 'denInc'
    noOutbreaks_opts = [True, False]  # including vs not including outbreak years
    outbreakYrs = [2017, 2019]
    
    # parameters related to modeling pipeline + dengue data
    startDate = '2007-01-01'  # start date (first year of dengue data)
    endDate   = '2020-12-31'  # end date (last year of downloaded MERRA/IMERG data)
    
    # parameters related to the crosstab element of interest
    # format is 'terc1-terc2' where terc1 is tercile of the variable, terc2 is tercile of dengue
    crosstab_elts = ['HI-HI', 'MED-HI', 'LO-HI', 'HI-LO', 'MED-LO', 'LO-LO']
    map_tercAbbr_to_tercName = {'HI': 'High', 'MED': 'Medium', 'LO': 'Low'} # map to names of rows/columns in crosstab csv file

    dataSrc_extended_main = f'{dataSrc_main}_lead{lead_time}_ensMean'
    dataSrc_extended_ref  = f'{dataSrc_ref}'
    dataSrcs_extended = f'{dataSrc_main}_lead{lead_time}_ensMean_VS_{dataSrc_ref}'
    for crosstab_elt in crosstab_elts:
        print('Crosstab element: ' + crosstab_elt)

        tercileStr_modelVar, tercileStr_den = crosstab_elt.split('-')
        dir_crosstabs = f'crosstabs_{tercileStr_modelVar}_{tercileStr_den}'
        tercileName_modelVar = map_tercAbbr_to_tercName.get(tercileStr_modelVar) # get corresponding row name in crosstabs csv file
        tercileName_den      = map_tercAbbr_to_tercName.get(tercileStr_den)      # get corresponding col name in crosstabs csv file
        colName_4_df = f'{tercileName_modelVar}-{tercileName_den} Value'         # (e.g., 'High-High Value')

        for loc in locs:
            print(loc)
            
            df_list_main = []
            df_name_list_main = []
            df_list_ref = []
            df_name_list_ref = []
            for i_lag, lagTime in enumerate(lagTimes):
                print(' Lag (mo): ' + str(lagTime))
                    
                for noOutbreaks in noOutbreaks_opts:
                    print('   noOutbreaks: ' + str(noOutbreaks))
                    if noOutbreaks:
                        fileStr_noOutbreaks = '_noOutbreaks'
                        exclYrs = outbreakYrs  # outbreak years (to exclude from plot)
                    else:
                        fileStr_noOutbreaks = ''
                        exclYrs = []
                    
                    den_colName = f'{den_colName_base}-lag{lagTime}' # this must match col name format set in create_df_withLags()
                    den_pltName = den_colName # just reusing the column name here for simplicity
                        
                    # read in crosstab hi-hi values
                    fileName_in_main   = f'crosstab_tercile_splitbyMon_{dataSrc_extended_main}_{tercileStr_modelVar}-{tercileStr_den}_{den_pltName}_{loc}{fileStr_noOutbreaks}'
                    filePath_norm_main = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc_main, dir_crosstabs, f'{fileName_in_main}_norm.csv')
                    df_crosstab_elt_norm_oneScenario_main = pd.read_csv(filePath_norm_main, index_col=0, header=0)
                    fileName_in_ref   = f'crosstab_tercile_splitbyMon_{dataSrc_extended_ref}_{tercileStr_modelVar}-{tercileStr_den}_{den_pltName}_{loc}{fileStr_noOutbreaks}'
                    filePath_norm_ref = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_PROCESSED_BASE, dataSrc_ref, dir_crosstabs, f'{fileName_in_ref}_norm.csv')
                    df_crosstab_elt_norm_oneScenario_ref = pd.read_csv(filePath_norm_ref, index_col=0, header=0)

                    # rename based on lag+outbreakYr scenario and add to list
                    df_name = f'lag{lagTime}{fileStr_noOutbreaks}'
                    df_crosstab_elt_norm_oneScenario_main.rename(columns={colName_4_df: df_name}, inplace=True)
                    df_list_main.append(df_crosstab_elt_norm_oneScenario_main)
                    df_crosstab_elt_norm_oneScenario_ref.rename(columns={colName_4_df: df_name}, inplace=True)
                    df_list_ref.append(df_crosstab_elt_norm_oneScenario_ref)

            df_crosstab_elt_norm_main = pd.concat(df_list_main, axis=1)
            df_crosstab_elt_norm_ref  = pd.concat(df_list_ref, axis=1)

            fileName_out = f'crosstab_tercile_splitbyMon_{dataSrcs_extended}_{tercileStr_modelVar}-{tercileStr_den}_{den_colName_base}_{loc}_norm'
            plot_scatter_crosstabs_compareDatasets(df_crosstab_elt_norm_main, df_crosstab_elt_norm_ref, loc, tercileName_modelVar, tercileName_den, fileName_out)




## Climatology line plots - comparing GEOS ensemble mean and MERRA/IMERG monthly climatologies
if False:
    
    # parameters related to modeling pipeline
    dataSrc_main  = 'GEOS'
    dataSrc_ref   = 'MERRA_IMERG'
    lead_time     = 0           # only relevant for GEOS
    ens_nums      = range(1, 5+1) # only relevant for GEOS
    locs          = ['Jaffna', 'Negombo', 'NuwaraEliya']
    years         = list(range(2001, 2020+1))
    months        = list(range(1, 12+1))
    intv          = 'intv00'
    initDays      = list(range(1, 7+1))
    initEggs      = 1000
    initLarvae    = 1000
    initPupae     = 1000
    survCurve     = 'Focks' # ['Focks', 'Eisen']

    # parameters related to plots
    month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    plot_rows = ['PRCP (mm)', 'TA', 'VPD (kPa)']
    var_plotInfo = {
        'PRCP (mm)': {'fullName': 'Total precipitation (mm)', 'ylim': (0, 450)},
        'TA': {'fullName': 'Min / mean / max 2-m air temperature (°C)','ylim': (10, 34)},
        'VPD (kPa)': {'fullName': 'Vapor pressure deficit (kPa)','ylim': (0, 1.5)}
    }

    font_size_colTitle    = 20
    font_size_monLabels   = 16
    font_size_yTicks      = 16
    font_size_yLabels     = 18
    font_size_legend      = 16

    fig, axes = plt.subplots(nrows=len(plot_rows), ncols=len(locs), figsize=(21, 16), 
                             sharex=True, gridspec_kw={'height_ratios': [1, 1.5, 1]})
    dataSrcs_extended = f'{dataSrc_main}_lead{lead_time}_VS_{dataSrc_ref}'

    for i_loc, loc in enumerate(locs):
        print(f'  Location: {loc}')

        # load model pipeline data at location of interest
        # also: convert from daily to monthly
        dfs_vectorPop_monthly_GEOS_ens = []
        for ens_num in ens_nums:
            df_vectorPop_dly_GEOS_ens     = get_df_vectorPop_dly_merged(dataSrc_main, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve, lead_time, ens_num)
            df_vectorPop_monthly_GEOS_ens = convert_df_vectorPop_dly2mon(df_vectorPop_dly_GEOS_ens)
            dfs_vectorPop_monthly_GEOS_ens.append(df_vectorPop_monthly_GEOS_ens)
        df_vectorPop_monthly_GEOS_ensMean = pd.concat(dfs_vectorPop_monthly_GEOS_ens).groupby(level=0).mean()

        df_vectorPop_dly_MERRAIMERG     = get_df_vectorPop_dly_merged(dataSrc_ref, loc, years, months, intv, initDays, initEggs, initLarvae, initPupae, survCurve)
        df_vectorPop_monthly_MERRAIMERG = convert_df_vectorPop_dly2mon(df_vectorPop_dly_MERRAIMERG)

        ## Copy all input dfs so we don't modify them
        df_GEOS_ensMean = df_vectorPop_monthly_GEOS_ensMean.copy()
        df_MERRAIMERG   = df_vectorPop_monthly_MERRAIMERG.copy()

        # Add month columns to support calculations
        df_GEOS_ensMean['Month'] = df_GEOS_ensMean.index.month
        df_MERRAIMERG['Month'] = df_MERRAIMERG.index.month

        # Compute climatologies for variables of interest
        df_GEOS_ensMean_monClim = df_GEOS_ensMean.groupby('Month')[['PRCP (mm)', 'TA_min (C)', 'TA (C)', 'TA_max (C)', 'VPD (kPa)']].mean()
        df_MERRAIMERG_monClim = df_MERRAIMERG.groupby('Month')[['PRCP (mm)', 'TA_min (C)', 'TA (C)', 'TA_max (C)', 'VPD (kPa)']].mean()

        for i_row, plot_row in enumerate(plot_rows):
            print(f' Plot row: {plot_row}')

            ax = axes[i_row, i_loc]

            # Plot precipitation
            if plot_row == 'PRCP (mm)':
                var_colName = 'PRCP (mm)'

                ax.plot(df_GEOS_ensMean_monClim[var_colName].index, df_GEOS_ensMean_monClim[var_colName], linewidth=6, color='blue', label='GEOS ensemble mean')
                ax.plot(df_MERRAIMERG_monClim[var_colName].index, df_MERRAIMERG_monClim[var_colName], linewidth=6, linestyle='--', color='orange', label='MERRA-2/IMERG')

            # Plot min/mean/max 2-m temperature
            elif plot_row == 'TA':
                var_colNames = ['TA_min (C)', 'TA (C)', 'TA_max (C)']

                for var_colName in var_colNames:
                    ax.plot(df_GEOS_ensMean_monClim[var_colName].index, df_GEOS_ensMean_monClim[var_colName], linewidth=5, color='blue')
                    ax.plot(df_MERRAIMERG_monClim[var_colName].index, df_MERRAIMERG_monClim[var_colName], linewidth=5, linestyle='--', color='orange')

            # Plot VPD
            elif plot_row == 'VPD (kPa)':
                var_colName = 'VPD (kPa)'

                ax.plot(df_GEOS_ensMean_monClim[var_colName].index, df_GEOS_ensMean_monClim[var_colName], linewidth=6, color='blue', label='GEOS ensemble mean')
                ax.plot(df_MERRAIMERG_monClim[var_colName].index, df_MERRAIMERG_monClim[var_colName], linewidth=6, linestyle='--', color='orange', label='MERRA-2/IMERG')

            if i_row == 0: # only add title to top row
                ax.set_title(loc, fontsize=font_size_colTitle, pad=15)

            ax.set_xlim(0.5, 12.5)
            ax.set_xticks(range(1, 12+1))
            if i_row == len(plot_rows) - 1: # only add x-axis tick labels to bottom row
                ax.set_xticklabels(month_labels, fontsize=font_size_monLabels)

            ax.set_ylim(var_plotInfo[plot_row]['ylim'])
            ax.tick_params(axis='y', labelsize=font_size_yTicks)
            if i_loc > 0: # only include y-axis tick labels on leftmost plot in each row
                ax.tick_params(axis='y', labelleft=False)

            ax.grid(True, alpha=0.3)

            axes[i_row, 0].set_ylabel(var_plotInfo[plot_row]['fullName'], fontsize=font_size_yLabels)  # Only put y-axis label on leftmost subplots

    axes[0, 0].legend(loc='upper left', fontsize=font_size_legend) # Put legend on first subplot

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.03) # reduce horizontal width between columns

    # Save the figure
    fName = f'compareMet_clim_{dataSrcs_extended}_allVars'
    dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/5_plotting/'
    os.makedirs(dir_out, exist_ok=True)
    plt.savefig(os.path.join(dir_out,  f'{fName}.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)




## Meteorology PDFs - comparing (GEOS and MERRA/IMERG) raw data and bias-corrected data with observations 
if False:

    regenerate_metData = False # whether to regenerate metData (True) or load from existing .pkl file (False)

    # Directories and filepaths
    pltDir = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/5_plotting/plots/plts_compareMet_pdfs/'
    fPath_metData = os.path.join(pltDir, 'metData_GEOS_MERRAIMERG_rawVSbcVSobs.pkl') # compiled met data, saved as .pkl to save time

    
    lead_time      = 0             # only relevant for GEOS
    ens_nums       = range(1, 5+1) # only relevant for GEOS
    varNames       = ['prcp', 'airTemp', 'soilTemp']
    
    # City metadata
    cities = {
        "Jaffna": {"lat": 9.6615, "lon": 80.0255, "elev": 5., "obs_file": PATH_WMO_JAFFNA},
        "Negombo": {"lat": 7.2008, "lon": 79.8737, "elev": 2., "obs_file": PATH_WMO_NEGOMBO},
        "NuwaraEliya": {"lat": 6.9497, "lon": 80.7891, "elev": 1868., "obs_file": PATH_WMO_NUWARAELIYA},
    }

    # If met data needs to be regenerated, do so and save it before plotting it.
    # If met data is already saved as .pkl, load it and plot it.

    if regenerate_metData:
        # Initialize data structure for storing all data
        metData = {}
        for varName in varNames:
            metData[varName] = {}
            for city in cities:
                metData[varName][city] = {
                    'obs': None,
                    'MERRAIMERG': {
                        'original': None,
                        'biasCorrected': None
                    },
                    'GEOS': {
                        'original': None,
                        'biasCorrected': None
                    }
                }

        # Load observation data
        for city, cityInfo in cities.items():
            df_obs_daily = load_obsData(cityInfo['obs_file'])
            metData['prcp'][city]['obs'] = df_obs_daily['PRCP'] # precipitation
            metData['airTemp'][city]['obs'] = df_obs_daily['TEMP'] # air temperature
            metData['soilTemp'][city]['obs'] = df_obs_daily['TEMP'] # soil temperature (use air temp observations)

        # MERRA/IMERG info
        fPath_IMERG_prcp_base,     _ = os.path.splitext(PATH_IMERG_PRECIP)
        fPath_MERRA_airTemp_base,  _ = os.path.splitext(PATH_MERRA_T)
        fPath_MERRA_soilTemp_base, _ = os.path.splitext(PATH_MERRA_TSOIL)
        MERRAIMERG_info = {
            'prcp': {'path': PATH_IMERG_PRECIP, 'path_base_bc': fPath_IMERG_prcp_base, 'var': 'precipitationCal', 'extract_func': extract_mainData_prcp},
            'airTemp': {'path': PATH_MERRA_T, 'path_base_bc': fPath_MERRA_airTemp_base, 'var': 'T2M', 'extract_func': extract_mainData_temp},
            'soilTemp': {'path': PATH_MERRA_TSOIL, 'path_base_bc': fPath_MERRA_soilTemp_base, 'var': 'TSOIL1', 'extract_func': extract_mainData_temp}
        }

        # Load MERRA/IMERG original data
        for varName, info in MERRAIMERG_info.items():
            ds = xr.open_dataset(info['path'])

            data = extract_city_data(ds, info['extract_func'], info['var'], cities, 'MERRAIMERG')
            for city in cities:
                metData[varName][city]['MERRAIMERG']['original'] = data[city]
            ds.close()

        # Load MERRA/IMERG bias-corrected data
        # Even though bc data is already city-specific, use extract_city_data() to apply unit conversion and daily aggregation
        for varName, info in MERRAIMERG_info.items():
            for city in cities:
                if varName == 'prcp':
                    fPath = f'{info["path_base_bc"]}_{city}_biasCorrected.nc'
                else:
                    fPath = f'{info["path_base_bc"]}.{city}.biasCorrected.nc'
                ds = xr.open_dataset(fPath)

                data = extract_city_data(ds, info['extract_func'], info['var'], cities, 'MERRAIMERG', city=city)
                metData[varName][city]['MERRAIMERG']['biasCorrected'] = data[city]
                ds.close()

        # Load GEOS original data
        GEOS_info = {
            'prcp': {'var': 'PRECTOT', 'extract_func': extract_mainData_prcp},
            'airTemp': {'var': 'T2M', 'extract_func': extract_mainData_temp},
            'soilTemp': {'var': 'TS', 'extract_func': extract_mainData_temp}
        }

        for varName, info in GEOS_info.items():
            ds_ens = []
            
            for ens_num in ens_nums:
                fPath = os.path.join(DIR_GEOS, f'{FILE_GEOS_PREFIX}_lead{lead_time}_ens{ens_num}_{info["var"]}.nc')
                ds_ens.append(xr.open_dataset(fPath))
            ds_ensMean = xr.concat(ds_ens, dim='ensemble').mean(dim='ensemble') # Calculate GEOS ensemble mean
            
            # Extract location-specific data from spatial grid
            data = extract_city_data(ds_ensMean, info['extract_func'], info['var'], cities, 'GEOS')
            
            for city in cities:
                metData[varName][city]['GEOS']['original'] = data[city]
            for ds in ds_ens:
                ds.close()

        # Load GEOS bias-corrected data
        # Even though bc data is already city-specific, use extract_city_data() to apply unit conversion and daily aggregation
        for varName, info in GEOS_info.items():
            for city in cities:
                ds_ens = []

                for ens_num in ens_nums:
                    fPath = os.path.join(DIR_GEOS, f'{FILE_GEOS_PREFIX}_lead{lead_time}_ens{ens_num}_{info["var"]}_{city}_biasCorrected.nc')
                    ds_ens.append(xr.open_dataset(fPath))
                ds_ensMean = xr.concat(ds_ens, dim='ensemble').mean(dim='ensemble') # Calculate GEOS ensemble mean

                data = extract_city_data(ds_ensMean, info['extract_func'], info['var'], cities, 'GEOS', city=city)
                metData[varName][city]['GEOS']['biasCorrected'] = data[city]

                for ds in ds_ens:
                    ds.close()

        # Save meteorology data as .pkl to save time in future
        with open(fPath_metData, 'wb') as f:
            pickle.dump(metData, f)
        print(f'Saved metData to: {fPath_metData}')

    else:

        # Load previously saved metData
        with open(fPath_metData, 'rb') as f:
            metData = pickle.load(f)
        print(f'Loaded metData from: {fPath_metData}')

    # Generate panel plot showing PDFs for all variables, cities, and data sources
    fileName = f'compareMet_pdf_GEOS_MERRAIMERG_rawVSbcVSobs'
    plot_pdf_compareMet_9panel(metData, varNames, cities, pltDir, fileName)


    # Assess number of small precip values for each dataset.
    # We do this because the above plotting showed GEOS has NO zero values in all locations.
    # This is perhaps because there are very small, nearly zero values for GEOS which we should
    # count as zeroes, so we need to find a threshold to call zero values.
    for city in cities:

        print(f'\n{"="*70}')
        print(f'{city}')
        print(f'{"="*70}')

        for dataSrc, label in [
            ('GEOS', 'GEOS original'),
            ('GEOS bc', 'GEOS bias corrected'),
            ('MERRAIMERG', 'MERRA/IMERG original'),
            ('MERRAIMERG bc', 'MERRA/IMERG bias corrected'),
            ('obs', 'Observations')
        ]:

            if dataSrc == 'GEOS':
                data = metData['prcp'][city]['GEOS']['original']
            elif dataSrc == 'GEOS bc':
                data = metData['prcp'][city]['GEOS']['biasCorrected']
            elif dataSrc == 'MERRAIMERG':
                data = metData['prcp'][city]['MERRAIMERG']['original']
            elif dataSrc == 'MERRAIMERG bc':
                data = metData['prcp'][city]['MERRAIMERG']['biasCorrected']
            elif dataSrc == 'obs':
                data = metData['prcp'][city]['obs']

            data = data.dropna()

            print(f'\n{label}:')
            print(f'  exactly 0 mm: {(data == 0).sum()} ({(data == 0).mean()*100:.2f}%)')
            print(f'  < 0.001 mm:   {(data < 0.001).sum()} ({(data < 0.001).mean()*100:.2f}%)')
            print(f'  < 0.01 mm:    {(data < 0.01).sum()} ({(data < 0.01).mean()*100:.2f}%)')
            print(f'  < 0.1 mm:     {(data < 0.1).sum()} ({(data < 0.1).mean()*100:.2f}%)')
            print(f'  < 1.0 mm:     {(data < 1.0).sum()} ({(data < 1.0).mean()*100:.2f}%)')