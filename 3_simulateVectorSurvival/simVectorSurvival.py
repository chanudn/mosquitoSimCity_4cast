#*******************************************************
# simVectorSurvival.py
#*******************************************************
# Revision history

# Import packages
import pandas as pd
import numpy as np
import random
import sys
import os

# Import custom functions and constants
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants") # directory of constants file
from CONST_FILE_DIR import DIR_WHATCHEM_OUTPUT_DATA_BASE
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/4_data_analysis") # directory of helper fns file
from fns_helper import WHATCHEM_specs_to_fpath, read_WHATCHEM_output, convert_hrly2dlyMeanMinMax, fillMissingVals_linReg


## calc_devRate()
#  Calculates hourly development rate for a mosquito immature given water temperature and life stage.
#  Args:
#    waterTemp_hrly - hourly water temperature [deg C]
#    lifeStage      - life stage of immatures (options: "Eggs", "Larvae", "Pupae")
#  Output:
#    r              - hourly development rate
def calc_devRate(waterTemp_hrly, lifeStage):

    # Convert temperature from deg C to K
    waterTemp_hrly = waterTemp_hrly + 273.15

    # Define universal gas constant
    R_const = 1.98720425864083 # [cal / K mol]

    # Define dataframe of development rate parameters.
    # These parameters are based on a temperature-dependent enzyme-kinetics developmental rate model.
    #
    # For more information, see:
    #   Table S1 of Magori et al. 2009 (Skeeter Buster paper)
    #   Table 2 and Eqns 4, 5 of Focks et al. 1993
    #
    # Parameters:
    #   rho_25C - Development rate at 25 deg C, assuming no temperature inactivation of the critical enzyme [hr^-1]
    #   delH_A  - Enthalpy of activation of the reaction catalyzed by the enzyme [cal/mol]
    #   delH_H  - Enthalpy change associated with high temperature inactivation of the enzyme [cal/mol]
    #   T_halfH - Temperature at which 50% of the enzyme is inactivated from high temperature [K]
    devParams_dict = {"rho_25C": [0.01066,  0.00873,   0.01610],
                      "delH_A":  [10798.18, 26018.51,  14931.94],
                      "delH_H":  [100000,   55990.75, -472379],
                      "T_halfH": [14184.5,  304.58,    148.45]}
    devParams_df = pd.DataFrame(devParams_dict, index=["Eggs", "Larvae", "Pupae"])

    # Look up development rate parameters from previously defined dataframe
    rho_25C = devParams_df.loc[lifeStage, "rho_25C"]
    delH_A  = devParams_df.loc[lifeStage, "delH_A"]
    delH_H  = devParams_df.loc[lifeStage, "delH_H"]
    T_halfH = devParams_df.loc[lifeStage, "T_halfH"]

    # Calculate hourly development rate
    r = ( rho_25C * (waterTemp_hrly/298.) * np.exp( (delH_A/R_const)*((1/298)-(1/waterTemp_hrly)) ) ) / ( 1 + np.exp( (delH_H/R_const)*((1/T_halfH)-(1/waterTemp_hrly)) ) )

    return r




## calc_survProb_temp_dly()
#  Calculates daily survival probability for a mosquito immature given water temperature and life stage.
#  Water temperature thresholds for survival are from Sec. S2.7 and Table S5 of Magori et al. 2009 (Skeeter Buster paper).
#  Args:
#    waterTemp_min - daily min water temperature [deg C]
#    waterTemp_max - daily max water temperature [deg C]
#    lifeStage     - life stage of immatures (options: "Eggs", "Larvae", "Pupae")
#  Output:
#    Pr_surv_temp  - daily survival probability given temperature
def calc_survProb_temp_dly(waterTemp_min, waterTemp_max, lifeStage):

    # Initialize output arrays
    Pr_surv_Tmin = np.full_like(waterTemp_min, np.nan)
    Pr_surv_Tmax = np.full_like(waterTemp_max, np.nan)

    # Define dataframe of water temperature thresholds for survival (in deg C).
    # The four thresholds (T0, T1, T2, T3) define daily survival probabilities (Pr_surv_Tmin, Pr_surv_Tmax) 
    # based on daily min/max water temps (waterTemp_min, waterTemp_max).
    # See Sec. S2.7 of Magori et al. 2009 for more details.
    tempThresholds_dict = {"T0": [-14, 5, 5],
                           "T1": [-6, 10, 10],
                           "T2": [30, 39, 39],
                           "T3": [47, 44, 44]}
    tempThresholds_df = pd.DataFrame(tempThresholds_dict, index=["Eggs", "Larvae", "Pupae"])

    # Look up temperature thresholds from previously defined dataframe
    T0 = tempThresholds_df.loc[lifeStage, "T0"]
    T1 = tempThresholds_df.loc[lifeStage, "T1"]
    T2 = tempThresholds_df.loc[lifeStage, "T2"]
    T3 = tempThresholds_df.loc[lifeStage, "T3"]

    # Determine survival probability due to min temperature
    Tmin_regime_low = waterTemp_min <= T0
    Tmin_regime_mid = (waterTemp_min > T0) & (waterTemp_min < T1)
    Tmin_regime_hi  = waterTemp_min >= T1
    Pr_surv_Tmin[Tmin_regime_low] = 0.05
    Pr_surv_Tmin[Tmin_regime_mid] = 0.05 + 0.95*((waterTemp_min[Tmin_regime_mid]-T0)/(T1-T0))
    Pr_surv_Tmin[Tmin_regime_hi]  = 1.00

    # Determine survival probability due to max temperature
    Tmax_regime_low = waterTemp_max <= T2
    Tmax_regime_mid = (waterTemp_max > T2) & (waterTemp_max < T3)
    Tmax_regime_hi  = waterTemp_max >= T3
    Pr_surv_Tmax[Tmax_regime_low] = 1.00
    Pr_surv_Tmax[Tmax_regime_mid] = 1.00 - 0.95*((waterTemp_max[Tmax_regime_mid]-T2)/(T3-T2))
    Pr_surv_Tmax[Tmax_regime_hi]  = 0.05

    # Calculate overall temperature-dependent survival probability
    Pr_surv_temp = Pr_surv_Tmin * Pr_surv_Tmax

    return Pr_surv_temp




## calc_survProb_desicc_dly()
#  Calculates daily survival probability for a mosquito immature given desiccation and life stage.
#  For mosquito immatures this ONLY applies if the container dries out.
#  Note: VPD (used here) is called SD (saturation deficit) in the Skeeter Buster documentation.
#  Args:
#    waterHeight      - daily mean water height [mm]
#    vpd              - daily mean vapor pressure deficit [kPa]
#    sunExp           - sun exposure fraction (constant for a given container)
#    lifeStage        - life stage of immatures (options: "Eggs", "Larvae", "Pupae")
#  Output:
#    Pr_surv_desicc   - daily survival probability given desiccation
def calc_survProb_desicc_dly(waterHeight, vpd, sunExp, lifeStage):

    # Initialize output array
    Pr_surv_desicc = np.full_like(vpd, np.nan)

    # Convert VPD from kPa (as it is in WHATCH'EM output file)
    # to mbar (as it is in Skeeter Buster equations).
    vpd = vpd*10.

    # Get desiccation-dependent survival probability (only applies if container dried out)
    if lifeStage == "Eggs": # survival impact of eggs depends on sunExp and VPD
        if sunExp > 0.85:
            Pr_surv_desicc[:] = 0.95
        else: # if sunExp <= 0.85
            vpd_regime_low = vpd <= 10
            vpd_regime_mid = (vpd > 10) & (vpd < 30)
            vpd_regime_hi  = vpd >= 30
            Pr_surv_desicc[vpd_regime_low] = 0.99
            Pr_surv_desicc[vpd_regime_mid] = 0.99 - 0.04*((vpd[vpd_regime_mid]-10)/(30-10))
            Pr_surv_desicc[vpd_regime_hi]  = 0.95
    elif lifeStage == "Larvae": # survival impact on larvae is a flat value
        Pr_surv_desicc[:] = 0.05
    elif lifeStage == "Pupae":  # no survival impact on pupae
        Pr_surv_desicc[:] = 1.00

    # For the days on which the container isn't dried out, overwrite the survival factor with 1.00
    minWaterHeight = 0.01   # [mm]; this seems like a reasonable threshold for the daily mean water height, as
                            # it's the smallest hourly water height that can be recorded in the WHATCH'EM output
    Pr_surv_desicc[waterHeight >= minWaterHeight] = 1.00

    return Pr_surv_desicc




## calc_survProb_dly()
#  Calculates overall daily survival probability for a mosquito immature given all environmental factors and life stage.
#  Args:
#    waterTemp_min - daily min water temperature [deg C]
#    waterTemp_max - daily max water temperature [deg C]
#    waterHeight   - daily water height [mm]
#    vpd           - daily vapor pressure deficit [kPa]
#    sunExp        - sun exposure fraction (constant for a given container)
#    lifeStage     - life stage of immatures (options: "Eggs", "Larvae", "Pupae")
#  Output:
#    Pr_surv       - daily survival probability
def calc_survProb_dly(waterTemp_min, waterTemp_max, waterHeight, vpd, sunExp, lifeStage):

    Pr_surv_nominal = 0.99   # default daily Pr(survival) given no other influences
    Pr_surv_temp    = calc_survProb_temp_dly(waterTemp_min, waterTemp_max, lifeStage) # daily Pr(survival) due to temperature
    Pr_surv_desicc  = calc_survProb_desicc_dly(waterHeight, vpd, sunExp, lifeStage)   # daily Pr(survival) due to desiccation

    Pr_surv = Pr_surv_nominal * Pr_surv_temp * Pr_surv_desicc
    return Pr_surv




## calc_monRecord_dev_surv()
#  Calculates daily development rates and survival factors for mosquito immatures given
#  a (month-long) timespan of WHATCH'EM data.
#  Args:
#    dataSrc     - source of climate data ("MERRA_IMERG" or "GEOS")
#    loc         - geographic location/city name ("Negombo", "Jaffna", or "NuwaraEliya")
#    year        - year of WHATCH'EM data
#    month       - month of WHATCH'EM data (0 to 12)
#    intv        - "intv00", "intvF0", "intv0E", or "intvFE" (0 is no intv, F is fill, E is empty)
#    lead_time   - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#    ens_num     - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
#  Outputs:
#    df_dev_hrly     - dataframe of hourly dev rates (also includes select WHATCH'EM inputs/outputs)
#    df_dev_surv_dly - dataframe of daily dev rates and survival factors (also includes select WHATCH'EM inputs/outputs)
def calc_monRecord_dev_surv(dataSrc, loc, year, month, intv, lead_time=None, ens_num=None):

    # Validate GEOS-specific args, then set string used in data source-specific directories/filenames
    if dataSrc == 'GEOS':
        # Error if GEOS-specific args not provided when they're expected
        if lead_time is None or ens_num is None:
            raise ValueError('For dataSrc="GEOS", both lead_time and ens must be provided.')
    else:
        # Warn if GEOS-specific args were provided but won't be used
        if lead_time is not None or ens_num is not None:
            print('Warning: provided lead_time and ens are ignored when dataSrc is not "GEOS".')

    filepath_in, _, _ = WHATCHEM_specs_to_fpath(dataSrc, loc, year, month, intv, lead_time, ens_num)
    df_in = read_WHATCHEM_output(filepath_in) # note: converts the -9999.0 missing values into nan values

    # Prepare the inputs for development rate/survival factor calculations and the variables we want to carry through
    # from the WHATCHEM output (the original climate inputs and container outputs of WHATCH'EM, which we include in the
    # same dfs as the dev rate/survival factor variables so that we can easily compare all these variables to one another.
    waterHeight_hrly = df_in["WH (mm)"]     # hourly water heights [mm]
    waterTemp_hrly   = df_in["TW (C)"]      # hourly water temperatures [deg C]
    prcp_hrly        = df_in["PRCP (mm)"]   # hourly precipitation [mm]
    airTemp_hrly     = df_in["TA (C)"]      # hourly air temperature [deg C]
    vpd_hrly         = df_in["VPD (kPa)"]   # hourly vapor pressure deficit [kPa]
    sunExp           = 1 - float(df_in.attrs["metadata"]["Shade"]) # sun exposure fraction based on shading fraction

    waterHeight = waterHeight_hrly.resample("D").mean()           # daily mean water height [mm]
    waterHeight_del_hrly = -1*waterHeight_hrly.diff(periods=-1)   # hourly change in water height [mm]
                                                                  # (shows how height WILL change in next time step, as this
                                                                  # aligns with the precip value at the current time step)
    waterHeight_del_hrly = waterHeight_del_hrly.rename("WH_del (mm)")
    waterHeight_del = waterHeight_del_hrly.resample("D").sum()    # daily change in water height [mm]
                                                                  # (shows how height WILL change in next time step, as this
                                                                  # aligns with the precip value at the current time step)
    prcp = prcp_hrly.resample("D").sum()                          # daily total precipitation [mm]
    vpd  = df_in["VPD (kPa)"].resample("D").mean()                # daily mean vapor pressure deficit [kPa]
    waterTemp, waterTemp_min, waterTemp_max = convert_hrly2dlyMeanMinMax(waterTemp_hrly, "TW_min (C)", "TW_max (C)")  # daily mean/min/max water temperature [deg C] (nan if any hourly value is nan)
    airTemp, airTemp_min, airTemp_max       = convert_hrly2dlyMeanMinMax(airTemp_hrly, "TA_min (C)", "TA_max (C)")    # daily mean/min/max air temperature [deg C] (nan if any hourly value is nan, but that shouldn't be the case for this var)

    # Create water temp variables where missing values are estimated and filled in.
    # This estimation is based on a linear regression of the existing water temp data
    # against air temp.
    waterTemp_nanFill = fillMissingVals_linReg(waterTemp, airTemp)
    waterTemp_nanFill = waterTemp_nanFill.rename("TW_nanFill (C)")
    waterTemp_min_nanFill = fillMissingVals_linReg(waterTemp_min, airTemp_min)
    waterTemp_min_nanFill = waterTemp_min_nanFill.rename("TW_min_nanFill (C)")
    waterTemp_max_nanFill = fillMissingVals_linReg(waterTemp_max, airTemp_max)
    waterTemp_max_nanFill = waterTemp_max_nanFill.rename("TW_max_nanFill (C)")
    # The process for creating the hourly water temp variable is more involved.
    # 1. Create an hourly variable waterTemp_hrly_ALLnanFill where ALL missing values are filled with the daily means from waterTemp_nanFill.
    # 2. Create a mask mask_lowWH_and_noTW for the specific nan values we actually want to replace (only those arising due to low water height).
    # 3. Use the mask to create waterTemp_hrly_nanFill, where missing hourly values corresponding to the mask are filled with values from waterTemp_hrly_ALLnanFill.
    waterTemp_hrly_ALLnanFill = waterTemp_nanFill.reindex(pd.date_range(waterTemp_nanFill.index[0], waterTemp_nanFill.index[-1] + pd.Timedelta(hours=23), freq='h'),
                                                          method='ffill')
    mask_lowWH_noTW = (waterHeight_hrly <= 15.00) & np.isnan(waterTemp_hrly)
    waterTemp_hrly_nanFill = waterTemp_hrly.copy(deep=True)
    waterTemp_hrly_nanFill[mask_lowWH_noTW] = waterTemp_hrly_ALLnanFill[mask_lowWH_noTW.index[mask_lowWH_noTW]]
    waterTemp_hrly_nanFill = waterTemp_hrly_nanFill.rename("TW_nanFill (C)")

    # Put all these variables neatly together into dfs.
    df_WHATCHEM_hrly = pd.concat([waterHeight_hrly, waterHeight_del_hrly, waterTemp_hrly, waterTemp_hrly_nanFill,
                                  prcp_hrly, airTemp_hrly, vpd_hrly], axis=1)
    df_WHATCHEM_dly  = pd.concat([waterHeight, waterHeight_del, waterTemp, waterTemp_min, waterTemp_max, 
                                  waterTemp_nanFill, waterTemp_min_nanFill, waterTemp_max_nanFill,
                                  prcp, airTemp, airTemp_min, airTemp_max, vpd], axis=1)

    lifeStages = ["Eggs", "Larvae", "Pupae"]


    # Create dataframe to store hourly development rate
    devCols = ["dev_E", "dev_L", "dev_P"]
    df_dev_hrly = pd.DataFrame(columns = ["Date"] + devCols)
    df_dev_hrly["Date"] = df_in.index
    df_dev_hrly.set_index("Date", inplace = True)

    # Compute hourly development rates and store these values
    devCols_dict = dict(zip(lifeStages, devCols))
    for lifeStage, colName in devCols_dict.items():
        df_dev_hrly[colName] = calc_devRate(waterTemp_hrly_nanFill, lifeStage)
    
    # Compute daily development rates from hourly values, ensuring NaN if any NaN value exists in a day
    df_dev_dly = df_dev_hrly.resample("D").apply(lambda x: np.nan if any(np.isnan(x)) else np.nansum(x))

    
    # Create dataframe to store daily survival factor
    survCols = ["surv_E", "surv_L", "surv_P"]
    df_surv_dly = pd.DataFrame(columns = ["Date"] + survCols)
    df_surv_dly["Date"] = pd.to_datetime(np.unique(df_in.index.date))
    df_surv_dly.set_index("Date", inplace = True)

    survCols_temp  = ["surv_temp_E", "surv_temp_L", "surv_temp_P"]
    df_surv_temp_dly = df_surv_dly.copy(deep=True)
    df_surv_temp_dly.rename(columns=dict(zip(df_surv_dly.columns, survCols_temp)), inplace=True)
    
    survCols_desicc = ["surv_desicc_E", "surv_desicc_L", "surv_desicc_P"]
    df_surv_desicc_dly = df_surv_dly.copy(deep=True)
    df_surv_desicc_dly.rename(columns=dict(zip(df_surv_dly.columns, survCols_desicc)), inplace=True)

    survCols_dict = dict(zip(lifeStages, survCols))
    for lifeStage, colName in survCols_dict.items():
        df_surv_dly[colName] = calc_survProb_dly(waterTemp_min_nanFill, waterTemp_max_nanFill, waterHeight, vpd, sunExp, lifeStage)
    survCols_temp_dict = dict(zip(lifeStages, survCols_temp))
    for lifeStage, colName in survCols_temp_dict.items():
        df_surv_temp_dly[colName] = calc_survProb_temp_dly(waterTemp_min_nanFill, waterTemp_max_nanFill, lifeStage)
    survCols_desicc_dict = dict(zip(lifeStages, survCols_desicc))
    for lifeStage, colName in survCols_desicc_dict.items():
        df_surv_desicc_dly[colName] = calc_survProb_desicc_dly(waterHeight, vpd, sunExp, lifeStage)

    # Concatenate all dfs along columns (axis=1)
    df_dev_hrly = pd.concat([df_dev_hrly, df_WHATCHEM_hrly], axis=1)
    df_dev_surv_dly = pd.concat([df_dev_dly, df_surv_dly, df_surv_temp_dly, df_surv_desicc_dly, df_WHATCHEM_dly], axis=1)

    # Add metadata from WHATCH'EM output to these dfs
    df_dev_hrly.attrs["metadata"]     = df_in.attrs["metadata"]
    df_dev_surv_dly.attrs["metadata"] = df_in.attrs["metadata"]

    return df_dev_hrly, df_dev_surv_dly




## save_monRecord_dev_surv()
#  Runs calc_monRecord_dev_surv() and saves the resulting dataframes as text files.
#  Args:
#    dataSrc     - source of climate data ("MERRA_IMERG" or "GEOS_S2S")
#    loc         - geographic location/city name ("Negombo", "Jaffna", or "NuwaraEliya")
#    year        - year of WHATCH'EM data
#    month       - month of WHATCH'EM data (0 to 12)
#    intv        - intv00, intvF0, intv0E, or intvFE (0 is no intv, F is fill, E is empty)
#    lead_time   - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#    ens_num     - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
def save_monRecord_dev_surv(dataSrc, loc, year, month, intv, lead_time=None, ens_num=None):

    # Validate GEOS-specific args
    if dataSrc == 'GEOS':
        # Error if GEOS-specific args not provided when they're expected
        if lead_time is None or ens_num is None:
            raise ValueError('For dataSrc="GEOS", both lead_time and ens must be provided.')
    else:
        # Warn if GEOS-specific args were provided but won't be used
        if lead_time is not None or ens_num is not None:
            print('Warning: provided lead_time and ens are ignored when dataSrc is not "GEOS".')

    df_dev_hrly, df_dev_surv_dly = calc_monRecord_dev_surv(dataSrc, loc, year, month, intv, lead_time, ens_num)
    _, subDir, fileStr_allSpecs = WHATCHEM_specs_to_fpath(dataSrc, loc, year, month, intv, lead_time, ens_num)
    
    dir_out = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_BASE, dataSrc, subDir)
    os.makedirs(dir_out, exist_ok=True)

    # Save hourly df
    with open(os.path.join(dir_out, f'df_dev_hrly_{fileStr_allSpecs}.txt'), 'w') as f:
        # Convert metadata to a string format suitable for a header
        metadata_str = ''
        for key, value in df_dev_hrly.attrs["metadata"].items():
            metadata_str += f"# {key}: {value}\n" # Each metadata item is on its own line that starts with a '#'
        # Save df with metadata
        f.write(metadata_str)
        df_dev_hrly.to_csv(f, sep="\t", na_rep='NaN', index=True)

    # Save daily  df
    with open(os.path.join(dir_out, f'df_dev_surv_dly_{fileStr_allSpecs}.txt'), 'w') as f:
        # Convert metadata to a string format suitable for a header
        metadata_str = ''
        for key, value in df_dev_surv_dly.attrs["metadata"].items():
            metadata_str += f"# {key}: {value}\n" # Each metadata item is on its own line that starts with a '#'
        # Save df with metadata
        f.write(metadata_str)
        df_dev_surv_dly.to_csv(f, sep="\t", na_rep='NaN', index=True)




## calc_monRecord_vectorPop()
#  Calculates daily population of mosquito immatures and adults given an initial seeding of 
#  immatures and a (month-long) timespan of WHATCH'EM data.
#  Stores results in txt file.
#  Args:
#    dataSrc     - source of climate data ("MERRA_IMERG" or "GEOS_S2S")
#    loc         - geographic location/city name ("Negombo", "Jaffna", or "NuwaraEliya")
#    year        - year of WHATCH'EM data
#    month       - month of WHATCH'EM data (0 to 12)
#    intv        - intv00, intvF0, intv0E, or intvFE (0 is no intv, F is fill, E is empty)
#    initDay     - day of month on which to introduce immatures (e.g., Jan 2 would be 2)
#    initEggs    - initial population of eggs
#    initLarvae  - initial population of larvae
#    initPupae   - initial population of pupae
#    lead_time   - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#    ens_num     - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
#  Outputs:
#    df_vectorPop_dly - dataframe of daily vector populations (also includes select WHATCH'EM inputs/outputs and dev rates/survival factors)
def calc_monRecord_vectorPop(dataSrc, loc, year, month, intv, initDay, initEggs, initLarvae, initPupae, lead_time=None, ens_num=None):

    # Validate GEOS-specific args, then set string used in data source-specific directories/filenames
    if dataSrc == 'GEOS':
        # Error if GEOS-specific args not provided when they're expected
        if lead_time is None or ens_num is None:
            raise ValueError('For dataSrc="GEOS", both lead_time and ens must be provided.')
    else:
        # Warn if GEOS-specific args were provided but won't be used
        if lead_time is not None or ens_num is not None:
            print('Warning: provided lead_time and ens are ignored when dataSrc is not "GEOS".')

    _, df_dev_surv_dly = calc_monRecord_dev_surv(dataSrc, loc, year, month, intv, lead_time, ens_num)
    containerHeight = 1000*float(df_dev_surv_dly.attrs["metadata"]["Height"])  # container height [mm]

    df_vectorPopCols = pd.DataFrame(0, index=df_dev_surv_dly.index, 
                                    columns=['pop_E', 'pop_L', 'pop_P', 'pop_A', 
                                             'pop_L(E)', 'pop_L(L)', 'pop_P(E)', 'pop_P(L)', 'pop_P(P)', 
                                             'pop_A(E)', 'pop_A(L)', 'pop_A(P)'])
    df_vectorPop_dly = pd.concat([df_vectorPopCols, df_dev_surv_dly], axis=1)

    # Define dataframe for tracking each mosquito individual
    vals_init_stage    = ['egg'] * initEggs + ['larva'] * initLarvae + ['pupa'] * initPupae
    vals_current_stage = vals_init_stage
    vals_is_alive      = [True] * (initEggs + initLarvae + initPupae)
    vals_cumul_dev     = [0.0] * (initEggs + initLarvae + initPupae)
    # Randomize initial egg heights within a band of heights above the water.
    # Allowed heights are in 2 mm increments above the water, as in Skeeter Buster.
    initDate = df_dev_surv_dly.index[initDay-1]
    initWaterHeight = df_dev_surv_dly.loc[initDate, 'WH (mm)']                         # initial water height in container [mm]
    eggHeight_max   = 40                                                               # max height above water for egg-laying [mm]
    eggBandSize     = 2                                                                # size of each height band for egg-laying [mm]
    eggHeightRange  = int(min(round(containerHeight-initWaterHeight), eggHeight_max))  # range of heights above water for egg-laying (can't go above top of container!) [mm]
    vals_egg_height = [initWaterHeight + random.choice(range(0, eggHeightRange, eggBandSize)) for _ in range(initEggs)] + ([np.nan] * (initLarvae + initPupae))
    vals_is_mature_egg = [False] * (initEggs + initLarvae + initPupae)

    df_vectorStates  = pd.DataFrame({'init_stage': vals_init_stage,
                                     'current_stage': vals_current_stage,
                                     'is_alive': vals_is_alive,
                                     'cumul_dev': vals_cumul_dev,
                                     'egg_height (mm)': vals_egg_height,
                                     'is_mature_egg': vals_is_mature_egg})

    # Iterate through the rows (i.e. dates).
    for idx_date, row_vectorPop in df_vectorPop_dly.iterrows():

        # first, check whether the day we're looking at is after/at the day on which the immatures are added to the container
        if idx_date >= initDate:
            # update the individual mosquitoes and then calculate populations based on those results
            # the population values for a given day are the values at the END of said day (so the first day's values are not necessarily the initial values)
            for idx_vectorNum, _ in df_vectorStates.iterrows():

                if df_vectorStates.at[idx_vectorNum, 'is_alive']: # if individual is alive

                    if df_vectorStates.at[idx_vectorNum, 'current_stage'] == 'egg': # if individual is an egg

                        # check whether individual dies
                        roll_4_dlySurvival = random.random() # generate random number in interval [0,1)

                        if roll_4_dlySurvival <= (1-row_vectorPop['surv_E']): # if individual dies
                            df_vectorStates.at[idx_vectorNum, 'is_alive'] = False # update alive status

                        else: # if individual survives

                            is_warm_water = row_vectorPop['TW (C)'] >= 22.0 # whether water is warm enough [in deg C] for hatching
                            is_submerged = row_vectorPop['WH (mm)'] > df_vectorStates.at[idx_vectorNum, 'egg_height (mm)'] # whether egg is submerged in water

                            if df_vectorStates.at[idx_vectorNum, 'is_mature_egg']: # if egg is mature

                                # Check for hatching based on mean daily water temp and submergence
                                if is_warm_water and is_submerged:
                                    roll_4_hatch = random.random() # generate random number in interval [0,1)
                                    if roll_4_hatch <= 0.596: # if egg hatches
                                        df_vectorStates.at[idx_vectorNum, 'current_stage'] = 'larva' # update vector stage
                                        df_vectorStates.at[idx_vectorNum, 'cumul_dev'] = 0.0         # reset development
                                        df_vectorStates.at[idx_vectorNum, 'is_mature_egg'] = False   # reset maturity tracking

                            else: # if egg not mature
                                
                                if df_vectorStates.at[idx_vectorNum, 'cumul_dev'] > 0.95: # if developed
                                    if is_warm_water:
                                        if is_submerged:
                                            # do egg hatching
                                            df_vectorStates.at[idx_vectorNum, 'current_stage'] = 'larva' # update vector stage
                                            df_vectorStates.at[idx_vectorNum, 'cumul_dev'] = 0.0         # reset development
                                            df_vectorStates.at[idx_vectorNum, 'is_mature_egg'] = False   # reset maturity tracking
                                        else: # if not submerged
                                            roll_4_hatch = random.random() # generate random number in interval [0,1)
                                            if roll_4_hatch <= 0.197: # if egg hatches
                                                df_vectorStates.at[idx_vectorNum, 'current_stage'] = 'larva' # update vector stage
                                                df_vectorStates.at[idx_vectorNum, 'cumul_dev'] = 0.0         # reset development
                                                df_vectorStates.at[idx_vectorNum, 'is_mature_egg'] = False   # reset maturity tracking
                                            else: # if egg doesn't hatch
                                                df_vectorStates.at[idx_vectorNum, 'is_mature_egg'] = True    # egg becomes mature
                                            
                                    else: # if water not warm
                                        df_vectorStates.at[idx_vectorNum, 'is_mature_egg'] = True # egg becomes mature (but no hatching)

                                else: # if not developed
                                    df_vectorStates.at[idx_vectorNum, 'cumul_dev'] += df_vectorPop_dly.at[idx_date, 'dev_E'] # update cumulative development

                    elif df_vectorStates.at[idx_vectorNum, 'current_stage'] == 'larva': # if individual is a larva

                        # check whether individual dies
                        roll_4_dlySurvival = random.random() # generate random number in interval [0,1)
                    
                        if roll_4_dlySurvival <= (1-row_vectorPop['surv_L']): # if individual dies
                            df_vectorStates.at[idx_vectorNum, 'is_alive'] = False # update alive status

                        else: # if individual survives
                            if df_vectorStates.at[idx_vectorNum, 'cumul_dev'] > 0.95: # if developed

                                # attempt to pupate
                                roll_4_pupation = random.random() # generate random number in interval [0,1)
                                if roll_4_pupation <= 0.95: # if pupation successful
                                    # do pupation
                                    df_vectorStates.at[idx_vectorNum, 'current_stage'] = 'pupa' # update vector stage
                                    df_vectorStates.at[idx_vectorNum, 'cumul_dev'] = 0.0        # reset development
                                else: # if pupation unsuccessful (death)
                                    df_vectorStates.at[idx_vectorNum, 'is_alive'] = False       # update alive status

                            else: # if not developed
                                df_vectorStates.at[idx_vectorNum, 'cumul_dev'] += df_vectorPop_dly.at[idx_date, 'dev_L'] # update cumulative development

                    elif df_vectorStates.at[idx_vectorNum, 'current_stage'] == 'pupa': # if individual is a pupa

                        # check whether individual dies
                        roll_4_dlySurvival = random.random() # generate random number in interval [0,1)
                    
                        if roll_4_dlySurvival <= (1-row_vectorPop['surv_P']): # if individual dies
                            df_vectorStates.at[idx_vectorNum, 'is_alive'] = False # update alive status

                        else: # if individual survives
                            if df_vectorStates.at[idx_vectorNum, 'cumul_dev'] > 0.95: # if developed

                                # attempt to eclose
                                roll_4_eclosion = random.random() # generate random number in interval [0,1)
                                if roll_4_eclosion <= 0.83: # if eclosion successful
                                    # do eclosion
                                    df_vectorStates.at[idx_vectorNum, 'current_stage'] = 'adult' # update vector stage
                                    df_vectorStates.at[idx_vectorNum, 'cumul_dev'] = 0.0         # reset development
                                else: # if eclosion unsuccessful (death)
                                    df_vectorStates.at[idx_vectorNum, 'is_alive'] = False        # update alive status
                                    
                            else: # if not developed
                                df_vectorStates.at[idx_vectorNum, 'cumul_dev'] += df_vectorPop_dly.at[idx_date, 'dev_P'] # update cumulative development
                    
                    else: # adult
                        continue

            # record total alive population in each stage (e.g., L is all larvae)
            df_vectorPop_dly.at[idx_date, 'pop_E'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & (df_vectorStates['current_stage'] == 'egg')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_L'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & (df_vectorStates['current_stage'] == 'larva')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_P'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & (df_vectorStates['current_stage'] == 'pupa')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_A'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & (df_vectorStates['current_stage'] == 'adult')].shape[0]

            # record alive population in each stage that originated from a specific initial stage (e.g., L(E) is larvae from egg)
            df_vectorPop_dly.at[idx_date, 'pop_L(E)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'larva') & (df_vectorStates['init_stage'] == 'egg')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_L(L)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'larva') & (df_vectorStates['init_stage'] == 'larva')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_P(E)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'pupa') & (df_vectorStates['init_stage'] == 'egg')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_P(L)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'pupa') & (df_vectorStates['init_stage'] == 'larva')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_P(P)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'pupa') & (df_vectorStates['init_stage'] == 'pupa')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_A(E)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'adult') & (df_vectorStates['init_stage'] == 'egg')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_A(L)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'adult') & (df_vectorStates['init_stage'] == 'larva')].shape[0]
            df_vectorPop_dly.at[idx_date, 'pop_A(P)'] = df_vectorStates[(df_vectorStates['is_alive'] == True) & 
                                                                        (df_vectorStates['current_stage'] == 'adult') & (df_vectorStates['init_stage'] == 'pupa')].shape[0]
    
    # Add metadata to df
    df_vectorPop_dly.attrs["metadata"] = df_dev_surv_dly.attrs["metadata"]

    return df_vectorPop_dly




## save_monRecord_vectorPop()
#  Runs calc_monRecord_vectorPop() and saves the resulting dataframe as a text file.
#  Args:
#    dataSrc     - source of climate data ("MERRA_IMERG" or "GEOS_S2S")
#    loc         - geographic location/city name ("Negombo", "Jaffna", or "NuwaraEliya")
#    year        - year of WHATCH'EM data
#    month       - month of WHATCH'EM data (0 to 12)
#    intv        - intv00, intvF0, intv0E, or intvFE (0 is no intv, F is fill, E is empty)
#    initDay     - day of month on which to introduce immatures (e.g., Jan 2 would be 2)
#    initEggs    - initial population of eggs
#    initLarvae  - initial population of larvae
#    initPupae   - initial population of pupae
#    lead_time   - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#    ens_num     - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
def save_monRecord_vectorPop(dataSrc, loc, year, month, intv, initDay, initEggs, initLarvae, initPupae, lead_time=None, ens_num=None):

    # Validate GEOS-specific args
    if dataSrc == 'GEOS':
        # Error if GEOS-specific args not provided when they're expected
        if lead_time is None or ens_num is None:
            raise ValueError('For dataSrc="GEOS", both lead_time and ens must be provided.')
    else:
        # Warn if GEOS-specific args were provided but won't be used
        if lead_time is not None or ens_num is not None:
            print('Warning: provided lead_time and ens are ignored when dataSrc is not "GEOS".')

    df_vectorPop_dly = calc_monRecord_vectorPop(dataSrc, loc, year, month, intv, initDay, initEggs, initLarvae, initPupae, lead_time, ens_num)
    _, subDir, fileStr_allSpecs = WHATCHEM_specs_to_fpath(dataSrc, loc, year, month, intv, lead_time, ens_num)
    
    # create string indicating init populations of immatures and the day they're introduced to container
    # (format = "dayW_XeYlZp", where W is day of introduction and X, Y, Z are # of eggs, larvae, pupae)
    popStr = f'day{initDay}_{initEggs}e{initLarvae}l{initPupae}p' 

    dir_out = os.path.join(DIR_WHATCHEM_OUTPUT_DATA_BASE, dataSrc, subDir)
    os.makedirs(dir_out, exist_ok=True)

    # Save df
    with open(os.path.join(dir_out, f'df_vectorPop_dly_{fileStr_allSpecs}-{popStr}.txt'), 'w') as f:
        # Convert metadata to a string format suitable for a header
        metadata_str = ''
        for key, value in df_vectorPop_dly.attrs["metadata"].items():
            metadata_str += f"# {key}: {value}\n" # Each metadata item is on its own line that starts with a '#'
        # Save df with metadata
        f.write(metadata_str)
        df_vectorPop_dly.to_csv(f, sep="\t", na_rep='NaN', index=True)










# Scripts

# create text files of each month of daily data on pop sizes, dev rates, survival factors, and 
# associated variables (TW_min, TW_max, TW, WH, TA, VPD)
# This code uses the temp-based survival curves from Focks et al. 1993
if True:
    dataSrcs   = ['GEOS'] #["MERRA_IMERG", "GEOS"]
    lead_times = [0]           # only relevant for GEOS
    ens_nums   = range(1, 5+1) # only relevant for GEOS
    locs       = ["Negombo", "Jaffna", "NuwaraEliya"]
    years      = range(2001, 2020+1)
    months     = range(1, 12+1)
    intv       = "intv00"
    initDays   = range(1, 7+1)
    initEggs   = 1000
    initLarvae = 1000
    initPupae  = 1000
    for dataSrc in dataSrcs:
        print(f'{dataSrc}')
        for loc in locs:
            print(f'{loc}')
            for year in years:
                print(f' {year}')
                for month in months:
                    print(f'  {month}')

                    # Next, run the code to create the met input files
                    if dataSrc == 'GEOS':
                        for lead_time in lead_times:
                            print(f'   Lead time: {lead_time}')
                            for ens_num in ens_nums:
                                print(f'    ens{ens_num}')
                                print(f'    Running dev/surv calculations...')
                                save_monRecord_dev_surv(dataSrc, loc, year, month, intv, lead_time, ens_num)
                                print(f'    Running vector population calculations...')
                                for initDay in initDays:
                                    print(f'     initDay: {initDay}')
                                    save_monRecord_vectorPop(dataSrc, loc, year, month, intv, initDay, initEggs, initLarvae, initPupae, lead_time, ens_num)
                    else: # if MERRA_IMERG
                        print(f'    Running dev/surv calculations...')
                        save_monRecord_dev_surv(dataSrc, loc, year, month, intv)
                        print(f'    Running vector population calculations...')
                        for initDay in initDays:
                            print(f'     initDay: {initDay}')
                            save_monRecord_vectorPop(dataSrc, loc, year, month, intv, initDay, initEggs, initLarvae, initPupae)





## TEST
#
#  Test variability among runs for each year (with a fixed immature introduction day of 1 and a fixed month) and show as boxplots
if False:
    numRuns =  10
    initPops = 1000
    initDay  = 1
    year_first = 2001
    year_last  = 2020
    cities = ["Negombo", "NuwaraEliya", "Jaffna"]
    months = [7, 6, 12] # corresponding to each city

    if False:
        listOfLists = []
        for runNum in range(numRuns):
            print('runNum: ' + str(runNum))
            numAdults = []
            for year in range(year_first, year_last+1):
                df_vectorPop_dly = calc_monRecord_vectorPop("MERRA_IMERG", cities[0], year, months[0], "intv00", initDay, initPops, initPops, initPops)
                numAdults.append(df_vectorPop_dly.iloc[-1]["pop_A"])
            listOfLists.append(numAdults)
        df_Negombo = pd.DataFrame(listOfLists, columns=[str(i) for i in range(year_first, year_last+1)])
        print(df_Negombo)
        
        listOfLists = []
        for runNum in range(numRuns):
            print('runNum: ' + str(runNum))
            numAdults = []
            for year in range(year_first, year_last+1):
                df_vectorPop_dly = calc_monRecord_vectorPop("MERRA_IMERG", cities[1], year, months[1], "intv00", initDay, initPops, initPops, initPops)
                numAdults.append(df_vectorPop_dly.iloc[-1]["pop_A"])
            listOfLists.append(numAdults)
        df_NuwaraEliya = pd.DataFrame(listOfLists, columns=[str(i) for i in range(year_first, year_last+1)])
        print(df_NuwaraEliya)

        listOfLists = []
        for runNum in range(numRuns):
            print('runNum: ' + str(runNum))
            numAdults = []
            for year in range(year_first, year_last+1):
                df_vectorPop_dly = calc_monRecord_vectorPop("MERRA_IMERG", cities[2], year, months[2], "intv00", initDay, initPops, initPops, initPops)
                numAdults.append(df_vectorPop_dly.iloc[-1]["pop_A"])
            listOfLists.append(numAdults)
        df_Jaffna = pd.DataFrame(listOfLists, columns=[str(i) for i in range(year_first, year_last+1)])
        print(df_Jaffna)

    data_Negombo     = [[1137.0, 752.0, 1560.0, 1515.0, 1435.0, 1637.0, 1328.0, 1551.0, 1551.0, 1380.0, 1526.0, 1390.0, 1522.0, 1566.0, 1535.0, 1520.0, 1566.0, 1481.0, 1148.0, 1552.0],
                        [1158.0, 815.0, 1558.0, 1526.0, 1440.0, 1646.0, 1409.0, 1576.0, 1544.0, 1387.0, 1484.0, 1400.0, 1515.0, 1534.0, 1497.0, 1548.0, 1507.0, 1461.0, 1160.0, 1564.0],
                        [1144.0, 813.0, 1573.0, 1581.0, 1419.0, 1628.0, 1355.0, 1544.0, 1552.0, 1351.0, 1506.0, 1443.0, 1522.0, 1539.0, 1560.0, 1553.0, 1548.0, 1456.0, 1161.0, 1511.0],
                        [1139.0, 818.0, 1561.0, 1548.0, 1450.0, 1607.0, 1382.0, 1566.0, 1568.0, 1314.0, 1500.0, 1432.0, 1540.0, 1548.0, 1548.0, 1533.0, 1548.0, 1482.0, 1171.0, 1516.0],
                        [1109.0, 807.0, 1571.0, 1553.0, 1439.0, 1591.0, 1387.0, 1536.0, 1518.0, 1339.0, 1514.0, 1442.0, 1496.0, 1541.0, 1544.0, 1504.0, 1529.0, 1425.0, 1137.0, 1569.0],
                        [1147.0, 814.0, 1548.0, 1573.0, 1466.0, 1597.0, 1410.0, 1537.0, 1508.0, 1358.0, 1473.0, 1416.0, 1527.0, 1546.0, 1522.0, 1552.0, 1518.0, 1478.0, 1168.0, 1561.0],
                        [1150.0, 837.0, 1575.0, 1571.0, 1439.0, 1558.0, 1394.0, 1564.0, 1520.0, 1339.0, 1497.0, 1441.0, 1543.0, 1528.0, 1557.0, 1525.0, 1540.0, 1466.0, 1164.0, 1539.0],
                        [1146.0, 830.0, 1531.0, 1594.0, 1432.0, 1586.0, 1387.0, 1561.0, 1556.0, 1376.0, 1514.0, 1390.0, 1555.0, 1524.0, 1558.0, 1571.0, 1515.0, 1468.0, 1168.0, 1541.0],
                        [1146.0, 826.0, 1584.0, 1513.0, 1458.0, 1573.0, 1408.0, 1538.0, 1603.0, 1304.0, 1489.0, 1464.0, 1556.0, 1550.0, 1548.0, 1550.0, 1542.0, 1446.0, 1168.0, 1582.0],
                        [1169.0, 846.0, 1510.0, 1552.0, 1471.0, 1620.0, 1385.0, 1554.0, 1587.0, 1373.0, 1486.0, 1435.0, 1530.0, 1544.0, 1544.0, 1555.0, 1508.0, 1488.0, 1161.0, 1566.0]
                       ]
    df_Negombo = pd.DataFrame(data_Negombo, columns=[str(i) for i in range(year_first, year_last+1)])

    data_NuwaraEliya = [[1385.0, 1396.0, 1444.0, 1433.0, 1450.0, 1438.0, 1451.0, 1431.0, 1456.0, 1438.0, 1160.0, 1410.0, 1398.0, 1435.0, 1415.0, 1462.0, 1411.0, 1410.0, 1466.0, 1448.0],
                        [1416.0, 1444.0, 1464.0, 1452.0, 1442.0, 1409.0, 1418.0, 1423.0, 1445.0, 1436.0, 1138.0, 1466.0, 1386.0, 1425.0, 1440.0, 1459.0, 1411.0, 1471.0, 1489.0, 1460.0],
                        [1419.0, 1396.0, 1480.0, 1424.0, 1438.0, 1443.0, 1450.0, 1411.0, 1408.0, 1433.0, 1142.0, 1448.0, 1419.0, 1449.0, 1421.0, 1445.0, 1426.0, 1444.0, 1435.0, 1429.0],
                        [1414.0, 1463.0, 1460.0, 1404.0, 1451.0, 1459.0, 1471.0, 1379.0, 1465.0, 1398.0, 1138.0, 1428.0, 1399.0, 1430.0, 1473.0, 1435.0, 1429.0, 1448.0, 1485.0, 1463.0],
                        [1415.0, 1424.0, 1469.0, 1440.0, 1471.0, 1395.0, 1424.0, 1413.0, 1456.0, 1421.0, 1167.0, 1472.0, 1416.0, 1438.0, 1426.0, 1457.0, 1446.0, 1442.0, 1444.0, 1432.0],
                        [1437.0, 1421.0, 1455.0, 1423.0, 1424.0, 1439.0, 1462.0, 1433.0, 1469.0, 1432.0, 1150.0, 1461.0, 1378.0, 1437.0, 1414.0, 1455.0, 1464.0, 1409.0, 1478.0, 1436.0],
                        [1405.0, 1432.0, 1460.0, 1415.0, 1447.0, 1469.0, 1448.0, 1422.0, 1420.0, 1452.0, 1143.0, 1448.0, 1393.0, 1437.0, 1438.0, 1407.0, 1439.0, 1423.0, 1445.0, 1438.0],
                        [1389.0, 1459.0, 1455.0, 1415.0, 1462.0, 1428.0, 1460.0, 1372.0, 1461.0, 1405.0, 1166.0, 1435.0, 1410.0, 1454.0, 1422.0, 1435.0, 1454.0, 1413.0, 1466.0, 1462.0],
                        [1443.0, 1448.0, 1402.0, 1430.0, 1441.0, 1406.0, 1439.0, 1427.0, 1474.0, 1434.0, 1139.0, 1431.0, 1385.0, 1414.0, 1429.0, 1427.0, 1433.0, 1420.0, 1467.0, 1411.0],
                        [1408.0, 1463.0, 1467.0, 1438.0, 1450.0, 1473.0, 1459.0, 1404.0, 1448.0, 1437.0, 1165.0, 1431.0, 1362.0, 1431.0, 1414.0, 1438.0, 1448.0, 1421.0, 1449.0, 1465.0]
                       ]
    df_NuwaraEliya = pd.DataFrame(data_NuwaraEliya, columns=[str(i) for i in range(year_first, year_last+1)])
    
    data_Jaffna      = [[1390.0, 1525.0, 1521.0, 1690.0, 1439.0, 1546.0, 1162.0, 1211.0, 1634.0, 1987.0, 1170.0, 1269.0, 1807.0, 1661.0, 1766.0, 1251.0, 1555.0, 1485.0, 1638.0, 1846.0],
                        [1411.0, 1505.0, 1511.0, 1655.0, 1463.0, 1520.0, 1241.0, 1249.0, 1640.0, 2001.0, 1179.0, 1296.0, 1826.0, 1644.0, 1788.0, 1256.0, 1513.0, 1470.0, 1627.0, 1909.0],
                        [1425.0, 1505.0, 1544.0, 1632.0, 1451.0, 1550.0, 1233.0, 1204.0, 1667.0, 2010.0, 1201.0, 1311.0, 1809.0, 1601.0, 1777.0, 1285.0, 1529.0, 1477.0, 1612.0, 1843.0],
                        [1423.0, 1492.0, 1508.0, 1676.0, 1481.0, 1530.0, 1168.0, 1236.0, 1637.0, 1956.0, 1191.0, 1306.0, 1841.0, 1632.0, 1762.0, 1275.0, 1533.0, 1530.0, 1641.0, 1888.0],
                        [1395.0, 1527.0, 1547.0, 1696.0, 1471.0, 1552.0, 1213.0, 1189.0, 1647.0, 1998.0, 1233.0, 1330.0, 1854.0, 1631.0, 1758.0, 1248.0, 1553.0, 1496.0, 1628.0, 1886.0],
                        [1409.0, 1520.0, 1533.0, 1706.0, 1479.0, 1518.0, 1232.0, 1221.0, 1655.0, 2011.0, 1179.0, 1336.0, 1883.0, 1661.0, 1801.0, 1268.0, 1561.0, 1513.0, 1649.0, 1915.0],
                        [1409.0, 1492.0, 1520.0, 1716.0, 1479.0, 1532.0, 1194.0, 1237.0, 1648.0, 2028.0, 1185.0, 1269.0, 1856.0, 1618.0, 1742.0, 1232.0, 1546.0, 1468.0, 1615.0, 1862.0],
                        [1355.0, 1481.0, 1538.0, 1694.0, 1469.0, 1560.0, 1206.0, 1234.0, 1653.0, 1957.0, 1208.0, 1303.0, 1807.0, 1640.0, 1823.0, 1260.0, 1542.0, 1427.0, 1661.0, 1880.0],
                        [1401.0, 1511.0, 1536.0, 1672.0, 1445.0, 1518.0, 1221.0, 1230.0, 1646.0, 1996.0, 1183.0, 1304.0, 1849.0, 1650.0, 1738.0, 1268.0, 1532.0, 1451.0, 1613.0, 1877.0],
                        [1416.0, 1518.0, 1531.0, 1656.0, 1471.0, 1508.0, 1203.0, 1210.0, 1641.0, 1986.0, 1215.0, 1317.0, 1847.0, 1647.0, 1764.0, 1222.0, 1554.0, 1467.0, 1664.0, 1881.0]
                       ]
    df_Jaffna = pd.DataFrame(data_Jaffna, columns=[str(i) for i in range(year_first, year_last+1)])

    # Plot box plots for each column
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    df_Negombo.boxplot(grid=False, boxprops=dict(color='k'), medianprops=dict(color='k'),
               whiskerprops=dict(color='k'), flierprops=dict(color='k'),
               return_type='both', patch_artist=True)
    plt.title('Variability between runs')
    plt.suptitle('\ninitPops = ' + str(initPops) + ', numRuns = ' + str(numRuns), x=0.85, y=0.92, ha='center', fontsize=10)
    plt.xlabel('Year')
    plt.ylabel('Population')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    # Plot box plots for each column
    plt.figure(figsize=(10, 6))
    df_NuwaraEliya.boxplot(grid=False, boxprops=dict(color='k'), medianprops=dict(color='k'),
               whiskerprops=dict(color='k'), flierprops=dict(color='k'),
               return_type='both', patch_artist=True)
    plt.title('Variability between runs')
    plt.suptitle('\ninitPops = ' + str(initPops) + ', numRuns = ' + str(numRuns), x=0.85, y=0.92, ha='center', fontsize=10)
    plt.xlabel('Year')
    plt.ylabel('Population')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    # Plot box plots for each column
    plt.figure(figsize=(10, 6))
    df_Jaffna.boxplot(grid=False, boxprops=dict(color='k'), medianprops=dict(color='k'),
               whiskerprops=dict(color='k'), flierprops=dict(color='k'),
               return_type='both', patch_artist=True)
    plt.title('Variability between runs')
    plt.suptitle('\ninitPops = ' + str(initPops) + ', numRuns = ' + str(numRuns), x=0.85, y=0.92, ha='center', fontsize=10)
    plt.xlabel('Year')
    plt.ylabel('Population')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
#
## TEST



## TEST
#
#  Test variability among runs for each immature introduction day (with a fixed year and month) and show as boxplots
if False:
    numRuns =  10
    initPops = 100
    initDay_first = 1
    initDay_last  = 7

    listOfLists = []
    for runNum in range(numRuns):
        print('runNum: ' + str(runNum))
        numAdults = []
        for initDay in range(initDay_first, initDay_last+1):
            df_vectorPop_dly = calc_monRecord_vectorPop("MERRA_IMERG", "Jaffna", 2020, 6, "intv00", initDay, initPops, initPops, initPops)
            numAdults.append(df_vectorPop_dly.iloc[-1]["pop_A"])
        listOfLists.append(numAdults)
    df = pd.DataFrame(listOfLists, columns=['day' + str(i) for i in range(initDay_first, initDay_last+1)])
    print(df)

    listOfLists = []
    for runNum in range(numRuns):
        print('runNum: ' + str(runNum))
        numAdults = []
        for initDay in range(initDay_first, initDay_last+1):
            df_vectorPop_dly = calc_monRecord_vectorPop("MERRA_IMERG", "Negombo", 2020, 6, "intv00", initDay, initPops, initPops, initPops)
            numAdults.append(df_vectorPop_dly.iloc[-1]["pop_A"])
        listOfLists.append(numAdults)
    df_Negombo = pd.DataFrame(listOfLists, columns=['day' + str(i) for i in range(initDay_first, initDay_last+1)])
    print(df_Negombo)

    listOfLists = []
    for runNum in range(numRuns):
        print('runNum: ' + str(runNum))
        numAdults = []
        for initDay in range(initDay_first, initDay_last+1):
            df_vectorPop_dly = calc_monRecord_vectorPop("MERRA_IMERG", "NuwaraEliya", 2020, 6, "intv00", initDay, initPops, initPops, initPops)
            numAdults.append(df_vectorPop_dly.iloc[-1]["pop_A"])
        listOfLists.append(numAdults)
    df_NuwaraEliya = pd.DataFrame(listOfLists, columns=['day' + str(i) for i in range(initDay_first, initDay_last+1)])
    print(df_NuwaraEliya)

    # Plot box plots for each column
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    df.boxplot(grid=False, boxprops=dict(color='k'), medianprops=dict(color='k'),
               whiskerprops=dict(color='k'), flierprops=dict(color='k'),
               return_type='both', patch_artist=True)
    plt.title('Variability between runs')
    plt.suptitle('\ninitPops = ' + str(initPops) + ', numRuns = ' + str(numRuns), x=0.85, y=0.92, ha='center', fontsize=10)
    plt.xlabel('Immatures introduction day')
    plt.ylabel('Population')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    # Plot box plots for each column
    plt.figure(figsize=(10, 6))
    df_Negombo.boxplot(grid=False, boxprops=dict(color='k'), medianprops=dict(color='k'),
               whiskerprops=dict(color='k'), flierprops=dict(color='k'),
               return_type='both', patch_artist=True)
    plt.title('Variability between runs')
    plt.suptitle('\ninitPops = ' + str(initPops) + ', numRuns = ' + str(numRuns), x=0.85, y=0.92, ha='center', fontsize=10)
    plt.xlabel('Immatures introduction day')
    plt.ylabel('Population')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

    # Plot box plots for each column
    plt.figure(figsize=(10, 6))
    df_NuwaraEliya.boxplot(grid=False, boxprops=dict(color='k'), medianprops=dict(color='k'),
               whiskerprops=dict(color='k'), flierprops=dict(color='k'),
               return_type='both', patch_artist=True)
    plt.title('Variability between runs')
    plt.suptitle('\ninitPops = ' + str(initPops) + ', numRuns = ' + str(numRuns), x=0.85, y=0.92, ha='center', fontsize=10)
    plt.xlabel('Immatures introduction day')
    plt.ylabel('Population')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
#
## TEST



## TEST
#
#  Test and compare various methods of estimating missing water temps.
if False:
    from sklearn.linear_model import LinearRegression
    import matplotlib.pyplot as plt

    filepath_in, _, _ = WHATCHEM_specs_to_fpath("MERRA_IMERG", "Jaffna", 2016, 4, "intv00") # 2020 6
    df_in = read_WHATCHEM_output(filepath_in) # note: converts the -9999.0 missing values into nan values

    # Compute daily survival factors and store these values
    waterTemp_hrly    = df_in["TW (C)"]
    groundTemp_hrly   = df_in["TG (C)"]
    airTemp_hrly      = df_in["TA (C)"]
    sunExp            = float(df_in.attrs["metadata"]["Shade"])

    airTemp, airTemp_min, airTemp_max          = convert_hrly2dlyMeanMinMax(airTemp_hrly, 'TA_min (C)', 'TA_max (C)')
    groundTemp, groundTemp_min, groundTemp_max = convert_hrly2dlyMeanMinMax(groundTemp_hrly, 'TG_min (C)', 'TG_max (C)')
    waterTemp, waterTemp_min, waterTemp_max    = convert_hrly2dlyMeanMinMax(waterTemp_hrly, 'TW_min (C)', 'TW_max (C)')

    # Estimate min/max water temps using Focks formulas
    waterTemp_min_est = 5.02 - 1.36*sunExp + 0.81*airTemp_min + 0.001*(airTemp_max**2)
    waterTemp_min_est = waterTemp_min_est.rename('TW_min_est (C)')
    waterTemp_max_est = 15.03 + 0.27*airTemp_min + 0.01*(airTemp_max**2) + 7.69*(sunExp**2)
    waterTemp_max_est = waterTemp_max_est.rename('TW_max_est (C)')
    

    # Fit linear regression model to hourly data
    df_4reg = pd.concat([airTemp_hrly, waterTemp_hrly], axis=1)
    df_4reg_noNaN = df_4reg.dropna() # remove rows that have nan values
    dates_missingData = df_4reg.index.difference(df_4reg_noNaN.index) # get dates that had nan values
    
    regression_model = LinearRegression()
    regression_model.fit(df_4reg_noNaN['TA (C)'].values.reshape(-1, 1), df_4reg_noNaN['TW (C)'].values.reshape(-1, 1))
    waterTemp_hrly_pred = regression_model.predict(airTemp_hrly.loc[dates_missingData].values.reshape(-1, 1)) # predict waterTemp for days with nan values
    waterTemp_hrly_adj = waterTemp_hrly.copy(deep=True)
    waterTemp_hrly_adj.loc[dates_missingData] = waterTemp_hrly_pred.flatten() # substitute in the predicted water temp values for the nan values
    waterTemp_hrly_adj = waterTemp_hrly_adj.rename('TW_adj (C)')


    # Fit linear regression model to daily mean/min/max data
    df_4reg_mean = pd.concat([airTemp, waterTemp], axis=1)
    df_4reg_mean_noNaN = df_4reg_mean.dropna() # remove rows that have nan values
    dates_missingData_mean = df_4reg_mean.index.difference(df_4reg_mean_noNaN.index) # get dates that had nan values
    regression_model = LinearRegression()
    regression_model.fit(df_4reg_mean_noNaN['TA (C)'].values.reshape(-1, 1), df_4reg_mean_noNaN['TW (C)'].values.reshape(-1, 1))
    waterTemp_pred = regression_model.predict(airTemp.loc[dates_missingData_mean].values.reshape(-1, 1)) # predict waterTemp for days with nan values
    waterTemp_adj = waterTemp.copy(deep=True)
    waterTemp_adj.loc[dates_missingData_mean] = waterTemp_pred.flatten() # substitute in the predicted water temp values for the nan values
    waterTemp_adj = waterTemp_adj.rename('TW_adj (C)')
    
    df_4reg_min = pd.concat([airTemp_min, waterTemp_min], axis=1)
    df_4reg_min_noNaN = df_4reg_min.dropna() # remove rows that have nan values
    dates_missingData_min = df_4reg_min.index.difference(df_4reg_min_noNaN.index) # get dates that had nan values
    regression_model = LinearRegression()
    regression_model.fit(df_4reg_min_noNaN['TA_min (C)'].values.reshape(-1, 1), df_4reg_min_noNaN['TW_min (C)'].values.reshape(-1, 1))
    waterTemp_min_pred = regression_model.predict(airTemp_min.loc[dates_missingData_min].values.reshape(-1, 1)) # predict waterTemp for days with nan values
    waterTemp_min_adj = waterTemp_min.copy(deep=True)
    waterTemp_min_adj.loc[dates_missingData_min] = waterTemp_min_pred.flatten() # substitute in the predicted water temp values for the nan values
    waterTemp_min_adj = waterTemp_min_adj.rename('TW_min_adj (C)')

    df_4reg_max = pd.concat([airTemp_max, waterTemp_max], axis=1)
    df_4reg_max_noNaN = df_4reg_max.dropna() # remove rows that have nan values
    dates_missingData_max = df_4reg_max.index.difference(df_4reg_max_noNaN.index) # get dates that had nan values
    regression_model = LinearRegression()
    regression_model.fit(df_4reg_max_noNaN['TA_max (C)'].values.reshape(-1, 1), df_4reg_max_noNaN['TW_max (C)'].values.reshape(-1, 1))
    waterTemp_max_pred = regression_model.predict(airTemp_max.loc[dates_missingData_max].values.reshape(-1, 1)) # predict waterTemp for days with nan values
    waterTemp_max_adj = waterTemp_max.copy(deep=True)
    waterTemp_max_adj.loc[dates_missingData_max] = waterTemp_max_pred.flatten() # substitute in the predicted water temp values for the nan values
    waterTemp_max_adj = waterTemp_max_adj.rename('TW_max_adj (C)')

    colors = ["blue", "orange", "brown", "cornflowerblue"]

    plt.figure(figsize=(10, 6))  # Set figure size if needed


    ## Estimating based on lin reg against daily mean/min/max air temp
    if True:
        df_mean = pd.concat([waterTemp, airTemp, groundTemp, waterTemp_adj], axis=1)
        df_min = pd.concat([waterTemp_min, airTemp_min, groundTemp_min, waterTemp_min_adj], axis=1)
        df_max = pd.concat([waterTemp_max, airTemp_max, groundTemp_max, waterTemp_max_adj], axis=1)
        for i, col in enumerate(df_mean.columns):
            if i == 3:
                plt.plot(df_mean.loc[dates_missingData_mean].index, df_mean[col].loc[dates_missingData_mean], label=f'{col}', color=colors[i])
            else:
                plt.plot(df_mean.index, df_mean[col], label=f'{col}', color=colors[i])
        for i, col in enumerate(df_min.columns):
            if i == 3:
                plt.plot(df_min.loc[dates_missingData_min].index, df_min[col].loc[dates_missingData_min], label=f'{col}', linestyle='--', color=colors[i])
            else:
                plt.plot(df_min.index, df_min[col], label=f'{col}', linestyle='--', color=colors[i])
        for i, col in enumerate(df_max.columns):
            if i == 3:
                plt.plot(df_max.loc[dates_missingData_max].index, df_max[col].loc[dates_missingData_max], label=f'{col}', linestyle='--', color=colors[i])
            else:
                plt.plot(df_max.index, df_max[col], label=f'{col}', linestyle='--', color=colors[i])


    ## Estimating based on lin reg against hrly air temp
    if False:
        df_hrly =pd.concat([waterTemp_hrly, airTemp_hrly, groundTemp_hrly, waterTemp_hrly_adj], axis=1)
        for i, col in enumerate(df_hrly.columns):
            if i == 3:
                plt.plot(df_hrly.loc[dates_missingData].index, df_hrly[col].loc[dates_missingData], label=f'{col}', color=colors[i])
            else:
                plt.plot(df_hrly.index, df_hrly[col], label=f'{col}', color=colors[i])


    ## Estimating based on Focks formulas
    if False:
        df_min = pd.concat([waterTemp_min, airTemp_min, groundTemp_min, waterTemp_min_est], axis=1)
        df_max = pd.concat([waterTemp_max, airTemp_max, groundTemp_max, waterTemp_max_est], axis=1)
        plt.figure(figsize=(10, 6))  # Set figure size if needed
        for i, col in enumerate(df_min.columns):
            plt.plot(df_min.index, df_min[col], label=f'{col}', color=colors[i])
        for i, col in enumerate(df_max.columns):
            plt.plot(df_max.index, df_max[col], label=f'{col}', color=colors[i])

    if True:
        # Add labels and title
        plt.xlabel('Date')
        plt.ylabel('Temperature')
        plt.title('')
        plt.legend()

        # Show the plot
        plt.show()

    #a = calc_survProb_temp_dly(waterTemp_hrly, groundTemp_hrly, waterHeight_hrly, "Eggs")
    #b = calc_survProb_desicc_dly(waterHeight_hrly, vpd_hrly, 0.5, "Eggs")
    #c = calc_survProb_dly(waterTemp_hrly, groundTemp_hrly, waterHeight_hrly, vpd_hrly, 0.5, "Eggs")
    #d = calc_devRate(waterTemp_hrly, groundTemp_hrly, waterHeight_hrly, "Eggs")
    #df_dev_hrly, df_dev_surv_dly = calc_monRecord_dev_surv("MERRA_IMERG", "Jaffna", 2020, 6, "intv00")
    #print(df_dev_surv_dly)
    #print(df_dev_surv_dly.index)
    #print(df_dev_surv_dly.attrs["metadata"])
    #save_monRecord_dev_surv("MERRA_IMERG", "Jaffna", 2020, 6, "intv00")
    #save_monRecord_vectorPop("MERRA_IMERG", "Jaffna", 2020, 6, "intv00", 100, 100, 100)
#   
## TEST



## TEST
#
# Check all monthly WHATCH'EM simulations and print out which ones have missing water temp data (and how much is missing).
if False:
    dataSrc   = "MERRA_IMERG"
    locs      = ["Negombo", "Jaffna", "NuwaraEliya"]
    yrStrt    = 2001
    yrStop    = 2020
    years     = range(yrStrt, yrStop+1) # note: first number is included, last number is not
    months    = range(1,13)             # note: first number is included, last number is not
    intv      = "intv00"
    initEggs   = 100
    initLarvae = 100
    initPupae  = 100
    numWithMissing = 0
    numWithMissing_min50p = 0
    for loc in locs:
        for year in years:
            for month in months:
                filepath_in, subDir, fileStr_allSpecs = WHATCHEM_specs_to_fpath(dataSrc, loc, year, month, intv)
                df_in = read_WHATCHEM_output(filepath_in) # note: converts the -9999.0 missing values into nan values
                missingValCnt = df_in['TW (C)'].isnull().sum()
                totalCnt = len(df_in['TW (C)'])
                if missingValCnt > 0:
                    missingValPcnt = round(100*missingValCnt/totalCnt)
                    print(loc + ' ' + str(year) + ' ' + str(month) + ": " + str(missingValCnt) + '/' + str(totalCnt) + ' (' + str(missingValPcnt) + ')')
                    numWithMissing += 1
                    if missingValPcnt >= 50:
                        numWithMissing_min50p += 1
    print('Number of months with any missing data: ' + str(numWithMissing))
    print('Number of months missing >= 50% of data: ' + str(numWithMissing_min50p))
#
## TEST


