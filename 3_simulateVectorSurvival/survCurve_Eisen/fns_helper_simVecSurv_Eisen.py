#*******************************************************
# fns_helper_sinVecSurv_Eisen.py
#*******************************************************
# Helper functions for processing the survival data from
# Eisen et al. 2014 into survival curves we can use for
# our modeling.
#
# Revision history
# 2025-04-17, CNY: Updated calc_survProb_temp_dly_eisenQuadFit()
#                  to clip values at 1 instead of scaling them

# Import packages
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit


#####################################################################################################
# New functions for processing Eisen et al. 2014 survival curves
#####################################################################################################


# This function converts the whole-stage survival rates from Eisen et al. 2014 (defined by the dictionaries above)
# to daily survival rates. It does so in a few steps.
# (1) For each temperature in the survival rates df, find the corresponding hourly development rate from Focks et al. 1993
#     (using calc_devRate() from simVectorSurvival.py)
# (2) Using the hourly development rates, estimate stage duration based on the assumption that a stage is complete when
#     cumulative develpment > 0.95.
# (3) Using the stage durations, convert the whole-stage survival rates to daily survival rates.
# The showSteps boolean is used to print all the intermediate calculated values to the console
def survCurve_Eisen_convert2dly(lifeStage, showSteps):

    # Define survival vs water temperatuee datapoints from Eisen et al. 2014.
    # (See photoshop and excel files for details on how these numbers were extracted from the paper's figures.)
    df_survVsTemp_eggs_lifeStage = pd.DataFrame({
        "TW (C)": [0.98, 6.98, 15.53, 15.91, 21.01, 21.91, 24.99, 26.64, 27.92, 30.92, 32.12, 34.97, 36.02, 37.75],
        "Pr_surv_T_lifeStage": [0.00, 0.00, 0.78, 0.81, 0.88, 0.94, 0.96, 0.91, 0.96, 0.83, 0.90, 0.48, 0.00, 0.00]
    })
    df_survVsTemp_larvae_lifeStage = pd.DataFrame({
        "TW (C)": [ 9.07, 10.15, 11.14, 12.05, 13.04, 14.19, 15.68, 16.09, 16.09, 16.09,
                    20.05, 20.05, 20.05, 20.05, 20.05, 24.02, 24.02, 24.51, 25.01, 25.01,
                    25.01, 25.50, 26.49, 26.99, 27.98, 27.98, 27.98, 29.46, 29.96, 29.96,
                    29.96, 29.96, 32.02, 32.02, 32.02, 32.52, 34.09, 34.09, 34.50, 34.99,
                    36.07, 36.07, 36.07, 37.97, 39.95, 43.25, 46.14, 48.78, 51.59],
        "Pr_surv_T_lifeStage": [0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 0.46, 0.10, 0.76, 0.77,
                                0.54, 0.56, 0.70, 0.74, 0.92, 0.88, 0.97, 0.94, 0.82, 0.86,
                                0.92, 0.73, 0.98, 0.93, 0.60, 0.96, 0.99, 0.97, 0.86, 0.89,
                                0.95, 1.00, 0.85, 0.96, 0.98, 0.92, 0.92, 0.95, 0.80, 0.06,
                                0.01, 0.37, 0.89, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
    })
    df_survVsTemp_pupae_lifeStage = pd.DataFrame({
        "TW (C)": [7.04, 9.04, 10.00, 11.04, 12.09, 12.09, 13.04, 16.00, 18.00, 20.09,
                23.04, 24.09, 25.65, 26.09, 27.04, 28.09, 28.09, 29.04, 30.09, 32.09,
                34.00, 36.00, 38.00, 43.30, 44.09, 45.04, 46.09, 47.04, 48.09, 49.04,
                51.74],
        "Pr_surv_T_lifeStage": [0.00, 0.00, 0.00, 0.01, 0.00, 0.54, 0.90, 0.92, 1.00, 0.97,
                                1.00, 0.84, 0.97, 1.00, 0.98, 0.98, 1.00, 1.00, 0.99, 0.79,
                                0.95, 0.86, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                0.00]
    })

    if lifeStage == "Eggs":
        df_survVsTemp_lifestage = df_survVsTemp_eggs_lifeStage
    elif lifeStage == "Larvae":
        df_survVsTemp_lifestage = df_survVsTemp_larvae_lifeStage
    elif lifeStage == "Pupae":
        df_survVsTemp_lifestage = df_survVsTemp_pupae_lifeStage

    # Extract water temperature and survival values from the relevant survival curve
    waterTemp_vals = df_survVsTemp_lifestage["TW (C)"].values
    surv_lifeStage_vals = df_survVsTemp_lifestage["Pr_surv_T_lifeStage"].values

    # Step 1: Get hourly development rates (from Focks et al. 1993)
    devRate_hourly = np.array([calc_devRate(waterTemp, lifeStage) for waterTemp in waterTemp_vals])

    # Step 2: Estimate duration in hours until cumulative development > 0.95
    # duration_hours = 0.95 / devRate_hourly, avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        duration_lifeStage = np.where(devRate_hourly > 0, 0.95 / devRate_hourly, np.nan) # duration in hrs
    duration_lifeStage = duration_lifeStage / 24.0 # convert from hours to days

    # Step 3: Convert whole-stage survival to daily survival using exponential decay model
    with np.errstate(divide='ignore', invalid='ignore'):
        surv_daily = np.where(duration_lifeStage > 0, surv_lifeStage_vals ** (1 / duration_lifeStage), np.nan)

    # Create and store the new DataFrame
    df_survVsTemp_daily = pd.DataFrame({
        "TW (C)": waterTemp_vals,
        "Pr_surv_T_daily": surv_daily
    })

    df_summary = pd.DataFrame({
        "TW (C)": waterTemp_vals,
        "Pr_surv_T_lifeStage": surv_lifeStage_vals,
        "devRate (hr^-1)": devRate_hourly,
        "duration_lifeStage (days)": duration_lifeStage,
        "Pr_surv_T_daily": surv_daily
    })

    if showSteps:
        print("\nSummary for " + lifeStage + ":")
        print(df_summary)

    return df_survVsTemp_daily




# Fit the Eisen survival curve data with a quadratic function of the form -c(T–T0)(T–Tm), 
# where T0 and Tm are the minimum and maximum temperature for non-zero survival and c is a
# positive rate constant.
# This is the kind of fit used in Mordecai et al. 2017 for fitting egg-to-adult survival data.
def getParams_quadFit_survCurve_Eisen(df_survVsTemp_daily, lifeStage, doPlot, ax=None):
    # Drop rows with NaN and sort by temperature
    df_all = df_survVsTemp_daily.dropna().sort_values("TW (C)")
    temps_all = df_all["TW (C)"].values
    surv_all = df_all["Pr_surv_T_daily"].values

    # Identify indices where survival is nonzero
    nonzero_inds = np.where(surv_all > 0)[0]
    if len(nonzero_inds) == 0:
        raise ValueError("No nonzero survival values found.")

    # Identify first and last nonzero survival indices
    first_nonzero = nonzero_inds[0]
    last_nonzero = nonzero_inds[-1]

    # T0 is the temperature just before the first nonzero value if that value is 0
    if first_nonzero > 0 and surv_all[first_nonzero - 1] == 0:
        T0 = temps_all[first_nonzero - 1]
    else:
        T0 = temps_all[first_nonzero]

    # Tm is the temperature just after the last nonzero value if that value is 0
    if last_nonzero < len(surv_all) - 1 and surv_all[last_nonzero + 1] == 0:
        Tm = temps_all[last_nonzero + 1]
    else:
        Tm = temps_all[last_nonzero]

    # Filter to only the nonzero survival values for fitting
    temps = temps_all[nonzero_inds]
    surv = surv_all[nonzero_inds]

    # Define the fit function: f(T) = -c*(T-T0)*(T-Tm)
    def quad_model(T, c):
        return -c*(T-T0)*(T-Tm)

    # Initial guess for c
    c0 = 0.01

    # Fit the curve
    popt, pcov = curve_fit(quad_model, temps, surv, p0=[c0])
    c_fit = popt[0]


    # Plotting
    if doPlot:

        print(f"Fitted c = {c_fit:.5f}")
        print(f"T0 = {T0:.2f}, Tm = {Tm:.2f}")

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        T_fit = np.linspace(T0, Tm, 300)
        surv_fit = quad_model(T_fit, c_fit)

        ax.plot(temps, surv, 'o', label='Empirical daily survival')
        ax.plot(T_fit, surv_fit, '-', label=f'Fitted: -{c_fit:.4f}(T - {T0:.1f})(T - {Tm:.1f})')
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Daily Survival Probability")
        ax.set_title(f"Quadratic Fit to Daily Survival Curve ({lifeStage})")
        ax.legend()
        ax.grid(True)

    return c_fit, T0, Tm




# Fit the Eisen survival curve data with the asymmetric Briere function of the form cT(T–T0)(Tm−T)^(1/2), 
# where T0 and Tm are the minimum and maximum temperature for non-zero survival and c is a
# positive rate constant.
# This is a kind of fit used in Mordecai et al. 2017 for fitting thermal response curves (but not the one
# they used for fitting the egg-to-survival data).
def getParams_briereFit_survCurve_Eisen(df_survVsTemp_daily, lifeStage, doPlot, ax=None):
    # Drop rows with NaN and sort by temperature
    df_all = df_survVsTemp_daily.dropna().sort_values("TW (C)")
    temps_all = df_all["TW (C)"].values
    surv_all = df_all["Pr_surv_T_daily"].values

    # Identify indices where survival is nonzero
    nonzero_inds = np.where(surv_all > 0)[0]
    if len(nonzero_inds) == 0:
        raise ValueError("No nonzero survival values found.")

    # Identify first and last nonzero survival indices
    first_nonzero = nonzero_inds[0]
    last_nonzero = nonzero_inds[-1]

    # T0 is the temperature just before the first nonzero value if that value is 0
    if first_nonzero > 0 and surv_all[first_nonzero - 1] == 0:
        T0 = temps_all[first_nonzero - 1]
    else:
        T0 = temps_all[first_nonzero]

    # Tm is the temperature just after the last nonzero value if that value is 0
    if last_nonzero < len(surv_all) - 1 and surv_all[last_nonzero + 1] == 0:
        Tm = temps_all[last_nonzero + 1]
    else:
        Tm = temps_all[last_nonzero]

    # Filter to only the nonzero survival values for fitting
    temps = temps_all[nonzero_inds]
    surv = surv_all[nonzero_inds]

    # Define the fit function: f(T) = c*T*(T-T0)*(Tm-T)^(1/2)
    def briere_model(T, c):
        return c*T*(T-T0)*(Tm-T)**(0.5)

    # Initial guess for c
    c0 = 0.01

    # Fit the curve
    popt, pcov = curve_fit(briere_model, temps, surv, p0=[c0])
    c_fit = popt[0]

    # Plotting
    if doPlot:

        print(f"Fitted c = {c_fit:.5f}")
        print(f"T0 = {T0:.2f}, Tm = {Tm:.2f}")

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        T_fit = np.linspace(T0, Tm, 300)
        surv_fit = briere_model(T_fit, c_fit)

        ax.plot(temps, surv, 'o', label='Empirical daily survival')
        ax.plot(T_fit, surv_fit, '-', label=f'Fitted: -{c_fit:.4f}(T - {T0:.1f})(T - {Tm:.1f})')
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Daily Survival Probability")
        ax.set_title(f"Briere Fit to Daily Survival Curve ({lifeStage})")
        ax.legend()
        ax.grid(True)

    return c_fit, T0, Tm




# Plot all the survival curves together to compare - quadratic fit edition
#  (1) Focks et al. 1993 curves
#  (2) Eisen et al. 2014 datapoints
#  (3) Eisen et al. 2014 curve fit (quadratic, clipped to max 1)
#  (4) Eisen et al. 2014 curve fit (quadratic, scaled to max 1)
def plt_compare_survCurves_eisenQuadFit(lifeStage):

    # Construct survival curve from Focks et al. 1993

    # Define dataframe of water temperature thresholds for survival (in deg C) from Focks et al.
    # The four thresholds (T0, T1, T2, T3) define daily survival probabilities (Pr_surv_Tmin, Pr_surv_Tmax) 
    # based on daily min/max water temps (waterTemp_min, waterTemp_max).
    # See Sec. S2.7 of Magori et al. 2009 for more details.
    tempThresholds_dict = {"T0": [-14, 5, 5],
                           "T1": [-6, 10, 10],
                           "T2": [30, 39, 39],
                           "T3": [47, 44, 44]}
    tempThresholds_df = pd.DataFrame(tempThresholds_dict, index=["Eggs", "Larvae", "Pupae"])
    T0, T1, T2, T3 = tempThresholds_df.loc[lifeStage] # get thresholds for specified life stage

    temps = np.arange(-15., 50.+0.25, 0.25) # create linearly spaced temp values

    surv = np.full_like(temps, 0.05, dtype=float)  # initialize survival array; default to 0.05

    # modify survival array based on masks for different temperature regimes
    mask1 = (temps > T0) & (temps < T1)   # Linearly increase from 0.05 to 1.0 between T0 and T1
    mask2 = (temps >= T1) & (temps <= T2) # Flat at 1.0 between T1 and T2
    mask3 = (temps > T2) & (temps < T3)   # Linearly decrease from 1.0 to 0.05 between T2 and T3
    surv[mask1] = 0.05 + 0.95*((temps[mask1] - T0) / (T1 - T0))
    surv[mask2] = 1.0
    surv[mask3] = 1.0 - 0.95*((temps[mask3] - T2) / (T3 - T2))

    df_survCurve_Focks = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv})

    # Construct survival curve from Eisen et al. 2014
    df_survCurve_Eisen = survCurve_Eisen_convert2dly(lifeStage, showSteps=False)

    # Construct quadratic fit to Eisen et al. 2014
    quadFit_c, quadFit_T0, quadFit_Tm = getParams_quadFit_survCurve_Eisen(df_survCurve_Eisen, lifeStage, doPlot=False)
    quadFit_lowerBoundT = max(0, quadFit_T0)   # as in Mordecai et al. 2017, we assume temps between 0C are lethal (regardless of what the fit says)
    quadFit_upperBoundT = min(45, quadFit_Tm)  # as in Mordecai et al. 2017, we assume temps between 45C are lethal (regardless of what the fit says)
    surv_quadFit = np.full_like(temps, 0.00, dtype=float)  # initialize survival array; default to 0.00
    mask_quad = (temps >= quadFit_lowerBoundT) & (temps <= quadFit_upperBoundT) # range over which we use the fit (outside of this range the surv val is zero)
    surv_quadFit[mask_quad] = -quadFit_c*(temps[mask_quad]-quadFit_T0)*(temps[mask_quad]-quadFit_Tm)

    # Do one version where we clip the min/max values to be between 0 and 1
    surv_quadFit_clipped = np.clip(surv_quadFit, 0, 1) # clip values so that they're between 0 and 1
    df_survCurve_eisenQuadFit_clipped = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv_quadFit_clipped})

    # Do another version where we scale the values so that the max value is 1
    surv_quadFit_maxVal = max(surv_quadFit)
    surv_quadFit_scaled = np.clip(surv_quadFit/surv_quadFit_maxVal, 0, 1) # scale to max value, but also clip in case the division causes floating point errors (values slightly below 0 or above 1)
    df_survCurve_eisenQuadFit_scaled = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv_quadFit_scaled})

    # Plot survival curves against each other
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot Focks curve
    ax.plot(df_survCurve_Focks["TW (C)"], df_survCurve_Focks["Pr_surv_T_daily"],
            label="Focks et al. 1993 (for daily min/max temps)", color='blue', linestyle='--')
    # Plot Eisen curve quadratic fit (clipped)
    ax.plot(df_survCurve_eisenQuadFit_clipped["TW (C)"], df_survCurve_eisenQuadFit_clipped["Pr_surv_T_daily"],
            label="Eisen et al. 2014 quadratic fit (clipped)", color='#CC0000', linestyle='-')
    # Plot Eisen curve quadratic fit (scaled)
    ax.plot(df_survCurve_eisenQuadFit_scaled["TW (C)"], df_survCurve_eisenQuadFit_scaled["Pr_surv_T_daily"],
            label="Eisen et al. 2014 quadratic fit (scaled)", color='#550000', linestyle='-')
    # Plot Eisen data points
    ax.plot(df_survCurve_Eisen["TW (C)"], df_survCurve_Eisen["Pr_surv_T_daily"],
            label="Eisen et al. 2014 data points", color='red', marker='o', linestyle='')

    # Labeling
    ax.set_title(f"Daily Survival Curve - {lifeStage}")
    ax.set_xlabel("Water Temperature (°C)")
    ax.set_ylabel("Daily Survival Probability")
    ax.set_ylim(0, 1.05)
    ax.grid(True)
    ax.legend()

    # Save figure
    fig.tight_layout()
    fig.savefig(f"survCurve_FocksVsEisen_quadFit_{lifeStage}.png", dpi=300)
    plt.close(fig)  # Close the figure to prevent display in some environments




# Plot all the survival curves together to compare - quadratic fit edition
#  (1) Focks et al. 1993 curves
#  (2) Eisen et al. 2014 datapoints
#  (3) Eisen et al. 2014 curve fit (quadratic, clipped to max 1)
def plt_compare_survCurves_eisenQuadFit_onlyClip(lifeStage):

    # Construct survival curve from Focks et al. 1993

    # Define dataframe of water temperature thresholds for survival (in deg C) from Focks et al.
    # The four thresholds (T0, T1, T2, T3) define daily survival probabilities (Pr_surv_Tmin, Pr_surv_Tmax) 
    # based on daily min/max water temps (waterTemp_min, waterTemp_max).
    # See Sec. S2.7 of Magori et al. 2009 for more details.
    tempThresholds_dict = {"T0": [-14, 5, 5],
                           "T1": [-6, 10, 10],
                           "T2": [30, 39, 39],
                           "T3": [47, 44, 44]}
    tempThresholds_df = pd.DataFrame(tempThresholds_dict, index=["Eggs", "Larvae", "Pupae"])
    T0, T1, T2, T3 = tempThresholds_df.loc[lifeStage] # get thresholds for specified life stage

    temps = np.arange(-15., 50.+0.25, 0.25) # create linearly spaced temp values

    surv = np.full_like(temps, 0.05, dtype=float)  # initialize survival array; default to 0.05

    # modify survival array based on masks for different temperature regimes
    mask1 = (temps > T0) & (temps < T1)   # Linearly increase from 0.05 to 1.0 between T0 and T1
    mask2 = (temps >= T1) & (temps <= T2) # Flat at 1.0 between T1 and T2
    mask3 = (temps > T2) & (temps < T3)   # Linearly decrease from 1.0 to 0.05 between T2 and T3
    surv[mask1] = 0.05 + 0.95*((temps[mask1] - T0) / (T1 - T0))
    surv[mask2] = 1.0
    surv[mask3] = 1.0 - 0.95*((temps[mask3] - T2) / (T3 - T2))

    df_survCurve_Focks = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv})

    # Construct survival curve from Eisen et al. 2014
    df_survCurve_Eisen = survCurve_Eisen_convert2dly(lifeStage, showSteps=False)

    # Construct quadratic fit to Eisen et al. 2014
    quadFit_c, quadFit_T0, quadFit_Tm = getParams_quadFit_survCurve_Eisen(df_survCurve_Eisen, lifeStage, doPlot=False)
    quadFit_lowerBoundT = max(0, quadFit_T0)   # as in Mordecai et al. 2017, we assume temps between 0C are lethal (regardless of what the fit says)
    quadFit_upperBoundT = min(45, quadFit_Tm)  # as in Mordecai et al. 2017, we assume temps between 45C are lethal (regardless of what the fit says)
    surv_quadFit = np.full_like(temps, 0.00, dtype=float)  # initialize survival array; default to 0.00
    mask_quad = (temps >= quadFit_lowerBoundT) & (temps <= quadFit_upperBoundT) # range over which we use the fit (outside of this range the surv val is zero)
    surv_quadFit[mask_quad] = -quadFit_c*(temps[mask_quad]-quadFit_T0)*(temps[mask_quad]-quadFit_Tm)

    # Do one version where we clip the min/max values to be between 0 and 1
    surv_quadFit_clipped = np.clip(surv_quadFit, 0, 1) # clip values so that they're between 0 and 1
    df_survCurve_eisenQuadFit_clipped = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv_quadFit_clipped})

    # Do another version where we scale the values so that the max value is 1
    surv_quadFit_maxVal = max(surv_quadFit)
    surv_quadFit_scaled = np.clip(surv_quadFit/surv_quadFit_maxVal, 0, 1) # scale to max value, but also clip in case the division causes floating point errors (values slightly below 0 or above 1)
    df_survCurve_eisenQuadFit_scaled = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv_quadFit_scaled})

    # Plot survival curves against each other
    fig, ax = plt.subplots(figsize=(8, 6))

    # Define color mapping
    color_vectBio_E = '#f394c5' # color for vector biology (eggs) variables
    color_vectBio_L = '#b9216e' # color for vector biology (larvae) variables
    color_vectBio_P = '#450c29' # color for vector biology (pupae) variables
    color_vectBio_P_lighter = mcolors.to_hex([c + (1 - c) * 0.1 for c in mcolors.to_rgb(color_vectBio_P)]) # make color slightly lighter
    hue_colors_surv = {'Eggs': color_vectBio_E, 'Larvae': color_vectBio_L, 'Pupae': color_vectBio_P_lighter}

    # Plot Eisen data points
    ax.scatter(df_survCurve_Eisen["TW (C)"], df_survCurve_Eisen["Pr_surv_T_daily"],
               label="Eisen et al. (2014) data", color=hue_colors_surv[lifeStage], marker='o', s=50, edgecolor='black', linewidth=0.8, zorder=3)
    # Plot Eisen curve quadratic fit (clipped)
    ax.plot(df_survCurve_eisenQuadFit_clipped["TW (C)"], df_survCurve_eisenQuadFit_clipped["Pr_surv_T_daily"],
            label="Eisen et al. (2014) fit", color=hue_colors_surv[lifeStage], linestyle='-', lw=4, zorder=2)
    # Plot Focks curve
    ax.plot(df_survCurve_Focks["TW (C)"], df_survCurve_Focks["Pr_surv_T_daily"],
            label="Focks et al. (1993)", color='#666666', linestyle='--', lw=3, zorder=1)

    # Adjust axes
    xmin_surv_temp  = -15
    xmax_surv_temp  = 50
    xstep_surv_temp = 0.25
    ax_x_ticks = np.round(np.arange(xmin_surv_temp+5, xmax_surv_temp + 0.01, 10),0)
    ax.set_xticks(ax_x_ticks)
    ax.tick_params(axis='both', which='minor', length=7, color='gray', direction='in', width=0.8, labelsize=16, left=True, right=True, bottom=True, top=True) # adjust minor ticks for y-axis
    ax.tick_params(axis='both', which='major', length=11, color='gray', direction='in', width=1.2, labelsize=16, left=True, right=True, bottom=True, top=True) # adjust major ticks for y-axis
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax_y_ticks = np.round(np.arange(0, 1+0.01, 0.2),1)
    ax.set_yticks(ax_y_ticks)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Labeling
    ax.set_title(f"{lifeStage}", fontsize=18)
    ax.set_xlabel('Water Temperature (°C)', fontsize=16)
    ax.set_ylabel('Survival factor (day$^{-1}$)', fontsize=16)
    ax.set_ylim(0, 1.05)
    ax.legend(loc='upper left', bbox_to_anchor=(0.03, 0.95))

    # Save figure
    fig.tight_layout()
    fig.savefig(f"survCurve_FocksVsEisen_quadFit_onlyClip_{lifeStage}.png", dpi=300)
    plt.close(fig)  # Close the figure to prevent display in some environments





# Plot all the survival curves together to compare - Briere fit edition
#  (1) Focks et al. 1993 curves
#  (2) Eisen et al. 2014 datapoints
#  (3) Eisen et al. 2014 curve fit (Briere, clipped to max 1)
#  (4) Eisen et al. 2014 curve fit (Briere, scaled to max 1)
def plt_compare_survCurves_eisenBriereFit(lifeStage):

    # Construct survival curve from Focks et al. 1993

    # Define dataframe of water temperature thresholds for survival (in deg C) from Focks et al.
    # The four thresholds (T0, T1, T2, T3) define daily survival probabilities (Pr_surv_Tmin, Pr_surv_Tmax) 
    # based on daily min/max water temps (waterTemp_min, waterTemp_max).
    # See Sec. S2.7 of Magori et al. 2009 for more details.
    tempThresholds_dict = {"T0": [-14, 5, 5],
                           "T1": [-6, 10, 10],
                           "T2": [30, 39, 39],
                           "T3": [47, 44, 44]}
    tempThresholds_df = pd.DataFrame(tempThresholds_dict, index=["Eggs", "Larvae", "Pupae"])
    T0, T1, T2, T3 = tempThresholds_df.loc[lifeStage] # get thresholds for specified life stage

    temps = np.arange(-15., 50.+0.25, 0.25) # create linearly spaced temp values

    surv = np.full_like(temps, 0.05, dtype=float)  # initialize survival array; default to 0.05

    # modify survival array based on masks for different temperature regimes
    mask1 = (temps > T0) & (temps < T1)   # Linearly increase from 0.05 to 1.0 between T0 and T1
    mask2 = (temps >= T1) & (temps <= T2) # Flat at 1.0 between T1 and T2
    mask3 = (temps > T2) & (temps < T3)   # Linearly decrease from 1.0 to 0.05 between T2 and T3
    surv[mask1] = 0.05 + 0.95*((temps[mask1] - T0) / (T1 - T0))
    surv[mask2] = 1.0
    surv[mask3] = 1.0 - 0.95*((temps[mask3] - T2) / (T3 - T2))

    df_survCurve_Focks = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv})

    # Construct survival curve from Eisen et al. 2014
    df_survCurve_Eisen = survCurve_Eisen_convert2dly(lifeStage, showSteps=False)

    # Construct Briere fit to Eisen et al. 2014
    briereFit_c, briereFit_T0, briereFit_Tm = getParams_briereFit_survCurve_Eisen(df_survCurve_Eisen, lifeStage, doPlot=False)
    briereFit_lowerBoundT = max(0, briereFit_T0)   # as in Mordecai et al. 2017, we assume temps between 0C are lethal (regardless of what the fit says)
    briereFit_upperBoundT = min(45, briereFit_Tm)  # as in Mordecai et al. 2017, we assume temps between 45C are lethal (regardless of what the fit says)
    surv_briereFit = np.full_like(temps, 0.00, dtype=float)  # initialize survival array; default to 0.00
    mask_briere = (temps >= briereFit_lowerBoundT) & (temps <= briereFit_upperBoundT) # range over which we use the fit (outside of this range the surv val is zero)
    surv_briereFit[mask_briere] = briereFit_c*temps[mask_briere]*(temps[mask_briere]-briereFit_T0)*(briereFit_Tm-temps[mask_briere])**(0.5)
    
    # Do one version where we clip the min/max values to be between 0 and 1
    surv_briereFit_clipped = np.clip(surv_briereFit, 0, 1) # clip values so that they're between 0 and 1
    df_survCurve_eisenBriereFit_clipped = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv_briereFit_clipped})

    # Do another version where we scale the values so that the max value is 1
    surv_briereFit_maxVal = max(surv_briereFit)
    surv_briereFit_scaled = np.clip(surv_briereFit/surv_briereFit_maxVal, 0, 1) # scale to max value, but also clip in case the division causes floating point errors (values slightly below 0 or above 1)
    df_survCurve_eisenBriereFit_scaled = pd.DataFrame({"TW (C)": temps, "Pr_surv_T_daily": surv_briereFit_scaled})

    # Plot survival curves against each other
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot Focks curve
    ax.plot(df_survCurve_Focks["TW (C)"], df_survCurve_Focks["Pr_surv_T_daily"],
            label="Focks et al. 1993 (for daily min/max temps)", color='blue', linestyle='--')
    # Plot Eisen curve Briere fit (clipped)
    ax.plot(df_survCurve_eisenBriereFit_clipped["TW (C)"], df_survCurve_eisenBriereFit_clipped["Pr_surv_T_daily"],
            label="Eisen et al. 2014 Briere fit (clipped)", color='#CC0000', linestyle='-')
    # Plot Eisen curve Briere fit (scaled)
    ax.plot(df_survCurve_eisenBriereFit_scaled["TW (C)"], df_survCurve_eisenBriereFit_scaled["Pr_surv_T_daily"],
            label="Eisen et al. 2014 Briere fit (scaled)", color='#550000', linestyle='-')
    # Plot Eisen data points
    ax.plot(df_survCurve_Eisen["TW (C)"], df_survCurve_Eisen["Pr_surv_T_daily"],
            label="Eisen et al. 2014 data points", color='red', marker='o', linestyle='')

    # Labeling
    ax.set_title(f"Daily Survival Curve - {lifeStage}")
    ax.set_xlabel("Water Temperature (°C)")
    ax.set_ylabel("Daily Survival Probability")
    ax.set_ylim(0, 1.05)
    ax.grid(True)
    ax.legend()

    # Save figure
    fig.tight_layout()
    fig.savefig(f"survCurve_FocksVsEisen_briereFit_{lifeStage}.png", dpi=300)
    plt.close(fig)  # Close the figure to prevent display in some environments






#####################################################################################################
# Old functions, some modified to adapt to Eisen et al. 2014 survival curves
#####################################################################################################

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




## calc_survProb_temp_dly_eisenQuadFit()
# Calculate daily survival probability based on the mean temperature for the day,
# based on the QUADRATIC fitted survival curve of Eisen et al. 2014.
# Following Mordecai et al. 2017, we assume temps <=0C or >=45C are lethal, regardless of what the fitted survival curve says.
# We also clip the survival curve so the maximum value is 1.
#  Args:
#    waterTemp_mean - daily mean water temperature [deg C]
#    lifeStage      - life stage of immatures (options: "Eggs", "Larvae", "Pupae")
#  Output:
#    Pr_surv_T      - daily survival probability given temperature
def calc_survProb_temp_dly_eisenQuadFit(waterTemp_mean, lifeStage):

    # Compute parameters of the quadratic fit and define the fit's functional form
    df_survCurve_Eisen = survCurve_Eisen_convert2dly(lifeStage, showSteps=False)
    quadFit_c, quadFit_T0, quadFit_Tm = getParams_quadFit_survCurve_Eisen(df_survCurve_Eisen, lifeStage, doPlot=False)
    def fn_quadFit(c, T, T0, Tm):
        return -c*(T-T0)*(T-Tm)
    
    # Set temperature limits
    quadFit_lowerBoundT = max(0, quadFit_T0)   # as in Mordecai et al. 2017, we assume temps below 0C are lethal (regardless of what the fit says)
    quadFit_upperBoundT = min(45, quadFit_Tm)  # as in Mordecai et al. 2017, we assume temps above 45C are lethal (regardless of what the fit says)

    # Create boolean arrays that indicate the temperature ranges that are lethal vs nonlethal
    T_regime_nonlethal = (waterTemp_mean >= quadFit_lowerBoundT) & (waterTemp_mean <= quadFit_upperBoundT)
    T_regime_lethal = ~T_regime_nonlethal

    # Initialize output array
    Pr_surv_T = np.full_like(waterTemp_mean, np.nan)

    # Set survival at lethal temperatures to zero
    Pr_surv_T[T_regime_lethal] = 0

    # Compute the survival probability within the nonlethal temperature range and clip unrealistic values
    Pr_surv_T[T_regime_nonlethal] = fn_quadFit(quadFit_c, waterTemp_mean[T_regime_nonlethal], quadFit_T0, quadFit_Tm)
    Pr_surv_T = np.clip(Pr_surv_T, 0, 1) # clip lower and upper values to be within 0 to 1

    return Pr_surv_T




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




## calc_survProb_dly_eisenQuadFit()
#  Calculates overall daily survival probability for a mosquito immature given all environmental factors and life stage.
#  Args:
#    waterTemp     - daily mean water temperature [deg C]
#    waterHeight   - daily water height [mm]
#    vpd           - daily vapor pressure deficit [kPa]
#    sunExp        - sun exposure fraction (constant for a given container)
#    lifeStage     - life stage of immatures (options: "Eggs", "Larvae", "Pupae")
#  Output:
#    Pr_surv       - daily survival probability
def calc_survProb_dly_eisenQuadFit(waterTemp, waterHeight, vpd, sunExp, lifeStage):

    Pr_surv_nominal = 0.99   # default daily Pr(survival) given no other influences
    Pr_surv_temp    = calc_survProb_temp_dly_eisenQuadFit(waterTemp, lifeStage) # daily Pr(survival) due to temperature
    Pr_surv_desicc  = calc_survProb_desicc_dly(waterHeight, vpd, sunExp, lifeStage)   # daily Pr(survival) due to desiccation

    Pr_surv = Pr_surv_nominal * Pr_surv_temp * Pr_surv_desicc
    return Pr_surv








if False:
    df_survVsTemp_daily_larvae = survCurve_Eisen_convert2dly("Larvae", False)


if False:
    for lifeStage in ["Eggs", "Larvae", "Pupae"]:
        survCurve_Eisen_convert2dly(lifeStage, True)


if False:
    for lifeStage in ["Eggs", "Larvae", "Pupae"]:
        plt_compare_survCurves_eisenQuadFit(lifeStage)
        plt_compare_survCurves_eisenBriereFit(lifeStage)

        
if False:
    for lifeStage in ["Eggs", "Larvae", "Pupae"]:
        plt_compare_survCurves_eisenQuadFit_onlyClip(lifeStage)