# plot_survDev_eqns.py
#
# Create plotted examples of the survival and development equations used in simVectorSurvival.py.


# import required packages
# ************************************************
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
from matplotlib.lines import Line2D

# import functions
# ************************************************
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/3_simulateVectorSurvival')  # directory of plotting files
from simVectorSurvival import calc_devRate, calc_survProb_temp_dly, calc_survProb_desicc_dly


lifeStages = ['Eggs', 'Larvae', 'Pupae']
lifeStage_properties = {
                        'Eggs':         {'color': '#f394c5', 'col_abbr': 'E'}, # color for vector biology (eggs) variables
                        'Larvae':       {'color': '#b9216e', 'col_abbr': 'L'}, # color for vector biology (larvae) variables
                        'Pupae':        {'color': '#450c29', 'col_abbr': 'P'}, # color for vector biology (pupae) variables
                        'Larvae/Pupae': {'color': '#7f174c', 'col_abbr': 'L'}  # midpoint of larvae and pupae colors
                       }

# Development rate
xmin_dev  = 0
xmax_dev  = 50
xstep_dev = 0.25
waterTemp = pd.Series(np.arange(xmin_dev, xmax_dev+xstep_dev, xstep_dev))
df_devRate = pd.DataFrame({'TW (C)': waterTemp,
                           'dev_E': calc_devRate(waterTemp, lifeStages[0]),
                           'dev_L': calc_devRate(waterTemp, lifeStages[1]),
                           'dev_P': calc_devRate(waterTemp, lifeStages[2])})
print(df_devRate)


# Survival rate (temperature-dependent)
xmin_surv_temp  = -15
xmax_surv_temp  = 50
xstep_surv_temp = 0.25
dummyTemp = 25 # temp that guarantees 100% survival (where we connect the min and max surv curves)
waterTemp_min = pd.Series(np.arange(xmin_surv_temp, dummyTemp+xstep_surv_temp, xstep_surv_temp))
waterTemp_max = pd.Series(np.arange(dummyTemp+xstep_surv_temp, xmax_surv_temp+xstep_surv_temp, xstep_surv_temp))
dummyTemps_4min = pd.Series([dummyTemp] * len(waterTemp_min)) 
dummyTemps_4max = pd.Series([dummyTemp] * len(waterTemp_max)) 
surv_E_temp_min = pd.Series(calc_survProb_temp_dly(waterTemp_min, dummyTemps_4min, lifeStages[0]))
surv_E_temp_max = pd.Series(calc_survProb_temp_dly(dummyTemps_4max, waterTemp_max, lifeStages[0]))
surv_E_temp = pd.concat([surv_E_temp_min, surv_E_temp_max])
surv_L_temp_min = pd.Series(calc_survProb_temp_dly(waterTemp_min, dummyTemps_4min, lifeStages[1]))
surv_L_temp_max = pd.Series(calc_survProb_temp_dly(dummyTemps_4max, waterTemp_max, lifeStages[1]))
surv_L_temp = pd.concat([surv_L_temp_min, surv_L_temp_max])
surv_P_temp_min = pd.Series(calc_survProb_temp_dly(waterTemp_min, dummyTemps_4min, lifeStages[2]))
surv_P_temp_max = pd.Series(calc_survProb_temp_dly(dummyTemps_4max, waterTemp_max, lifeStages[2]))
surv_P_temp = pd.concat([surv_P_temp_min, surv_P_temp_max])

df_surv_temp = pd.DataFrame({'TW_extreme (C)': pd.concat([waterTemp_min, waterTemp_max]),
                             'surv_E_temp': surv_E_temp,
                             'surv_L_temp': surv_L_temp,
                             'surv_P_temp': surv_P_temp})
print(df_surv_temp)


# Survival rate (desiccation-dependent)
waterHeight = 0 # [mm] just a value corresponding to 'dried out'
sunExp = 0.5 # same as our model
xmin_surv_desicc  = 0
xmax_surv_desicc  = 3.5
xstep_surv_desicc = 0.025
vpd = pd.Series(np.arange(xmin_surv_desicc, xmax_surv_desicc+xstep_surv_desicc, xstep_surv_desicc)) # [kPa]

surv_E_desicc = calc_survProb_desicc_dly(waterHeight, vpd, sunExp, lifeStages[0])
surv_L_desicc = calc_survProb_desicc_dly(waterHeight, vpd, sunExp, lifeStages[1])
surv_P_desicc = calc_survProb_desicc_dly(waterHeight, vpd, sunExp, lifeStages[2])

df_surv_desicc = pd.DataFrame({'VPD (kPa)': vpd,
                               'surv_E_desicc': surv_E_desicc,
                               'surv_L_desicc': surv_L_desicc,
                               'surv_P_desicc': surv_P_desicc})
print(df_surv_desicc)



## Plotting
dir_out = '/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/5_plotting/' # TEMPORARY!


# Plot!
if True:
    fig, ax0 = plt.subplots(figsize=(16, 12))

    # Plot development rate
    for lifeStage in lifeStages:
        colName = 'dev_' + lifeStage_properties[lifeStage]['col_abbr']
        color = lifeStage_properties[lifeStage]['color']
        ax0.plot(df_devRate['TW (C)'], df_devRate[colName], label=lifeStage, marker='', lw=4, color=color, alpha=0.8)
    ax0.set_xlabel('Water temperature (°C)', fontsize=18)
    ax0.set_ylabel('Development rate (hr$^{-1}$)', fontsize=18)

    # Manually set the position of the inset for temperature-dependent survival factor
    ax1_pos_left = 0.20
    ax1_pos_bot  = 0.57
    ax1_width    = 0.29
    ax1_height   = 0.25
    ax1 = fig.add_axes([ax1_pos_left, ax1_pos_bot, ax1_width, ax1_height])  # position/size in fractions of figure size

    # Plot temperature-dependent survival factor
    # Note: this is identical for larvae and pupae, so we represent both with one curve
    for lifeStage in ['Eggs', 'Larvae/Pupae']:
        colName = 'surv_' + lifeStage_properties[lifeStage]['col_abbr'] + '_temp'
        color = lifeStage_properties[lifeStage]['color']
        ax1.plot(df_surv_temp['TW_extreme (C)'], df_surv_temp[colName], label=lifeStage, marker='', lw=4, color=color, alpha=0.8)
    ax1.set_xlabel('Water temperature extreme (°C)', fontsize=16)
    ax1.set_ylabel('Survival factor (day$^{-1}$)', fontsize=16)

    # Manually set the position of the inset for desiccation-dependent survival factor
    ax2_pos_left = ax1_pos_left + ax1_width
    ax2_pos_bot  = ax1_pos_bot
    ax2_width    = ax1_width*0.75
    ax2_height   = ax1_height
    ax2 = fig.add_axes([ax2_pos_left, ax2_pos_bot, ax2_width, ax2_height])  # position/size in fractions of figure size

    # Plot desiccation-dependent survival factor
    for lifeStage in lifeStages:
        colName = 'surv_' + lifeStage_properties[lifeStage]['col_abbr'] + '_desicc'
        color = lifeStage_properties[lifeStage]['color']
        ax2.plot(df_surv_desicc['VPD (kPa)'], df_surv_desicc[colName], label=lifeStage, marker='', lw=4, color=color, alpha=0.8)
    ax2.set_xlabel('Vapor pressure deficit (kPa)', fontsize=16)
    
    # Adjust axes for the main plot
    ax0_x_ticks = np.round(np.arange(xmin_dev, xmax_dev + 0.01, 10),0)
    ax0.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax0.set_xticks(ax0_x_ticks)
    ax0.tick_params(axis='both', which='minor', length=7, color='gray', direction='in', width=0.8, labelsize=16, left=True, right=True, bottom=True, top=True) # adjust minor ticks for y-axis
    ax0.tick_params(axis='both', which='major', length=11, color='gray', direction='in', width=1.2, labelsize=16, left=True, right=True, bottom=True, top=True) # adjust major ticks for y-axis
    ax0.xaxis.set_minor_locator(AutoMinorLocator())
    ax0.yaxis.set_minor_locator(AutoMinorLocator())
    ax0.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Adjust axes for inset plots
    ax1_x_ticks = np.round(np.arange(xmin_surv_temp+5, xmax_surv_temp + 0.01, 10),0)
    ax2_x_ticks = np.round(np.arange(xmin_surv_desicc, xmax_surv_desicc + 0.01, 1),1)
    ax_x_ticks_list = [ax1_x_ticks, ax2_x_ticks]
    for i, ax in enumerate([ax1, ax2]):
        ax.set_xticks(ax_x_ticks_list[i])
        ax.tick_params(axis='both', which='minor', length=4, color='gray', direction='in', width=0.5, labelsize=14, left=True, right=True, bottom=True, top=True) # adjust minor ticks for y-axis
        ax.tick_params(axis='both', which='major', length=7, color='gray', direction='in', width=0.7, labelsize=14, left=True, right=True, bottom=True, top=True) # adjust major ticks for y-axis
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax1_y_ticks = np.round(np.arange(0, 1+0.01, 0.2),1)
    ax1.set_yticks(ax1_y_ticks)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax2.tick_params(axis='y', which='both', labelleft=False, labelright=False)  # Hide left ticks

    fileName = 'plt_survDevEqns'
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)



# Backup
# This one plots the three as separate plots in one panel
if False:
    fig, ax = plt.subplots(2, 2, figsize=(16, 16), gridspec_kw={'height_ratios': [2, 1]})

    # Merge the top-left and top-right subplots into one plot
    gs = ax[0, 0].get_gridspec()  # Get the gridspec from ax
    for a in ax[0, 0:]:  # Remove all ax in the first row
        a.remove()
    ax0 = fig.add_subplot(gs[0, :]) # create axis spanning both columns for the top plot
    ax1 = ax[1,0]
    ax2 = ax[1,1]

    # Plot development rate
    for lifeStage in lifeStages:
        colName = 'dev_' + lifeStage_properties[lifeStage]['col_abbr']
        color = lifeStage_properties[lifeStage]['color']
        df_devRate.plot(x='TW (C)', y=colName, label=lifeStage, marker='', lw=4, color=color, alpha=0.8, ax=ax0)

    # Plot temperature-dependent survival factor
    # Note: this is identical for larvae and pupae, so we represent both with one curve
    for lifeStage in ['Eggs', 'Larvae/Pupae']:
        colName = 'surv_' + lifeStage_properties[lifeStage]['col_abbr'] + '_temp'
        color = lifeStage_properties[lifeStage]['color'],
        df_surv_temp.plot(x='TW_extreme (C)', y=colName, label=lifeStage, marker='', lw=4, color=color, alpha=0.8, ax=ax1)

    # Plot desiccation-dependent survival factor
    for lifeStage in lifeStages:
        colName = 'surv_' + lifeStage_properties[lifeStage]['col_abbr'] + '_desicc'
        color = lifeStage_properties[lifeStage]['color']
        df_surv_desicc.plot(x='VPD (kPa)', y=colName, label=lifeStage, marker='', lw=4, color=color, alpha=0.8, ax=ax2)

    # Adjust axes for all the plots
    ax0.set_xlabel('Water temperature (°C)', fontsize=16)
    ax0.set_ylabel('Development rate (hr$^{-1}$)', fontsize=16)
    ax1.set_xlabel('Water temperature extreme (°C)', fontsize=16)
    ax1.set_ylabel('Survival factor', fontsize=16)
    ax2.set_xlabel('Vapor pressure deficit (kPa)', fontsize=16)
    
    # Adjust axes for the two survival factor plots
    for ax in [ax1, ax2]:
        ax.tick_params(axis='y', which='minor', length=4, color='gray', direction='in', width=0.5, labelsize=12, left=True, right=True) # adjust minor ticks for y-axis
        ax.tick_params(axis='y', which='major', length=7, color='gray', direction='in', width=0.7, labelsize=12, left=True, right=True) # adjust major ticks for y-axis
    ax2.tick_params(axis='y', which='both', labelleft=False, labelright=False)  # Hide left ticks

    # Adjust the subplot parameters to overlap
    ax0.set_position([0.05, 0.55, 0.90, 0.40])  # Adjust the position for ax0 to cover more area
    ax1.set_position([0.05, 0.10, 0.45, 0.40])  # Bottom left
    ax2.set_position([0.50, 0.10, 0.45, 0.40])  # Bottom right

    fileName = 'plt_survDevEqns'
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)
    plt.savefig(dir_out + fileName + '.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

