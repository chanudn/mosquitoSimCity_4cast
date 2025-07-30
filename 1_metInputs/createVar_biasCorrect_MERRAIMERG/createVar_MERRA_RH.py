# createVar_MERRA_RH.py
#
# WHATCH'EM requires relative humidity, RH, for calculating vapor pressure.
# However our MERRA-2 satellite data is specific humidity.
# So we calculate RH from MERRA-2...
# - 2-m specific humidity (QV2M)
# - 2-m air temperature (T2M)
# - surface pressure (PS)

import xarray as xr
import numpy as np

def createVar_MERRA_RH(q_xrda, T_xrda, p_xrda, fpath_out_RH):

    ## Calculate RH
    #  ***********************************************************************************
    # We will be calculating RH as
    #   RH = 100*(w/w_s)
    # where
    #   w   - mass mixing ratio of water vapor to dry air [dimensionless]
    #   w_s - saturation mass mixing ratio of water vapor to dry air [dimensionless]

    # First, define constants for RH calculation
    e_s0 = 611      # [Pa] at T_0 = 0C
    L_v  = 2.442e6  # [J/kg] at 25C
    R_v  = 461.5    # [J/kg/K]
    T_0  = 273.15   # [K]
    R_d  = 287.1    # [J/kg/K]

    # Next, we calculate w.
    # We can do so from specific humidity, q (i.e. the mass mixing
    # ratio of water vapor to TOTAL air):
    #      q = w / (w+1)
    #   -> w = q / (1-q)
    w = q_xrda / (1 - q_xrda)  # [dimensionless]

    # Now we calculate w_s.
    # But before we do so we need to calculate saturation vapor pressure, e_s.
    #   e_s = e_s0*exp[ (L_v/R_v)*( (1/T_0)-(1/T) ) ]
    # where
    #   e_s0 - saturation vapor pressure at T_0 [Pa]
    #   L_v  - specific enthalpy of vaporization [J/kg]
    #   R_v  - specific gas constant for water vapor [J/kg/K]
    #   T_0  - reference temperature [K] (we use 273.15 K)
    #   T    - temperature [K]
    #
    # Note: We get L_v at 25C from https://www.engineeringtoolbox.com/water-properties-d_1573.html.
    #       It's a fn of T (2.500e6 at 0C, 2.442e6 at 25C, 2.256e6 at 100C), but we use a constant value.
    e_s = e_s0 * np.exp((L_v / R_v) * ((1 / T_0) - (1 / T_xrda)))  # [Pa]

    # Now we can calcuate w_s as
    #   w_s = (e_s*R_d)/(R_v*(p-e_s))
    # where
    #   e_s  - saturation vapor pressure [Pa]
    #   R_d  - specific gas constant for dry air [J/kg/K]
    #   R_v  - specific gas constant for water vapor [J/kg/K]
    #   p    - pressure [Pa]
    w_s = (e_s * R_d) / (R_v * (p_xrda - e_s))  # [dimensionless]

    # Finally we can calculate relative humidity RH
    RH_xrda = 100 * (w / w_s)  # [%]

    # Create a new dataset for RH
    RH_xrda.name = 'RH'
    RH_xrda.attrs.update(q_xrda.attrs) # copy over attributes, then overwrite some of them
    RH_xrda.attrs['units'] = '%'
    RH_xrda.attrs['long_name'] = 'relative_humidity'
    RH_xrda.attrs['standard_name'] = 'relative_humidity'
    RH_xrds = RH_xrda.to_dataset()

    # Save data as NetCDF file
    RH_xrds.to_netcdf(fpath_out_RH)




