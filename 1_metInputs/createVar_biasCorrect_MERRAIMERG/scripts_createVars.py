# Script to run variable creating functions.

# import required packages
# ************************************************
import xarray as xr
import sys

# import required functions
# ************************************************
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/1_metInputs")
from createVar_MERRA_RH import createVar_MERRA_RH


# create RH variable from MERRA T2M, QV2M, PS. (No bias correction of T2M)

# define directories and filepaths
fDir_MERRA = '/mnt/redwood/local_drive/chanud/RawData/MERRA_1980-2021_SL/'
fPath_in_MERRA_QV2M = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.QV2M.nc'
fPath_in_MERRA_T2M  = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.T2M.nc'
fPath_in_MERRA_PS   = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.PS.nc'
fPath_out_MERRA_RH  = fDir_MERRA + 'MERRA_hourly_SLV/MERRA2_tavg1_2d_slv_Nx.19800101_20211231.SUB.RH.nc'

## Import and extract MERRA data
# Load the MERRA reanalysis data as xarray datasets
QV2M_xrds = xr.open_dataset(fPath_in_MERRA_QV2M)  # 2-m specific humidity [kg/kg]
T2M_xrds  = xr.open_dataset(fPath_in_MERRA_T2M)   # 2-m temperature [K]
PS_xrds   = xr.open_dataset(fPath_in_MERRA_PS)    # surface pressure [Pa]

# Convert to data arrays
QV2M_xrda = QV2M_xrds['QV2M']
T2M_xrda  = T2M_xrds['T2M']
PS_xrda   = PS_xrds['PS']

# Close the original datasets.
QV2M_xrds.close()
T2M_xrds.close()
PS_xrds.close()

# run function
createVar_MERRA_RH(QV2M_xrda, T2M_xrda, PS_xrda, fPath_out_MERRA_RH)