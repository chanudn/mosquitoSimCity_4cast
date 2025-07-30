# scripts_createMetInput.py
#
# Assorted blocks of code for creating meteorological input files for WHATCH'EM.
# Essentially just runs the create_met_input() function from createMetInput.py
# with different arguments.


# import required packages
# ************************************************
import calendar
import sys
import os

# import constants
# ************************************************
sys.path.append('/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/0_constants')
from CONST_FILE_DIR import (DIR_WHATCHEM_INPUT_DATA_BASE)
sys.path.append("/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/1_metInputs")
from createMetInput import create_met_input

# Define parameters
dataSrcs   = ['GEOS'] #['MERRA_IMERG', 'GEOS']
lead_times = [0]           # only relevant for GEOS
ens_nums   = range(1, 5+1) # only relevant for GEOS
locs       = ['Negombo', 'Jaffna', 'NuwaraEliya']
years      = range(2001, 2020+1)
months     = range(1, 12+1)

## Create met inputs for each month of our time range of interest (2001 to 2020)
if True:
    for dataSrc in dataSrcs:
        print(f'{dataSrc}')
        fdir_out = os.path.join(DIR_WHATCHEM_INPUT_DATA_BASE, dataSrc) # define output directory based on data source
        for loc in locs:
            print(f'{loc}')
            for year in years:
                print(f' {year}')
                for month in months:
                    print(f'  {month}')

                    # Define start and end times in YYYYMMDDHH format
                    days_in_month = calendar.monthrange(year, month)[1]
                    t_start = year * 1000000 + month * 10000 + 100
                    t_end = year * 1000000 + month * 10000 + days_in_month * 100 + 23


                    # Next, run the code to create the met input files
                    if dataSrc == 'GEOS':
                        for lead_time in lead_times:
                            print(f'   Lead time: {lead_time}')
                            for ens_num in ens_nums:
                                print(f'    ens{ens_num}')
                                create_met_input(dataSrc, loc, t_start, t_end, fdir_out, lead_time, ens_num)
                    else: # if MERRA_IMERG
                        create_met_input(dataSrc, loc, t_start, t_end, fdir_out)

