# mosquitoSimCity_4cast
Files associated with the AGU GeoHealth journal article "Can Seasonal Meteorological Forecasts Support Mechanistic Predictions of Aedes aegypti Abundance? Evaluating GEOS-S2S Hindcasts for Climate-Informed Dengue Early Warning".

This repository primarily contains code for data analysis and plotting, organized into the following folders:
- 0_constants: Defines directories and filepaths used throughout the code. (NOTE: take care to update these before running the code yourself)
- 1_metInputs: Prepares meteorological inputs (from MERRA-2/IMERG or GEOS-S2S) for Step 1 of the modeling pipeline.
- 2_container_modeling: Runs Step 1 of the modeling pipeline--container habitat dynamics modeling using WHATCH'EM (this is a modified version of the WHATCH'EM code).
- 3_simulateVectorSurvival: Runs Step 2 and 3 of the modeling pipeline--vector biology and vector population dynamics simulations.
- 4_data_analysis: Analyzes data output by the modeling pipeline.
- 5_plotting: Creates plots of data.

The repository also includes the meteorological inputs for the modeling pipeline, found in 000_data/input_data/. These files were produced by the code in 1_metInputs/ and are the inputs for the code in 2_container_modeling/.
