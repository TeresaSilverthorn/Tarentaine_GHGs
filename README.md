# Tarentaine_GHGs
R script and data files associated with analysis and visualization of data from the manuscript "Small water impoundments have longitudinal and local impacts on carbon cycling along a river network"

This repository includes: 
- an R scirpt used to analyze and visualize the data (Data_analysis_Tarentaine.R)
- an R script used to calculate the CO2 and CH4 fluxes (GHG_Tarentaine.R)
- an R scirpt used to calculate the leaf litter decopmosition rate (Leaf_litter_decomp.R)
- an R script used to calculate the organic matter stocks and fluxes (OM_stock_flux.R)
- an R script used to calculate the air and water temperature from continuous ibutton datalogger measurements (iButton_air_temp.R)
- a .csv file file with the CO2 and CH4 fluxes (OLD VERSION: CO2.CH4.fluxes.csv; NEW VERSION includes 4 previously missing values: CO2.CH4.fluxes.2026.csv): use the CO2 and CH4 flux values from this file as they are the most up to date and calculated based on site specific pressure values based on elevation and temperature
- a .csv file with all of the data (dat.csv): these are the flux values used in the article, they are based on a pressure value of 1 atm for all sites, and 4 values are missing
- a .csv file with the site mean data (dat_means.csv)

Data dictionary

ID unique: unique ID for each sampling point at each site at each sampling date

campaign: sampling campaign number (1,2, or 3)

date: date of the sampling

site: unique site ID

flow_state: state of flow (flowing or standing water)

transect: sampling transect at each site, 1 starts at the most downstream portion of the sampling reach

start_time: start time of the GHG measurement

end_time: end_time of the GHG measurement

water_temp_C: water temperature point measurement at the time of sampling

water_pH: water pH point measurement at the time of sampling

water_conductivity_us_cm: water conductivity point measurement at the time of sampling

DO_mg_L: dissolved oxygen point measurement at the time of sampling in mg/L

DO_.: dissolved oxygen point measurement at the time of sampling in %

canopy_cover_.: visual estimation of the canopy cover (%)

wetted_width_m: wetted width of the stream in m

discharge_m3_s: calculated discharge in m3/s for the transect

mean_velocity_m_s: mean velocity across the transect

mean_depth_cm: mean depth across the transect

LML_coarse: leaf mass loss in the coarse mesh bags

LML_fine: leaf mass loss in the fine mesh bags

k_day_coarse: leaf litter decomposition rate in coarse mesh bags

k_day_coarse: leaf litter decomposition rate in fine mesh bags

k_dday_coarse: leaf litter decomposition rate per degree day in coarse mesh bags

k_dday_fine: leaf litter decomposition rate per degree day in fine mesh bags

OM_stock_g_m2: riverine organic matter stock average in g/m2

OM_flux_g_m2_s: riverine organic matter flux average in g/m2/s

sed_OM_percent: percent of orgnic matter in a composite sediment sample

CO2_C_mg_m2_h: CO2 flux rate in mg/m2/h

CH4_C_mg_m2_h: CH4 flux rate in mg/m2/h

temp_C: average daily air temperature from the ibuttons

DailyMeanWaterTemp_C: average daily water temperature from the ibuttons

lat: latitude of the study site

lon: longitude of the study site

masl: elevation above sea level in m at the study site

Distance_to_source_km: distance to the source of the stream in km

dist_ds_dam_km: distance downstream from the nearest dam in km

bedrock: percent bedrock substrate

boulder: percent boulder substrate

cobble: percent cobble substrate

pebble_gravel: percent pebble and/or gravel substrate

fine_substrate: percent fine substrate

season: season of sampling (spring, summer, or fall) 

position: position upstream or downstream of the two largest dams in the catchment

position_d: position upstream or downstream, or at the two largest dams in the catchment






















