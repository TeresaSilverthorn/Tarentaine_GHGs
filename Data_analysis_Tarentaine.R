#### Data analysis and visualization of the OM and GHG dynamics in a river network (Tarentaine) impacted by damming ####

#Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(purrr)
library(sf)
library(randomForest)
library(caret) #for making predictions with RF
library(caTools)
library(ggpubr)
library(olsrr) #for VIF and tolerance
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(cAIC4) #step AIC
library(FactoMineR) #for PCA
library(factoextra) #for visualizing CPA
library(buildmer) #for model selection
library(data.table)
library(vegan)
library(gam)
library(mgcv) #gamm
library(itsadug) #plotting gam
library(RVAideMemoire) #for pairwise comparison of permanova
library(geosphere)
library(pls)
library(emmeans)
library(scales)

# Set wd for figures

setwd("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Tarentaine_GHGs/Figures")

#Load and merge the 1. ancillary data, 2. GHG flux data, 3. OM stock and flux data, 4. leaf litter decomposition data, 5. daily temperature data, 6. Sediment OM content. 7. Latitude and Longitude. 8. Import network distances. Check for quality control 

#### 1. Load ancillary data ####

ancil_dat <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Ancillary data/Data_entry_Tarentaine_2022_2023-08-14.csv", header=T) #Make sure it's the latest version from the Google Drive

str(ancil_dat) #327 obs of 42 vars

#Drop the useless columns, get the depth and velocity average per transect
ancil_dat <- ancil_dat %>%
  select(-CO2_start_offset_mins, -CO2_end_offset_mins, -CH4_start_offset_mins, -CH4_end_offset_mins, -ibutton_air, -ibutton_water, -X, -chamber_area, -chamber_volume, -OM_flux_net_time_mins) %>%
  rowwise() %>%
  select(-starts_with("depth"), -starts_with("velocity"))

#Add a "ID_unique" column with this format : "2022-05-30_TA01_1"

#rename 
ancil_dat <- ancil_dat %>%
  rename(bedrock=X._bedrock, 
         boulder=X._boulder, 
         cobble=X._cobble, 
         pebble_gravel=X._pebble_gravel, 
         fine_substrate=X._fine_substrate)

############################################################################
#### 2. Load in the GHG flux data ####

GHG <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Tarentaine_GHGs/CO2.CH4.fluxes.csv")

str(GHG) #323 obs. of  51 variables

#Subset just the useful columns

GHG <- GHG %>% 
select(ID_unique, CO2_C_mg_m2_h, CH4_C_mg_m2_h)


############################################################################

#### 3. Load in the OM flux and stock data ####

OM <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/OM_flux_stock/OM_stock_flux.csv")

#Drop the useless columns
OM <- OM %>%
  select(-X)

##########################################################################

#### 4. Load in leaf litter decomposition data ####

decomp <- read.csv( "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Leaf_packs/Leaf_pack_decomposition.csv")

str(decomp) #115 obs. of  5 variables

#reshape the data so coarse and fine mesh LML become their own columns per site/campaign
decomp1 <- decomp %>%
  pivot_wider(names_from = Mesh_size, values_from = LML) %>%
  rename(LML_coarse = Large, LML_fine = Small) %>%
  group_by(campaign, site) %>%   #to solve for the NA in alternating rows, group
  summarize(LML_coarse = mean(LML_coarse, na.rm = TRUE),   #then take the mean
            LML_fine = mean(LML_fine, na.rm = TRUE)) %>%
  mutate(across(everything(), ~replace(., is.nan(.), NA)))    #replace NaN with NA

str(decomp1) #58 obs of 4 vars

#do the same for k_day and kk_day
decomp2 <- decomp %>%
  pivot_wider(names_from = Mesh_size, values_from = k_day) %>%
  rename(k_day_coarse = Large, k_day_fine = Small) %>%
  group_by(campaign, site) %>%   #to solve for the NA in alternating rows, group
  summarize(k_day_coarse = mean(k_day_coarse, na.rm = TRUE),   #then take the mean
            k_day_fine = mean(k_day_fine, na.rm = TRUE)) %>%
  mutate(across(everything(), ~replace(., is.nan(.), NA)))    #replace NaN with NA


decomp3 <- decomp %>%
  pivot_wider(names_from = Mesh_size, values_from = k_dday) %>%
  rename(k_dday_coarse = Large, k_dday_fine = Small) %>%
  group_by(campaign, site) %>%   #to solve for the NA in alternating rows, group
  summarize(k_dday_coarse = mean(k_dday_coarse, na.rm = TRUE),   #then take the mean
            k_dday_fine = mean(k_dday_fine, na.rm = TRUE)) %>%
  mutate(across(everything(), ~replace(., is.nan(.), NA)))    #replace NaN with NA

#merge decomp1, decomp2 and decomp 3 by campaign and site

decomp <- inner_join(decomp1, decomp2, by = c("campaign", "site")) %>%
  inner_join(decomp3, by = c("campaign", "site"))

decomp$k_ratio_f_c <- decomp$k_dday_fine/decomp$k_dday_coarse

########################################################################

#### 5. Import the daily temperature data ####

# We only have site specific daily air temperature for the first campaign. 
#we have water temperature for the leaf pack incubation period 


air_temp <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/ibutton_air_temp.csv") 

str(air_temp) #327 obs. of  3 variables


water_temp <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Daily_mean_water_temp.csv")

str(water_temp) #2824 obs. of  4 variables


###################################################################################

#### 6. Sediment OM content ####

sed_OM <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Sediment/Sediment_OM_Tarentaine.csv") 

#subset only the useful columns
sed_OM <- sed_OM %>%
  select(sed_OM_percent, Site, Campaign)

#Remove the weird outlying value of sed_OM_percent "14.5095801"

sed_OM <- sed_OM %>%
  filter(sed_OM_percent <= 14)

#Rename columns to match
sed_OM <- sed_OM %>%
  rename(site = Site,
         campaign = Campaign)


############################################################################
#### 7. add latitude and longitudem and altitude ####

coords <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Site selection/GPS_coordinates_site_selection_2022-08-14.csv")

# subset useful columns

coords <- coords %>%
  select(Site, Lat, Lon, masl)  %>%
  rename(site=Site, 
         lat=Lat, 
         lon=Lon)
  
############################################################################

#### 8. Add network distances #####

dist <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Network distances/Tarentaine_distance_to_source.csv")

str(dist)

dist <- dist  %>%
  rename(site=Site, , 
         dist_ds_dam_km = "Distance_downstream_dam_5m._km", 
         dist_to_dam_km = "Distance_to_dam_5m.", 
         dist_to_source_km = "Distance_to_source_km", 
         dist_ds_weir_km = "Distance_downstream_weir_or_dam_km", 
         dist_to_weir_km = "Distance_to_weir_or_dam_km", 
         percent_open_upstream_weir = "percent_open_upstream_weir_dam")

############################################################################

#### Merge the dataframes together ####

# Join multiple dataframes
df_list = list(ancil_dat, decomp, OM, sed_OM)  #join the ones that have site and campaign first
dat <- df_list %>%
  reduce(full_join, by = c("campaign", "site")) 

#Add an ID_unique column to dat

dat <- dat  %>%
  mutate(ID_unique = paste(date, site, transect, sep = "_"))%>%
  select(ID_unique, everything())

#Now join the ones with ID_unique
df_list = list(dat, GHG, air_temp)  
dat <- df_list %>%
  reduce(full_join, by = c("ID_unique")) 

#Now join water temp by date and site

water_temp <- water_temp %>%
  rename(date = Date)  #rename

df_list = list(dat, water_temp) 
dat <- df_list %>%
  reduce(left_join, by = c("site", "date"))  #use a left join here

#merge coords just by site
dat <- dat %>%
  left_join(coords, by = "site")

#merge dist by site
dat <- dat %>%
  left_join(dist, by = "site")


#remove unnecessary columns

str(dat)

dat <- dat %>%
  select(-X.x, -X.y, -Notes) 

#########################################################################

#### Fill any missing data ####

## add the substrate type to all of the campaigns

columns_to_fill <- c("bedrock", "boulder", "cobble", "pebble_gravel", "fine_substrate")

# Fill missing values in specified columns based on site
dat <- dat %>%
  group_by(site) %>%
  fill(!!!syms(columns_to_fill), .direction = "downup") %>%
  ungroup()

#Need to fill in values for TA02 and TA06, based on photos
#TA02: 70% boulder, 20% cobble, 10% gravel
#TA06: 15% boulders, 65% cobble, 15% gravel, 5% fine.

fill_values_TA02 <- data.frame(
  site = "TA02",
  bedrock= 0,
  boulder = 70,
  cobble = 20,
  pebble_gravel = 10, 
  fine_substrate = 0 )

fill_values_TA06 <- data.frame(
  site = "TA06",
  bedrock = 0,
  boulder = 15,
  cobble = 65,
  pebble_gravel = 15,
  fine_substrate = 5 )

# Combine fill values for both sites
all_fill_values <- bind_rows(fill_values_TA02, fill_values_TA06)

# Replace "TA02" and "TA06" rows with fill values
dat <- dat %>%
  left_join(all_fill_values, by = "site") %>%
  mutate(
    bedrock = if_else(is.na(bedrock.y), bedrock.x, bedrock.y),
    boulder = if_else(is.na(boulder.y), boulder.x, boulder.y),
    cobble = if_else(is.na(cobble.y), cobble.x, cobble.y),
    pebble_gravel = if_else(is.na(pebble_gravel.y), pebble_gravel.x, pebble_gravel.y),
    fine_substrate = if_else(is.na(fine_substrate.y), fine_substrate.x, fine_substrate.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>%
  arrange(campaign, site, transect)

###############################################################################

### Many statistical analyses don't allow NA in the data, so we need to fill using data imputation with means

# Print the rows with NA values in certain columns
print(dat[(which(is.na(dat$water_temp_C))), ]) # 3 rows
print(dat[(which(is.na(dat$water_pH))), ]) # 5 rows
#canopy cover had 9 missing values which I filled with the next closest transect in the excel file
#16 missing values for wetted width replaced with average from the other 2 campaigns, or other transects at the same site, noted on excel 
print(dat[(which(is.na(dat$discharge_m3_s))), ]) #16 values
#Fill NA's in discharge with 0 for the reservoirs (TA05 and TA10), the rest fill in from other transects
dat <- dat %>%
  mutate(discharge_m3_s = case_when(
    site %in% c("TA10", "TA05") & is.na(discharge_m3_s) ~ 0,
    TRUE ~ discharge_m3_s
  ))

print(dat[(which(is.na(dat$mean_depth_cm))), ]) #10 cases for TA10 and TA05 in C3, so use the mean of C1 and C2, in the google sheet

print(dat[(which(is.na(dat$LML_coarse))), ])    #5 cases all at TA05 replace with mean 

dat <- dat %>%
  mutate(LML_coarse = ifelse(is.na(LML_coarse) & site == "TA05" & campaign == "3", mean(LML_coarse,na.rm=T), LML_coarse)) %>%
  mutate(k_day_coarse = ifelse(is.na(k_day_coarse) & site == "TA05" & campaign == "3", mean(k_day_coarse,na.rm=T), k_day_coarse))%>%
  mutate(k_dday_coarse = ifelse(is.na(k_dday_coarse) & site == "TA05" & campaign == "3", mean(k_dday_coarse,na.rm=T), k_dday_coarse))

print(dat[(which(is.na(dat$OM_stock_g_m2))), ]) #3 NAs, replace with mean

dat <- dat %>%
  mutate(OM_stock_g_m2 = ifelse(is.na(OM_stock_g_m2) & site == "TA11" & campaign == "3", mean(OM_stock_g_m2,na.rm=T), OM_stock_g_m2))

print(dat[(which(is.na(dat$OM_flux_g_m2_s))), ]) #3 NAs, replace with mean

dat <- dat %>%
  mutate(OM_flux_g_m2_s = ifelse(is.na(OM_flux_g_m2_s) & site == "TA11" & campaign == "3", mean(OM_flux_g_m2_s,na.rm=T), OM_flux_g_m2_s))

print(dat[(which(is.na(dat$sed_OM_percent))), ]) #8 NAs, replace with mean

dat <- dat %>%
  mutate(sed_OM_percent = ifelse(is.na(sed_OM_percent) & campaign == "3", mean(sed_OM_percent,na.rm=T), sed_OM_percent))

#replace with the mean
dat$water_temp_C <- ifelse(is.na(dat$water_temp_C), mean(dat$water_temp_C, na.rm=T), dat$water_temp_C)
dat$water_pH <- ifelse(is.na(dat$water_pH), mean(dat$water_pH, na.rm=T), dat$water_pH)

print(dat[(which(is.na(dat$CO2_C_mg_m2_h))), ]) #4 NAs, for C1 TA10 replace with the mean of transect 1-5
#FOR TA14 in C3, replace the NA in T2 and T5 with the mean from 1,3, and 4

mean_CO2_TA10 <- mean(dat %>% filter(campaign == "1" & site == "TA10" & transect %in% 1:5) %>% pull(CO2_C_mg_m2_h), na.rm = TRUE)
mean_CO2_TA14 <- mean(dat %>% filter(campaign == "3" & site == "TA14" & transect %in% c(1, 3, 4))  %>% pull(CO2_C_mg_m2_h), na.rm = TRUE)
mean_CO2_TA01 <- mean(dat %>% filter(campaign == "3" & site == "TA01" & transect %in% c(1,2,3,5,6)) %>% pull(CO2_C_mg_m2_h), na.rm = TRUE)

dat$CO2_C_mg_m2_h <- ifelse(is.na(dat$CO2_C_mg_m2_h) & dat$campaign == "1" & dat$site == "TA10" & dat$transect == 6, mean_CO2_TA10, dat$CO2_C_mg_m2_h)
dat$CO2_C_mg_m2_h <- ifelse(is.na(dat$CO2_C_mg_m2_h) & dat$campaign == "3" & dat$site == "TA14", mean_CO2_TA14, dat$CO2_C_mg_m2_h)
dat$CO2_C_mg_m2_h <- ifelse(is.na(dat$CO2_C_mg_m2_h) & dat$campaign == "3" & dat$site == "TA01", mean_CO2_TA01, dat$CO2_C_mg_m2_h)


print(dat[(which(is.na(dat$CH4_C_mg_m2_h))), ])  #4 NAs, same as CO2

mean_CH4_TA10 <- mean(dat %>% filter(campaign == "1" & site == "TA10" & transect %in% 1:5) %>% pull(CH4_C_mg_m2_h), na.rm = TRUE)
mean_CH4_TA14 <- mean(dat %>% filter(campaign == "3" & site == "TA14" & transect %in% c(1, 3, 4))  %>% pull(CH4_C_mg_m2_h), na.rm = TRUE)
mean_CH4_TA01 <- mean(dat %>% filter(campaign == "3" & site == "TA01" & transect %in% c(1,2,3,5,6)) %>% pull(CH4_C_mg_m2_h), na.rm = TRUE)

dat$CH4_C_mg_m2_h <- ifelse(is.na(dat$CH4_C_mg_m2_h) & dat$campaign == "1" & dat$site == "TA10" & dat$transect == 6, mean_CH4_TA10, dat$CH4_C_mg_m2_h)
dat$CH4_C_mg_m2_h <- ifelse(is.na(dat$CH4_C_mg_m2_h) & dat$campaign == "3" & dat$site == "TA14", mean_CH4_TA14, dat$CH4_C_mg_m2_h)
dat$CH4_C_mg_m2_h <- ifelse(is.na(dat$CH4_C_mg_m2_h) & dat$campaign == "3" & dat$site == "TA01", mean_CH4_TA01, dat$CH4_C_mg_m2_h)

###############################################################################

#Calculate CO2e
#the GWP of methane is 28 over a 100 yr horizon (Myhre et al 2013)
dat$CO2e <- dat$CO2_C_mg_m2_h + (dat$CH4_C_mg_m2_h*28)

###############################################################################

#make relevant columns as factors 
 
dat <- dat  %>%
  mutate(site = as.factor(site),
         campaign = as.factor(campaign), 
         flow_state = as.factor(flow_state), 
         transect = as.factor(transect))

# Make campaign a factor
dat$campaign <- as.factor(dat$campaign)
dat$water_temp_C <- as.numeric(dat$water_temp_C)


#make a new column where campaign is season

dat <- dat %>%
  mutate(season = case_when(
    campaign == 1 ~ "Spring",
    campaign == 2 ~ "Summer",
    campaign == 3 ~ "Fall",
    TRUE ~ "Other"
  ))

dat$season<- as.factor(dat$season)

#set spring as reference level
dat$season <- relevel(dat$season, ref = "Spring")

#Reorder factor levels
dat$season <- factor(dat$season, levels = c("Spring", "Summer", "Fall"))


#Add a variable for network position upstream/dowsntream
dat <- dat %>%
  mutate(position = case_when(
    dist_ds_dam_km < 0 ~ "upstream",
    dist_ds_dam_km > 0 ~ "downstream",
  ))

dat$position<- as.factor(dat$position)

dat$position<- relevel(dat$position, ref = "upstream")

#Add a variable for network position which also includes dam
dat <- dat %>%
  mutate(position_d = case_when(
    flow_state == "standing" ~ "reservoir",
    dist_ds_dam_km < 0 ~ "upstream",
    dist_ds_dam_km > 0 ~ "downstream",
  ))

dat$position_d<- as.factor(dat$position_d)

dat$position_d <- relevel(dat$position_d, ref = "upstream")

#Reorder factor levels
dat$position_d <- factor(dat$position_d, levels = c("upstream", "reservoir", "downstream"))

#calculate the shannon diversity of the substrate types to represent habitat diversity/complexity per:
#Barnes, J. B., Vaughan, I. P., & Ormerod, S. J. (2013). Reappraising the effects of habitat structure on river macroinvertebrates. Freshwater Biology, 58(10), 2154-2167.

dat <- dat %>%
  mutate(substrate_complexity = diversity(select(., bedrock, boulder, cobble, pebble_gravel, fine_substrate), index = "shannon"))


#### Save as .csv ####

str(dat) #327 × 54 variables

write.csv(dat, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/dat.csv")

###############################################################################

##### Calculate site averages ####

dat_means <- dat %>% 
  dplyr::select(-ID_unique, -transect, -start_time, -end_time) %>% # use select with -var_name to eliminate columns 
  dplyr::group_by(site, campaign, date, flow_state, season, position, position_d) %>% # we group by the two values
  dplyr::summarise(across(.cols = everything(), .fns = mean, .names = "{.col}"))
   
dat_means <- as.data.frame(dat_means)

str(dat_means) #58 obs of 40 vars


#calculate averages of all variables #For Table S2 in Supplementary Materials

data_summary <- dat %>%
  group_by(season) %>% #position_d , season
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  t()

mean(dat$CH4_C_mg_m2_h[dat$position_d == "downstream" & dat$site != "TA04" & dat$site != "TA09"],)
#0.1983259 0.3957998

mean(dat$CH4_C_mg_m2_h[dat$position_d == "reservoir"],)

subset_data <- dat[dat$position_d == "downstream", ]

# Printing the site column for the subset
print(subset_data$site)


write.csv(data_summary, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/dat_summary.csv")

data_summary2 <- dat %>%
  summarise(across(where(is.numeric), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE))))

data_summary_total <- data_summary %>%
  summarise(across(-campaign, ~ mean(.x)))

# Bind the total mean row to the campaign summary
data_sum <- bind_rows(data_summary, data_summary_total)

print(data_sum)

unique_sites_per_camp <- dat %>%
  select(site, campaign) %>%
  distinct() %>%
  group_by(campaign)


common_sites <- dat %>%
  group_by(site) %>%
  summarise(num_campaigns = n_distinct(campaign)) %>%
  filter(num_campaigns == n_distinct(dat$campaign)) %>%
  pull(site)

print(common_sites)
length(common_sites) #18

#summary up upstream reservoir downstream 

site_summary <- dat %>%
  filter(site %in% c("TA09", "TA10", "TA11" )) %>%  ##   "TA04", "TA05", "TA06"  #"TA09", "TA10", "TA11"  
  group_by(season, site, position_d) %>% #add season for ANOVA #otherwise site, position_d
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

write.csv(site_summary, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/site_summary.csv")


#ANOVA on upstream, reservoir, downstream differences in env vars
anova <- aov(OM_g_m3 ~ site, data = site_summary)
summary(anova)
TukeyHSD(anova)

sd(dat$OM_g_m3, na.rm=T)

sd(dat$CO2e)
###############################################################################
#### Make preliminary plots ####

#Can separate sites for plotting

#Subset just the upstream sites, TA06, TA11, TA07, TA08, TA12 and higher
upstream_sites <- c("TA06", "TA05", "TA10", "TA11", "TA07", "TA08", "TA12", "TA13", "TA14", "TA15", "TA17", "TA20", "TA21", "TA22", "TA24", "TA02R")
#15 sites
downstream_sites <- c("TA04", "TA09", "TA01", "TA02", "TA03")
#7 sites


#### Make a correlation plot ####

#subset just the numeric columns
dat_numeric <- dat %>%
  select_if(is.numeric)

dat_numeric <- na.omit(dat_numeric) #remove NAs

#remove irrelevant variables
dat_numeric <- dat_numeric   %>%
  select(-LML_coarse, -LML_fine, -k_day_fine, -k_day_coarse, -lat, -lon)

M <-cor(dat_numeric)

#corrplot(M, type="upper", order="hclust",
         #tl.cex = 0.7, tl.srt = 45,
        # col=brewer.pal(n=8, name="RdYlBu"))

#Another option with the correlation coefficients
corrplot(cor(dat_numeric), method = "number", tl.cex = 0.7, tl.srt = 45)

#Hard to see, so in table form
correlation_matrix <- cor(dat_numeric)
correlation_df <- as.data.frame(as.table(correlation_matrix))

#### Make plots for each reservoir ####

#subset the data from each reservoir
eau_verte <- dat %>%
  filter(site %in% c("TA04", "TA05", "TA06"))  %>% 
group_by(site, season, position_d) %>%   #add site if you want site means
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))

tarentaine <- dat %>%
  filter(site %in% c("TA09", "TA10", "TA11"))   %>% 
group_by(site, season, position_d) %>%   #add site if you want site means
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))


eau_verte_sum <- dat %>%
  filter(site %in% c("TA04", "TA05", "TA06")) %>%  
  group_by(site, season, position_d) %>%   #add site if you want site means
    summarise(
      CO2_sd = sd(CO2_C_mg_m2_h, na.rm = TRUE),
      CO2_C_mg_m2_h = mean(CO2_C_mg_m2_h, na.rm = TRUE)) 
  

tarentaine_sum <- dat %>%
  filter(site %in% c("TA09", "TA10", "TA11"))  %>% 
  group_by(site, season, position_d) %>%  #add site if you want site means
    summarise(
      CO2_sd = sd(CO2_C_mg_m2_h, na.rm = TRUE),
      CO2_C_mg_m2_h = mean(CO2_C_mg_m2_h, na.rm = TRUE)) 


eau_verte_sum2 <- dat %>%
  filter(site %in% c("TA04", "TA05", "TA06")) %>%  
  group_by(site, season, position_d) %>%   #add site if you want site means
  summarise(
    CH4_sd = sd(CH4_C_mg_m2_h, na.rm = TRUE),
    CH4_C_mg_m2_h = mean(CH4_C_mg_m2_h, na.rm = TRUE)) 


tarentaine_sum2 <- dat %>%
  filter(site %in% c("TA09", "TA10", "TA11"))  %>% 
  group_by(site, season, position_d) %>%  #add site if you want site means
  summarise(
    CH4_sd = sd(CH4_C_mg_m2_h, na.rm = TRUE),
    CH4_C_mg_m2_h = mean(CH4_C_mg_m2_h, na.rm = TRUE)) 

sd(dat$OM_flux_g_m2_s)


# OM per reservoir

OM_tarentaine<- ggplot(data=tarentaine, aes(position_d, y=OM_stock_g_m2)) +
  geom_point(aes(colour=position_d, shape=season), alpha=0.5, size=4) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black")) +    ylab(expression("OM stock (g m"^"-2"*")"))  + scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "a) TAR") + ylim(0,101) +
  annotate("text", x = 1.8, y = 100, label=sprintf('\u2191'), size=2.5 , colour=(values=c("#26867c"))) +
annotate("text", x = 2.0, y = 99.5, label="866", size=2.5 , colour=(values=c("#26867c")))  + guides(color = FALSE)
OM_tarentaine

OM_eauverte<- ggplot(data=eau_verte, aes(position_d, y=OM_stock_g_m2)) +
  geom_point(aes(colour=position_d, shape=season), size=4, alpha=0.5) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black")) +    ylab(expression("OM stock (g m"^"-2"*")"))  + scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) + labs(title = "b) EVE") + guides(color = FALSE) + ylim(0,100)
OM_eauverte



# k coarse per reservoir

kcoarse_tarentaine<- ggplot(data=tarentaine, aes(position_d, y=k_day_coarse)) +
  geom_point(aes(colour=position_d, shape=season), alpha=0.5, size=4) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black")) +   scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "a) TAR - total") + ylab(expression(paste(italic(k), " (", day^-1, ")"))) + guides(color = FALSE) + ylim(0,0.20)
kcoarse_tarentaine

kcoarse_eauverte<- ggplot(data=eau_verte, aes(position_d, y=k_day_coarse)) +
  geom_point(aes(colour=position_d, shape=season), alpha=0.5, size=4) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black")) +   scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "b) EVE - total") + ylab(expression(paste(italic(k), " (", day^-1, ")"))) + guides(color = FALSE) + ylim(0,0.20)
kcoarse_eauverte


# k fine per reservoir

kfine_tarentaine<- ggplot(data=tarentaine, aes(position_d, y=k_day_fine)) +
  geom_point(aes(colour=position_d, shape=season), alpha=0.5, size=4) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black")) +   scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "c) TAR -  microbes") + ylab(expression(paste(italic(k), " (", day^-1, ")"))) + guides(color = FALSE)  + ylim(0,0.04)
kfine_tarentaine

kfine_eauverte<- ggplot(data=eau_verte, aes(position_d, y=k_day_fine)) +
  geom_point(aes(colour=position_d, shape=season), alpha=0.5, size=4) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black")) +   scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "d) EVE - microbes") + ylab(expression(paste(italic(k), " (", day^-1, ")"))) + guides(color = FALSE) + ylim(0,0.04)
kfine_eauverte


# CO2 per reservoir
CO2_tarentaine<- ggplot(data=tarentaine, aes(position_d, y=CO2_C_mg_m2_h)) +
  geom_pointrange(aes(ymin = CO2_C_mg_m2_h-CO2_sd, ymax = CO2_C_mg_m2_h+CO2_sd, colour=position_d, shape=season), position = position_dodge(0.8), alpha=0.5, size=1, data = tarentaine_sum) +
  theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black"))  + scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "c) TAR")  + guides(color = FALSE)  + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) 
CO2_tarentaine

CO2_eauverte<- ggplot(data=subset(dat, site=="TA05" | site=="TA04" | site=="TA06") , aes(position_d, y=CO2_C_mg_m2_h)) +
  geom_pointrange(aes(ymin = CO2_C_mg_m2_h-CO2_sd, ymax = CO2_C_mg_m2_h+CO2_sd, colour=position_d, shape=season), position = position_dodge(0.8), alpha=0.5, size=1, data = eau_verte_sum) + theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black"))  + scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "d) EVE")   + guides(color = FALSE) + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) 
CO2_eauverte

# CH4 per reservoir
CH4_tarentaine<- ggplot(data=tarentaine, aes(position_d, y=CH4_C_mg_m2_h)) +
  geom_pointrange(aes(ymin = CH4_C_mg_m2_h-CH4_sd, ymax = CH4_C_mg_m2_h+CH4_sd, colour=position_d, shape=season), position = position_dodge(0.8), alpha=0.5, size=1, data = tarentaine_sum2) + 
  theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black"))  + scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "e) TAR")  + guides(color = FALSE) + ylab(expression(mg~CH[4]*`-C`~m^-2~h^-1))  
CH4_tarentaine

CH4_eauverte<- ggplot(data=eau_verte, aes(position_d, y=CH4_C_mg_m2_h)) +
  geom_pointrange(aes(ymin = CH4_C_mg_m2_h-CH4_sd, ymax = CH4_C_mg_m2_h+CH4_sd, colour=position_d, shape=season), position = position_dodge(0.8), alpha=0.5, size=1, data = eau_verte_sum2) +
       theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 10), axis.text = element_text(size = 10, colour="black"))  + scale_colour_manual(values=c("#317cae","#26867c" ,"#87dad1")) +labs(title = "f) EVE")  + guides(color = FALSE) + ylab(expression(mg~CH[4]*`-C`~m^-2~h^-1))  
CH4_eauverte



tiff("OM_byreservoir", units="in", width=9, height=4, res=300)
OM_byreservoir <- ggarrange(OM_tarentaine, OM_eauverte,ncol = 2, nrow = 1, align="v", common.legend = T)
OM_byreservoir
dev.off()

tiff("k_byreservoir", units="in", width=6, height=5.5, res=300)
k_byreservoir <- ggarrange(kcoarse_tarentaine, kcoarse_eauverte,
                           kfine_tarentaine, kfine_eauverte,
                           ncol = 2, nrow = 2, align="v", common.legend = T)
k_byreservoir
dev.off()


tiff("CO2_byreservoir", units="in", width=9, height=4, res=300)
CO2_byreservoir <- ggarrange(CO2_tarentaine, CO2_eauverte, ncol = 2, nrow = 1, align="v", common.legend = T)
CO2_byreservoir
dev.off()


tiff("CH4_byreservoir", units="in", width=9, height=4, res=300)
CH4_byreservoir <- ggarrange(CH4_tarentaine, CH4_eauverte,ncol = 2, nrow = 1, align="v", common.legend = T)
CH4_byreservoir
dev.off()

#combine all
tiff("vars_byreservoir.tiff", units="in", width=6, height=7, res=300)
vars_byreservoir <- ggarrange(OM_tarentaine  + rremove("x.text") , OM_eauverte + rremove("x.text"), 
                             # kcoarse_tarentaine + rremove("x.text"),  kcoarse_eauverte + rremove("x.text"), 
                              CO2_tarentaine + rremove("x.text"), CO2_eauverte + rremove("x.text") , 
                              CH4_tarentaine, CH4_eauverte,ncol = 2, nrow = 3, align="v", common.legend = T, heights= c(1, 1, 1.3)) 
vars_byreservoir
dev.off()


#### GHG by site and campaign ####
CO2_site<- ggplot(data=dat, aes(x=as.factor(site), y=CO2_C_mg_m2_h, fill=as.factor(season))) +
  geom_bar(stat="identity", position=position_dodge())+  theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) 
CO2_site

CH4_site<- ggplot(data=dat, aes(x=as.factor(site), y=CH4_C_mg_m2_h, fill=as.factor(season))) +
  geom_bar(stat="identity", position=position_dodge())+  theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) 
CH4_site

#### Get Pearson R2  ####

# List of predictor variable names
predictor_vars <- dat %>% 
  select(CH4_C_mg_m2_h, water_pH, water_conductivity_us_cm, DO_mg_L, DO_.,
         canopy_cover_., wetted_width_m, discharge_m3_s,
         mean_velocity_m_s, DailyMeanWaterTemp_C, temp_C,
         substrate_complexity, fine_substrate, mean_depth_cm,
         masl, dist_to_source_km, OM_stock_g_m2)

# Calculate Pearson correlation coefficients (R) and p-values for selected variables
correlation_results <- data.frame(
  Pearson_R = sapply(predictor_vars, function(var) cor(var, predictor_vars$CH4_C_mg_m2_h)),
  p_value = sapply(predictor_vars, function(var) cor.test(var, predictor_vars$CH4_C_mg_m2_h)$p.value))






#### GHG by position ####

tiff("CO2_position_d", units="in", width=6, height=4, res=300)
CO2_position<- ggplot(data=dat, aes(position_d, y=CO2_C_mg_m2_h, fill=season)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) + ylim(-40, 350) + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) 
CO2_position
dev.off()

tiff("CH4_position_d", units="in", width=6, height=4, res=300)
CH4_position<- ggplot(data=dat, aes(position_d, y=CH4_C_mg_m2_h, fill=season)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),   axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) + ylim(-.01, 20) + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))  
CH4_position
dev.off()

#### OM stock by position ####

OMstock_position<- ggplot(data=dat, aes(position_d, y=OM_stock_g_m2, fill=season)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +    theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black"))  
OMstock_position

#### GHG by altitude (masl) ####
CO2_masl <- ggplot(dat, aes(masl, CO2_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_masl

CH4_masl <- ggplot(dat, aes(masl, CH4_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_masl


#### GHG and distance to source ####
CO2_distance <- ggplot(subset(dat, position=="downstream"), aes(dist_to_source_km, CO2_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) + geom_smooth(method = "lm", se = F, size = 1)
CO2_distance


#relevel
dat$position <- factor(dat$position, levels = c("upstream", "downstream"))


tiff("CO2_distancebyposition", units="in", width=8.5, height=4, res=300)
CO2_distancebyposition <- ggplot(dat, aes(dist_to_source_km, CO2_C_mg_m2_h)) + geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(), legend.position = "top", legend.justification = "center", panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + scale_colour_manual(values=c("#91278e","#26867c", "#f7921e")) + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) + xlab("Distance to the source (km)") + scale_linetype_manual(values = c("dashed", "solid")) +facet_wrap(~season) +
geom_line(data = subset(dat, season=="Fall"), stat="smooth",method = "lm", colour = "black", linewidth = 1, linetype="solid", alpha=0.5)  + ylim(-35, 600) 
CO2_distancebyposition 
dev.off()

CH4_distance <- ggplot(dat, aes(dist_to_source_km, CH4_C_mg_m2_h)) + geom_point(aes(colour=site), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_distance


tiff("CH4_distancebyposition", units="in", width=9, height=4, res=300)
CH4_distancebyposition <- ggplot(dat, aes(dist_to_source_km, CH4_C_mg_m2_h)) + geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  legend.position = "top", legend.justification = "center", axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + scale_colour_manual(values=c("#91278e","#26867c" ,"#f7921e")) + ylab(expression(mg~CH[4]*`-C`~m^-2~h^-1)) + xlab("Distance to the source (km)") + facet_wrap(~season)
CH4_distancebyposition
dev.off()

#weird that the site upstream the dam, TA11 has high values for all 3 campaigns


#total C in CO2 equivalents and distance to source
CO2e_distance <- ggplot(dat, aes(dist_to_source_km, CO2e/1000)) + geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  legend.position = "top", legend.justification = "center", axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + scale_colour_manual(values=c("#91278e","#26867c" ,"#f7921e")) +  facet_wrap(~season) + xlab("Distance to the source (km)") + ylab(expression(g~CO[2]*"-eq m"^-2*" h"^-1))  #ylab(expression("g CO"[2]*"-eq~m^-2~d^-1")) 
CO2e_distance

#combine CO2 and CH4 and CO2e distance to source by position
tiff("CO2_CH4_distancebyposition.tiff", units="in", width=9, height=11, res=300)

CO2_CH4_distancebyposition <- ggarrange(CO2_distancebyposition + theme(axis.title.x = element_blank()) ,
                                        CH4_distancebyposition + theme(axis.title.x = element_blank()),
                                        CO2e_distance, 
ncol = 1, nrow = 3, align="hv",common.legend = T,legend="top", labels = c("(a)", "(b)",  "(c)")) 
CO2_CH4_distancebyposition

dev.off()


#### GHG and distance to dam ####
CO2_distancedam <- ggplot(dat, aes(dist_ds_dam_km, CO2_C_mg_m2_h)) + geom_point(aes(colour=position), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_distancedam



CO2_distancedam <- ggplot(dat, aes(dist_to_source_km, DailyMeanWaterTemp_C)) + geom_point(aes(colour=position_d), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_distancedam

distancedam <- ggplot(dat_means, aes(dist_to_dam_km, mean_velocity_m_s)) + geom_point(aes(colour=position_d), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
distancedam


CH4_distancedam <- ggplot(dat, aes(dist_ds_dam_km, CH4_C_mg_m2_h)) + geom_point(aes(colour=site), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_distancedam

min(dat$CO2_C_mg_m2_h)

#### GHG and percent openness ####
tiff("CO2_openness", units="in", width=3, height=3, res=300)
CO2_openness <- ggplot(dat, aes(percent_open_upstream, log(CO2_C_mg_m2_h+35.92515) )) + geom_point(aes(colour=position, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(),  legend.position = "top",  legend.justification = "left", panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +
#geom_smooth(method = "lm", formula = y~x, se = F, size = 1, colour="black") + 
scale_colour_manual(values=c("#91278e", "#f7921e")) +  ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) + xlab("Upstream openness (%)")
CO2_openness
dev.off()

#### GHG and width/depth ####

CO2_width <- ggplot(dat, aes(wetted_width_m, CO2_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", formula = y ~ log(x), se = F, size = 1, colour="#1e6b63") + xlab("Wetted width (m)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_width

# Fit the linear model and log models
linear_model <- lm(CO2_C_mg_m2_h ~ wetted_width_m, data = dat)
log_model <- lm(log(CO2_C_mg_m2_h+39.37704) ~ wetted_width_m, data = dat)

# Calculate R-squared values
linear_r_squared <- summary(linear_model)$r.squared #0.06011128
log_r_squared <- summary(log_model)$r.squared #0.1331127 #so the log model fits better

tiff("CO2_depth", units="in", width=4, height=4, res=300)
CO2_depth <- ggplot(subset(dat, site %in% upstream_sites), aes(mean_depth_cm, CO2_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") + xlab("Mean depth (cm)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_depth
dev.off()

# Fit the linear model and log models
linear_model <- lm(CO2_C_mg_m2_h ~ mean_depth_cm, data = dat) #0.01478489
log_model <- lm(log(CO2_C_mg_m2_h+39.37704) ~ mean_depth_cm, data = dat) #0.04608046 
log2_model <- lm(CO2_C_mg_m2_h ~ log(mean_depth_cm), data = dat) #0.008470237 
sqrt_model <- lm(CO2_C_mg_m2_h ~ sqrt(mean_depth_cm), data = dat) #0.01242048

# Calculate R-squared values
linear_r_squared <- summary(linear_model)$r.squared 
log_r_squared <- summary(log_model)$r.squared 
log2_r_squared <- summary(log2_model)$r.squared 
sqrt_r_squared <- summary(sqrt_model)$r.squared 

CH4_width <- ggplot(dat, aes(wetted_width_m, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) +   facet_wrap(~ season, ncol = 1)  # Facet by campaign
CH4_width

tiff("CH4_depth", units="in", width=4, height=4, res=300)
CH4_depth <- ggplot(subset(dat, site %in% upstream_sites), aes(mean_depth_cm, CH4_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab("Mean depth (cm)")
CH4_depth
dev.off()


tiff("CH4_canopy", units="in", width=4, height=4, res=300)
CH4_canopy <- ggplot(subset(dat, site %in% upstream_sites), aes(canopy_cover_., CH4_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") + xlab("Canopy cover (%)") +  ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))
CH4_canopy
dev.off()

tiff("depth_season", units="in", width=3, height=4, res=300)
depth <- ggplot(dat, aes(season, mean_depth_cm)) + geom_boxplot(colour="#26867c") + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) 
depth
dev.off()

#### GHG temperature ####

#Did the temperature differ between reservoir upstream and downstream?

temp <- lm(DailyMeanWaterTemp_C ~ position_d, (data=subset(dat_means, site=="TA09" | site=="TA10" | site=="TA11" | site=="TA0" | site=="TA05"| site=="TA0")) )
summary(temp)

temp2 <- lm(DailyMeanWaterTemp_C ~ position_d, data=dat_means)
summary(temp2)


CO2_airtemp <- ggplot(dat, aes(temp_C, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_airtemp

tiff("CH4_airtemp", units="in", width=4, height=4, res=300)
CH4_airtemp <- ggplot(subset(dat, site %in% upstream_sites), aes(temp_C, CH4_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63")  + xlab("Air temperature (C)") +  ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))
CH4_airtemp
dev.off()

CO2_watertemp <- ggplot(dat, aes(DailyMeanWaterTemp_C, CO2_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") + xlab("Water temperature (C)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_watertemp

CH4_watertemp <- ggplot(dat, aes(DailyMeanWaterTemp_C, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_watertemp

watertemp <- ggplot(dat_means, aes(dist_to_source_km, DailyMeanWaterTemp_C)) + geom_point(aes(colour=position_d), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) + geom_text(aes(label = site), size = 2, hjust = -0.1, vjust = -0.2, color = "black") +  facet_wrap(~ season, ncol = 1) 
watertemp

#### GHG and velocity /discharge ####

tiff("discharge_season", units="in", width=3, height=4, res=300)
discharge_season <- ggplot(dat, aes(season, discharge_m3_s)) + geom_boxplot(colour="#26867c") + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) 
discharge_season
dev.off()

tiff("CO2_veloc", units="in", width=4, height=4, res=300)
CO2_veloc <- ggplot(subset(dat, site %in% upstream_sites), aes(mean_velocity_m_s, CO2_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63")+ xlab("Mean velocity (m/s)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_veloc
dev.off()

tiff("CO2_veloc", units="in", width=4, height=4, res=300)
CO2_veloc <- ggplot(subset(dat, site %in% downstream_sites), aes(mean_velocity_m_s, CO2_C_mg_m2_h)) + geom_point(colour="#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#50161c")+ xlab("Mean velocity (m/s)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_veloc
dev.off()


tiff("CH4_veloc", units="in", width=4, height=4, res=300)
CH4_veloc <- ggplot(dat_means, aes(discharge_m3_s, CH4_C_mg_m2_h)) + geom_point(aes(colour= fine_substrate) ,size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + 
  #geom_smooth(method = "lm", se = F, size = 1, colour="#50161c") + 
  ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab("Mean velocity (m/s)")
CH4_veloc
dev.off()

site_veloc <- ggplot(dat_means, aes(site, mean_velocity_m_s)) + geom_boxplot() + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) 
site_veloc


tiff("CO2_discharge", units="in", width=4, height=4, res=300)
CO2_discharge <- ggplot(subset(dat, site %in% downstream_sites), aes(discharge_m3_s, CO2_C_mg_m2_h)) + geom_point(colour="#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + xlab("Discharge (m3/s)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) +geom_smooth(method = "lm", se = F, size = 1, colour="#50161c")
CO2_discharge
dev.off()

tiff("CH4_discharge", units="in", width=4, height=4, res=300)
CH4_discharge <- ggplot(subset(dat, site %in% downstream_sites), aes(discharge_m3_s, CH4_C_mg_m2_h)) + geom_point(colour= "#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#50161c")+ ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab("Discharge (m3/s")
CH4_discharge
dev.off()

#### GHG and OM stock ####
tiff("CH4OMstock", units="in", width=5, height=4, res=300)
CH4OMstock <- ggplot(dat, aes(log(OM_stock_g_m2), CH4_C_mg_m2_h)) + geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + #+ geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") 
xlab("log OM stock (g/m2)") + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))
CH4OMstock
dev.off()

#### GHG and sed OM ####

CO2_sedOM <- ggplot(dat, aes(sed_OM_percent, CO2_C_mg_m2_h)) + geom_point(colour="#26867c",size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size =12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") + xlab("Sediment OM (%)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_sedOM

tiff("CH4_sedOM", units="in", width=4, height=4, res=300)
CH4_sedOM <- ggplot(subset(dat, site %in% upstream_sites), aes(sed_OM_percent, CH4_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63")  + xlab("Sediment OM (%)") + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))
CH4_sedOM
dev.off()

#### GHG fine substrate ####

tiff("CH4_finesub", units="in", width=5, height=4, res=300)
CH4_finesub <- ggplot(dat, aes(fine_substrate, CH4_C_mg_m2_h)) + geom_point(aes(colour=position_d), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + 
  #geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63")  
  xlab("Fine substrate (%)") + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))
CH4_finesub
dev.off()

#### GHG water chemistry ####

CO2_DO <- ggplot(dat, aes(DO_mg_L, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_DO

CH4_DO <- ggplot(dat, aes(DO_mg_L, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_DO

tiff("CO2_DO_sat", units="in", width=4, height=4, res=300)
CO2_DO_sat <- ggplot(dat, aes(DO_., CO2_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + xlab("Dissolved oxygen (%)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) 
  #geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") 
CO2_DO_sat
dev.off()

tiff("CH4_DO_sat", units="in", width=4, height=4, res=300)
CH4_DO_sat <- ggplot(subset(dat, site %in% downstream_sites), aes(DO_., CH4_C_mg_m2_h)) + geom_point(colour= "#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#50161c") + xlab("Dissolved oxygen (%)") + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))
CH4_DO_sat
dev.off()

CO2_pH <- ggplot(dat, aes(water_pH, CO2_C_mg_m2_h)) + geom_point(colour= "#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63")  + xlab("Water pH") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_pH

CH4_pH <- ggplot(dat, aes(water_pH, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_pH

CO2_cond <- ggplot(dat, aes(water_conductivity_us_cm, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_cond

CH4_cond <- ggplot(dat, aes(water_conductivity_us_cm, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_cond

tiff("dist_cond", units="in", width=4, height=4, res=300)
dist_cond <- ggplot(dat, aes(dist_to_source_km, water_conductivity_us_cm, )) + geom_point(colour= "#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + xlab("Distance to source (km)") + ylab("Water conductivity (us/cm)")
dist_cond
dev.off()

#### GHG and canopy cover ####
tiff("CO2_canopy", units="in", width=4, height=4, res=300)
CO2_canopy <- ggplot(subset(dat, site %in% downstream_sites), aes(canopy_cover_., CO2_C_mg_m2_h)) + geom_point(colour="#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#50161c") + xlab("Canopy cover (%)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1)) 
CO2_canopy
dev.off()

tiff("CH4_canopy", units="in", width=4, height=4, res=300)
CH4_canopy <- ggplot(subset(dat, site %in% upstream_sites), aes(canopy_cover_., CH4_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") + xlab("Canopy cover (%)") +  ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))
CH4_canopy
dev.off()

#### Decomposition #### 

tiff("k_fine", units="in", width=9, height=4, res=300)
k_fine <- ggplot(dat_means, aes(dist_to_source_km, k_dday_fine) ) +
  geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.position="top",  panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + ylab(expression(paste(italic(k), " (", dday^-1, ")")))  +  xlab("Distance to the source (km)") + labs(title = "(b) microbes") + facet_wrap(~season) +
scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + facet_wrap(~season) +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") 
k_fine
dev.off()

tiff("CO2_kcoarse", units="in", width=4, height=4, res=300)
CO2_kcoarse <- ggplot(subset(dat, site %in% downstream_sites), aes(k_dday_coarse, CO2_C_mg_m2_h)) + geom_point(colour="#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#50161c")+ xlab("k coarse (dday-1)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_kcoarse
dev.off()

tiff("CO2_kfine", units="in", width=4, height=4, res=300)
CO2_kfine <- ggplot(subset(dat, site %in% downstream_sites), aes(k_dday_fine, CO2_C_mg_m2_h)) + geom_point(colour="#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", se = F, size = 1, colour="#50161c")+ xlab("k fine (dday-1)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_kfine
dev.off()


kratio <-  ggplot(dat_means, aes(dist_to_source_km, k_ratio_f_c) ) +
  geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.position="top",  panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + ylab(expression("kfine:kcoarse"))  +  xlab("Distance to the source (km)") + facet_wrap(~season) 
kratio


tiff("k_coarse", units="in", width=9, height=4, res=300)
k_coarse <- ggplot(dat_means, aes(dist_to_source_km, k_dday_coarse) ) +
  geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.position="top",  axis.ticks.x=element_blank(), legend.title=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + ylab(expression(paste(italic(k), " (", dday^-1, ")")))  +  xlab("Distance to the source (km)") + facet_wrap(~season) + labs(title = "(a) invertebrates + microbes") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1) +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed")  
k_coarse
dev.off()


k_coarse_summer <- ggplot(subset(dat_means, season=="Summer" & k_dday_coarse <=0.013), aes(dist_to_source_km, k_dday_coarse) ) +   geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.position="top",  axis.ticks.x=element_blank(), legend.title=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + ylab(expression(paste(italic(k), " (", dday^-1, ")")))  +  xlab("Distance to the source (km)")
k_coarse_summer

tiff("k_fineOM", units="in", width=5, height=4, res=300)
k_fineOM <- ggplot(dat, aes(sed_OM_percent, k_dday_fine) ) +
  geom_point(aes(colour=position), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + scale_colour_manual(values=c("#91278e", "#f7921e")) + ylab(expression(paste(italic(k), " fine (", dday^-1, ")"))) + geom_smooth(data = subset(dat_means, position %in% c("upstream", "downstream")), method = "lm", se = F, aes(group = position, colour = position), size = 1)  + xlab("Sediment OM (%)")
k_fineOM
dev.off()

tiff("kcoarse_cond", units="in", width=5, height=4, res=300)
kcoarse_cond <- ggplot(dat, aes(water_conductivity_us_cm, k_dday_coarse)) + geom_point(aes(colour=position), size=3.5, alpha=0.3) +  theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + ylab("k coarse (dday-1)") + scale_colour_manual(values=c("#91278e", "#f7921e")) + geom_smooth(data = subset(dat_means, position %in% c("upstream", "downstream")), method = "lm", se = F, aes(group = position, colour = position), size = 1) + xlab("Water conductivity (us/cm)")
kcoarse_cond
dev.off()

masl_cond <- ggplot(dat, aes(masl, k_dday_coarse)) + geom_point(aes(colour=position), size=3.5, alpha=0.3) +  theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + ylab("k coarse (dday-1)") + scale_colour_manual(values=c("#91278e", "#f7921e")) + geom_smooth(data = subset(dat_means, position %in% c("upstream", "downstream")), method = "lm", se = F, aes(group = position, colour = position), size = 1) + xlab("Water conductivity (us/cm)") 
masl_cond

lm_model <- lm(k_dday_coarse ~ masl, data = dat)
summary(lm_model)

kcoarse_dist <- ggplot(dat, aes(dist_to_dam_km, k_dday_coarse)) + geom_point(aes(colour=position), size=3.5, alpha=0.3) +  theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + ylab("k coarse (dday-1)") + scale_colour_manual(values=c("#91278e", "#f7921e")) + geom_smooth(data = subset(dat_means, position %in% c("upstream", "downstream")), method = "lm", se = F, aes(group = position, colour = position), size = 1) 
kcoarse_dist

## Combine the k plots ##

tiff("k_fine_coarse", units="in", width=9, height=8, res=300)

k_fine_coarse <- ggarrange(k_coarse + theme(axis.title.x = element_blank() ),  
                           k_fine, 
                           ncol = 1, nrow = 2, align="hv",common.legend = T,legend="top")
k_fine_coarse

dev.off()


#### For Fanny: k per day (not dday) ####


kday_fine <- ggplot(dat_means, aes(dist_to_source_km, k_day_fine) ) +
  geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.position="top",  panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + ylab(expression(paste(italic(k), " (", day^-1, ")")))  +  xlab("Distance to the source (km)") + labs(title = "(b) microbes") +    scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + facet_wrap(~season) +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="solid") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") 
kday_fine

kday_coarse <- ggplot(dat_means, aes(dist_to_source_km, k_day_coarse) ) +
  geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.title=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.position="top",  panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + ylab(expression(paste(italic(k), " (", day^-1, ")")))  +  xlab("Distance to the source (km)") + labs(title = "(a) invertebrates + microbes") +    scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + facet_wrap(~season) +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1,  linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") 
kday_coarse



## Combine the k day plots ##

tiff("kday_fine_coarse.tiff", units="in", width=9, height=8, res=300)

kday_fine_coarse <- ggarrange(kday_coarse + theme(axis.title.x = element_blank() ),  
                              kday_fine, 
                           ncol = 1, nrow = 2, align="hv",common.legend = T,legend="top")
kday_fine_coarse

dev.off()





####OM : Decomposition#### 
#new variable that is the ratio of OM stock to decomp rate
dat_means$coarse_OM <- (dat_means$OM_stock_g_m2/dat_means$k_day_coarse)/1000
dat_means$fine_OM <- (dat_means$OM_stock_g_m2/dat_means$k_day_fine)/1000

kfineOM_dist <- ggplot(dat_means, aes(dist_to_source_km, fine_OM)) + geom_point(aes(colour=position_d), size=3.5, alpha=0.6) +  theme_bw() + theme(legend.position="top", legend.title=element_blank(), axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + ylab("OM stock : k microbes") + scale_colour_manual(values=c("#91278e", "#1e6b63","#f7921e")) + facet_wrap(~season) + ylim(0,4.5) + xlab("Distance to source (km)") +
 geom_text(data = subset(dat_means, position_d == "reservoir"), aes(label = site), size = 2.5, vjust = 1.5, hjust = -0.5) +
  geom_text(data = subset(dat_means, season == "Summer"), aes(x = 16.620201, y = 4, label = sprintf('\u2191')), colour="#1e6b63") +
  geom_text(data = subset(dat_means, season == "Summer"), aes(x = 20, y = 4, label = "TA10"), size=2.5)  + 
  labs(title = "(b) microbes") 
kfineOM_dist # Note that the ratio is divided by 1000


kcoarseOM_dist <- ggplot(dat_means, aes(dist_to_source_km, coarse_OM)) + geom_point(aes(colour=position_d), size=3.5, alpha=0.6) +  theme_bw() + theme(legend.position="top", legend.title=element_blank(), axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + ylab("OM stock : k total") + scale_colour_manual(values=c("#91278e", "#1e6b63","#f7921e")) + facet_wrap(~season) + ylim(0,4.5)  + xlab("Distance to source (km)") +
  geom_text(data = subset(dat_means, position_d == "reservoir"), aes(label = site), size = 2.5, vjust = 1.5, hjust = -0.5) +
  geom_text(data = subset(dat_means, season == "Summer"), aes(x = 16.620201, y = 3.8, label = sprintf('\u2191')), colour="#1e6b63")  +
 geom_text(data = subset(dat_means, season == "Summer"), aes(x = 20, y = 3.8, label = "TA10"), size=2.5) + labs(title = "(a) microbes + invertebrates") 
kcoarseOM_dist 

#combine

tiff("k_OM_ratio", units="in", width=6.5, height=7, res=300)

k_fine_coarse <- ggarrange(kcoarseOM_dist + theme(axis.title.x = element_blank() ),  
                           kfineOM_dist, 
                           ncol = 1, nrow = 2, align="hv",common.legend = T,legend="top")
k_fine_coarse

dev.off()

#### OM #### 

#CO2 vs OM flux #26867c

tiff("CO2_OMflux", units="in", width=4, height=4, res=300)
CO2_OMflux <- ggplot(subset(dat, site %in% upstream_sites), aes(OM_flux_g_m2_s, CO2_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + geom_smooth(method = "lm", formula = y ~ x, se = F, size = 1, colour="#1e6b63") + xlab("OM flux (g/m2/s)") + ylab(expression(mg~CO[2]*`-C`~m^-2*~h^-1))
CO2_OMflux
dev.off()

tiff("OMflux_season", units="in", width=3, height=4, res=300)
OMflux_season <- ggplot(dat, aes(season, OM_flux_g_m2_s)) + geom_boxplot(colour="#26867c") + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) 
OMflux_season
dev.off()

tiff("CH4_OMstock", units="in", width=4, height=4, res=300)
CH4_OMstock <- ggplot(subset(dat, site %in% upstream_sites), aes(OM_stock_g_m2, CH4_C_mg_m2_h)) + geom_point(colour="#26867c", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +geom_smooth(method = "lm", se = F, size = 1, colour="#1e6b63") + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab("OM stock (g/m2)")
CH4_OMstock
dev.off()

tiff("CH4_OMstock", units="in", width=4, height=4, res=300)
CH4_OMstock <- ggplot(subset(dat, site %in% downstream_sites), aes(OM_stock_g_m2, CH4_C_mg_m2_h)) + geom_point(colour="#862630", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +geom_smooth(method = "lm", se = F, size = 1, colour="#50161c") + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab("OM stock (g/m2)")
CH4_OMstock
dev.off()

tiff("OM_stock", units="in", width=9, height=4, res=300)
dist_OMstock <- ggplot(dat_means, aes(dist_to_source_km, OM_stock_g_m2 )) + geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), legend.position="top",  legend.title=element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + ylab(expression("OM stock (g m"^{-2} * ")")) + xlab("Distance to source (km)") + scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) + facet_wrap(~season) +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1) +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") + scale_y_log10(labels = label_number()) 
dist_OMstock
dev.off()

tiff("OM_stock", units="in", width=6, height=4, res=300)
OM_stock <-  ggplot(dat_means, aes(dist_to_source_km, log(OM_stock_g_m2) )) + geom_point(aes(colour=position, shape=season), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + scale_colour_manual(values=c("#91278e", "#f7921e"))  + ylab(expression("log OM stock (g m"^{-2} * ")")) +
xlab("Distance to the source (km)") + 
  geom_smooth(aes(group=season, linetype=season), data = subset(dat_means, position == "upstream"),
    method = "gam",  formula = y ~ s(x),  # Gam smooth relationship
    se = FALSE, colour = "#91278e", linewidth = 1) + 
geom_smooth(aes(group=season,  linetype=season), data = subset(dat_means, position == "downstream"), method = "lm",   formula = y ~ poly(x, 2),  # Linear regression for downstream
    se = FALSE,colour = "#f7921e", linewidth = 1)
OM_stock
dev.off()

tiff("OM_stock", units="in", width=6, height=4, res=300)
OM_stock <-  ggplot(dat_means, aes(dist_to_source_km, log(OM_stock_g_m2) )) + geom_point(aes(colour=position_d), size=3.5, alpha=0.6) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title=element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + scale_colour_manual(values=c("#186a84", "#ef4345", "#186a84"))  + ylab(expression("log OM stock (g m"^{-2} * ")")) +
  xlab("Distance to the source (km)") + 
  geom_smooth(data = subset(dat_means, position == "upstream"),colour="#ef4345",
              method = "gam",  formula = y ~ s(x),  # Gam smooth relationship
              se = FALSE, linewidth = 1) + 
  geom_smooth(data = subset(dat_means, position == "downstream"), colour="#ef4345", method = "lm",   formula = y ~ poly(x, 2),  # Linear regression for downstream
              se = FALSE, linewidth = 1)
OM_stock
dev.off()


OM_stock_up <-  ggplot(subset(dat, site %in% upstream_sites), aes(dist_to_source_km, log(OM_stock_g_m2) )) + geom_point(aes( shape=season), colour="#91278e", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))+ ylab(expression("log OM stock (g m"^{-2} * ")")) +
  xlab("Distance to the source (km)") + 
  geom_smooth(data = subset(dat_means, position == "upstream"),
              method = "lm",  formula = y ~ poly(x, 2),  # Polynomial regression for upstream
              se = FALSE, colour = "#91278e", linewidth = 1  ) + stat_regline_equation(label.y = 6) + stat_cor(label.y = 5) 
OM_stock_up

OM_stock_up <-  ggplot(subset(dat, site %in% upstream_sites), aes(dist_to_source_km, log(OM_stock_g_m2) )) + geom_point(aes( shape=season), colour="#91278e", size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))+ ylab(expression("log OM stock (g m"^{-2} * ")")) +
  xlab("Distance to the source (km)") + 
  geom_smooth(data = subset(dat_means, position == "upstream"),
              method = "lm",  se = FALSE,  formula = y ~ x, colour = "#91278e", linewidth = 1  ) + stat_regline_equation(label.y = 6) + stat_cor(label.y = 5) 
OM_stock_up

tiff("OM_stock_open_d", units="in", width=6, height=4, res=300)
OM_stock_open_d <-  ggplot(dat_means, aes(percent_open_downstream, log(OM_stock_g_m2) )) + geom_point(aes(colour=season, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))+ ylab(expression("log OM stock (g m"^{-2} * ")")) +
  xlab("open downstream (%)")+
  geom_smooth(aes(linetype=season, colour=season), data = dat_means, method = "lm",  formula = y ~ x,  se = FALSE, linewidth = 1  ) + scale_linetype_manual(values = c("dotted","dashed", "solid"))
OM_stock_open_d
dev.off()

tiff("OM_stock_open_up", units="in", width=6, height=4, res=300)
OM_stock_open_up <-  ggplot(dat_means, aes(percent_open_upstream, log(OM_stock_g_m2) )) + geom_point(aes(colour=season, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))+ ylab(expression("log OM stock (g m"^{-2} * ")")) +
  xlab("open upstream (%)")+
  geom_smooth(aes(linetype=season, colour=season), data = dat_means, method = "lm",  formula = y ~ x,  se = FALSE, linewidth = 1  ) + scale_linetype_manual(values = c("dotted","dashed", "twodash"))
OM_stock_open_up
dev.off()



dist_OMtransport <-  ggplot(dat_means, aes(dist_to_source_km, OM_g_m3) ) + geom_point(aes(colour=position_d, shape=season), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.title=element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), legend.position="top", axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  +  scale_colour_manual(values=c("#91278e", "#1e6b63", "#f7921e")) +ylab(expression(paste("OM transport (g m"^{-3} * ")"))) +  xlab("Distance to the source (km)") +  facet_wrap(~ season) +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Spring"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Summer"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "upstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#91278e", linewidth = 1, linetype="dashed") +
  geom_smooth(data = subset(dat_means, position_d == "downstream" & season=="Fall"), method = "lm",  se = FALSE,colour = "#f7921e", linewidth = 1, linetype="dashed") +  scale_y_log10(labels = label_number()) 
dist_OMtransport

## Combine the OM plots ##

tiff("OM_stock_flux.tiff", units="in", width=9, height=8, res=300)

OM_stock_flux <- ggarrange(dist_OMstock + theme(axis.title.x = element_blank()),  
                           dist_OMtransport, 
                          ncol = 1, nrow = 2, align="hv",common.legend = T,legend="top",  
                          labels = c("(a)", "(b)"))
OM_stock_flux

dev.off()

## Combine plots for CO2 drivers from RF ##

tiff("CO2_inf_vars", units="in", width=4, height=6, res=300)

CO2_inf_vars <- ggarrange(CO2_DO_sat + theme(axis.title.x = element_blank()),  
                          CO2_pH + theme(axis.title.x = element_blank()), 
                          CO2_sedOM + theme(axis.title.x = element_blank()), 
                          CO2_watertemp + theme(axis.title.x = element_blank()), 
                          ncol = 1, nrow = 2, align="hv",common.legend = T,legend="top",  
                          labels = c("(a)", "(b)"))
CO2_inf_vars

dev.off()


############# Make network maps with variables ###############
#### Script for making Konza STIC map

# load in shapefile
waterways <- st_read("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Site selection/waterways_clipped/tarentaine.shp")

catchment_area <- st_read("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/ArcMap/05_Adour-Garonne_BassinVersantTopographique.shp/05_Adour-Garonne_BassinVersantTopographique.shp")

#subset sub catchment
# IDs to filter
selected_ids <- c("05B0000002150458641", "05B0000002150459669", "05B0000002150460843", "05B0000002150460858", "05B0000002150456990")

# Filter catchment areas based on selected IDs
tarentaine_catch <- catchment_area[catchment_area$CdOH %in% selected_ids, ]

tarentaine <- st_intersection(waterways, tarentaine_catch)

sites <- st_as_sf(dat_means, coords = c("lon", "lat"), crs = 4326)


# Subset C1
C1 <- subset(dat_means, campaign==1)
sites_C1 <- st_as_sf(C1, coords = c("lon", "lat"), crs = 4326)
# Subset C2
C2 <- subset(dat_means, campaign==2)
sites_C2 <- st_as_sf(C2, coords = c("lon", "lat"), crs = 4326)
# Subset C3
C3 <- subset(dat_means, campaign==3)
sites_C3 <- st_as_sf(C3, coords = c("lon", "lat"), crs = 4326)

#Plot CO2 
CO2_map <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites, aes(colour = CO2_C_mg_m2_h),  size = 3, alpha = 0.3) +
  scale_colour_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1) +
  theme_classic()+ ggtitle("CO2") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
CO2_map

#Plot CH4 
CH4_map <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites, aes(colour = CH4_C_mg_m2_h),  size = 3, alpha = 0.3) +
  scale_colour_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1, trans = "log") +
  theme_classic()+ ggtitle("CH4") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
CH4_map #note the log scale


CO2spring <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites_C1, aes(colour = CO2_C_mg_m2_h),  size = 3, alpha = 0.3) +
  scale_colour_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1) +
  theme_classic()+ ggtitle("CO2 stock Spring") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
CO2spring

CO2summer <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites_C2, aes(colour = CO2_C_mg_m2_h),  size = 3, alpha = 0.3) +
  scale_colour_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1) +
  theme_classic()+ ggtitle("CO2 stock Summer") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
CO2summer

CO2fall <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites_C3, aes(colour = CO2_C_mg_m2_h),  size = 3, alpha = 0.3) +
  scale_colour_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1) +
  theme_classic()+ ggtitle("CO2 stock Fall") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
CO2fall


#plot 
OMspring <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites_C1, aes(fill = OM_stock_g_m2, shape=flow_state),  size = 3, alpha = 0.3, colour="black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1) +
  theme_classic()+ ggtitle("OM stock Spring") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
OMspring



#plot 
OMsummer <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites_C2, aes(fill = OM_stock_g_m2, shape=flow_state, size=OM_stock_g_m2) , alpha = 0.3, colour="black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1, trans = "log") +
  theme_classic()+ ggtitle("OM stock Summer") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
OMsummer



#plot 
OMfall <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites_C3, aes(fill = OM_stock_g_m2, shape=flow_state),  size = 4, alpha = 0.3, colour="black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1) +
  theme_classic()+ ggtitle("OM stock Fall") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 
OMfall

##############################################################################

##### Run Random Forest Model on the data ####

# Goal: to use the sites upstream of the 2 main dams, find the drivers of GHG fluxes and OM, and then using that relationship predict the downstream values, then compare those to the actual impacted values to see what effect the dams have. 



#### run RF for CO2 ####

# Subset 
dat_rf <- dat %>% 
 # filter(site %in% upstream_sites)  %>%
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel)

#str(dat_upstream)
str(dat_rf)


# Convert the tibble to a data frame
#dat_upstream <- as.data.frame(dat_upstream)
dat_rf <- as.data.frame(dat_rf)

#The training data is used for building a model, while the testing data is used for making predictions. This means after fitting a model on the training data set, finding of the errors and minimizing those error, the model is used for making predictions on the unseen data which is the test data.

#split <- sample.split(dat_upstream, SplitRatio = 0.8) 
#split 

#can also trying using the entire data set

split <- sample.split(dat_rf, SplitRatio = 0.8) 
split 

#The split method splits the data into train and test datasets with a ratio of 0.8 This means 80% of our dataset is passed in the training dataset and 20% in the testing dataset.

data_train <- subset(dat_rf, split == "TRUE") 
data_test <- subset(dat_rf, split == "FALSE") 

any_na <- any(is.na(data_train))

na_columns <- colnames(data_train)[apply(data_train, 2, anyNA)]

dim(data_train)

str(data_train)


#Find optimized value of number of random variables
bestmtry <- tuneRF(data_train,data_train$CO2_C_mg_m2_h,stepFactor = 1.2, improve = 0.01, trace=T, plot= T)  #9

#create RF model: 
model <- randomForest(CO2_C_mg_m2_h~.,data= data_train, ntree = 1000, mtry = 9, importance = TRUE) #need to determine what is the best value for mtry and ntree
model 

plot(model)

imp <- randomForest::importance(model)  # returns the importance of the variables

varImpPlot(model)  # visualizing the importance of variables of the model.

impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]

# Select the top 6 variables
top_variables <- impvar[1:6]

op <- par(mfrow=c(2, 3))

for (i in seq_along(top_variables)) {
  partialPlot(model, data_train, impvar[i], xlab=impvar[i],
              main=paste("Partial Dependence on", impvar[i]),
              ylim=c(30, 70))
}
par(op)

# Visualize variable importance 

# Get variable importance from the model fit
ImpData <- as.data.frame(importance(model))
ImpData$Var.Names <- row.names(ImpData)

ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


#Make predictions

pred_test <- predict(model, newdata = data_test)

table(pred_test, data_test$CO2_C_mg_m2_h)

#The issue is that RF cannot extrapolate outside of it's dataset, so a linear regression might be better


#Evaluation of predictions using RMSE
sqrt(mean((pred_test - data_test$CO2_C_mg_m2_h)^2))
#63.63705

#Goodness of fit of your regression model, you can calculate the R-squared value:
1 - sum((data_test$CO2_C_mg_m2_h - pred_test)^2) / sum((data_test$CO2_C_mg_m2_h - mean(data_test$CO2_C_mg_m2_h))^2)
#0.5937679



#### run RF for CH4 ####

# Subset 
dat_rf <- dat %>% 
  # filter(site %in% upstream_sites)  %>%
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel)

#str(dat_upstream)
str(dat_rf)


# Convert the tibble to a data frame
#dat_upstream <- as.data.frame(dat_upstream)
dat_rf <- as.data.frame(dat_rf)

#The training data is used for building a model, while the testing data is used for making predictions. This means after fitting a model on the training data set, finding of the errors and minimizing those error, the model is used for making predictions on the unseen data which is the test data.

#split <- sample.split(dat_upstream, SplitRatio = 0.8) 
#split 

#can also trying using the entire data set

split <- sample.split(dat_rf, SplitRatio = 0.8) 
split 

#The split method splits the data into train and test datasets with a ratio of 0.8 This means 80% of our dataset is passed in the training dataset and 20% in the testing dataset.

data_train <- subset(dat_rf, split == "TRUE") 
data_test <- subset(dat_rf, split == "FALSE") 

any_na <- any(is.na(data_train))

na_columns <- colnames(data_train)[apply(data_train, 2, anyNA)]

dim(data_train)

str(data_train)

#Find optimized value of number of random variables
bestmtry <- tuneRF(data_train,data_train$CH4_C_mg_m2_h,stepFactor = 1.2, improve = 0.01, trace=T, plot= T)  #9

#create RF model for CH4
modelCH4 <- randomForest(CH4_C_mg_m2_h~.,data= data_train, ntree = 400, mtry = 9, importance = TRUE) #need to determine what is the best value for mtry and ntree
modelCH4 

plot(modelCH4)

imp <- randomForest::importance(modelCH4)  # returns the importance of the variables

varImpPlot(modelCH4)  # visualizing the importance of variables of the model.

impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]

# Select the top 6 variables
top_variables <- impvar[1:6]

op <- par(mfrow=c(2, 3))

for (i in seq_along(top_variables)) {
  partialPlot(model, data_train, impvar[i], xlab=impvar[i],
              main=paste("Partial Dependence on", impvar[i]),
              ylim=c(30, 70))
}
par(op)

# Visualize variable importance 

# Get variable importance from the model fit
ImpData <- as.data.frame(sw(modelCH4))
ImpData$Var.Names <- row.names(ImpData)

ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )


################################################################################
#### Test collinearity ####

#Before running an LMM, need to test the predictors for multicollinearity

# Check colinearity for CO2 
#Subset just the upstream sites, TA06, TA11, TA07, TA08, TA12 and higher
upstream_sites <- c("TA06", "TA10", "TA05", "TA11", "TA07", "TA08", "TA12", "TA13", "TA14", "TA15", "TA17", "TA20", "TA21", "TA22", "TA24", "TA02R")
#15 sites
downstream_sites <- c("TA04", "TA09", "TA01", "TA02", "TA03")
#7 sites

#subset relevant columns and filter just the upstream/downstream sites
dat_lm_CO2_up <- dat %>% 
  filter(site %in% upstream_sites)  %>%
  select(-site, -season, -date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position)

dat_lm_CH4_up <- dat %>% 
  filter(site %in% upstream_sites)  %>%
  select(-site, -season, -date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position)


dat_lm_CO2_down <- dat %>% 
  filter(site %in% downstream_sites)  %>%
  select(-site, -season, -date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -bedrock) #note bedrock is all 0 downstream, so remove

dat_lm_CH4_down <- dat %>% 
  filter(site %in% downstream_sites)  %>%
  select(-site, -season, -date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -bedrock) #note bedrock is all 0 downstream, so remove

dat_lm_OMstock_up <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-site, -season, -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position) #remove some of the water chem data as it does not make sense ecologically 

dat_lm_OMstock_down <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select(-site, -season, -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -bedrock) #note to remove bedrock

dat_lm_OMflux_up <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-site, -season, -date, -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position)

dat_lm_OMflux_down <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select(-site, -season, -date, -campaign,  -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -bedrock) #note to remove bedrock

dat_lm_kcoarse_up <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-site, -season, -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position)

dat_lm_kcoarse_down <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select(-site, -season, -date, -campaign,  -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -bedrock, -DO_mg_L, -dist_ds_dam_km) #too many variables, need to remove some


#create a linear model with all of the data. 
CO2_model_collinearity <- lm(CO2_C_mg_m2_h~., data = dat_lm_CO2_up) 
summary(CO2_model_collinearity)

CH4_model_collinearity <- lm(CH4_C_mg_m2_h~., data = dat_lm_CH4_up) 
summary(CH4_model_collinearity)

CO2_model_collinearity_d <- lm(CO2_C_mg_m2_h~., data = dat_lm_CO2_down) 
summary(CO2_model_collinearity_d)

CH4_model_collinearity_d <- lm(CH4_C_mg_m2_h~., data = dat_lm_CH4_down) 
summary(CH4_model_collinearity_d)

OMstock_model_collinearity <- lm(OM_stock_g_m2 ~., data = dat_lm_OMstock_up) 
summary(OMstock_model_collinearity)

OMstock_model_collinearity_d <- lm(OM_stock_g_m2 ~., data = dat_lm_OMstock_down) 
summary(OMstock_model_collinearity_d)

OMflux_model_collinearity <- lm(OM_flux_g_m2_s ~., data = dat_lm_OMflux_up) 
summary(OMflux_model_collinearity)

OMflux_model_collinearity_d <- lm(OM_flux_g_m2_s ~., data = dat_lm_OMflux_down) 
summary(OMflux_model_collinearity_d)

k_coarse_model_collinearity <- lm(k_dday_coarse ~., data = dat_lm_kcoarse_up) 
summary(k_coarse_model_collinearity)

k_coarse_model_collinearity_d <- lm(k_dday_coarse ~., data = dat_lm_kcoarse_down) 
summary(k_coarse_model_collinearity_d)


#Test tolerance and VIF. Tolerance measures the percent of the variance in the independent variable that cannot be accounted for by the other independent variables. The Variance Inflation Factor (VIF) measures the inflation in the coefficient of the independent variable due to the collinearities among the other independent variables (VIF>10 indicates multicollinearity).

vif_df <-as.data.frame(ols_vif_tol(CO2_model_collinearity))
vif_df_d <-as.data.frame(ols_vif_tol(CO2_model_collinearity_d))
vif_df_CH4 <-as.data.frame(ols_vif_tol(CH4_model_collinearity))
vif_df_d_CH4 <-as.data.frame(ols_vif_tol(CH4_model_collinearity_d))
vif_df_OMstock <-as.data.frame(ols_vif_tol(OMstock_model_collinearity))
vif_df_OMstock_d <-as.data.frame(ols_vif_tol(OMstock_model_collinearity_d))
vif_df_OMflux <-as.data.frame(ols_vif_tol(OMflux_model_collinearity))
vif_df_OMflux_d <-as.data.frame(ols_vif_tol(OMflux_model_collinearity_d))
vif_df_kcoarse <-as.data.frame(ols_vif_tol(k_coarse_model_collinearity))
vif_df_kcoarse_d <-as.data.frame(ols_vif_tol(k_coarse_model_collinearity_d))

# Order the data frame by VIF values in decreasing order

vif_df <- vif_df[order(vif_df$VIF, decreasing = TRUE), ]
#Variables with VIF >10 are DO_mg_L, water_temp_C, dist_ds_dam_km, DO_., dist_to_source_km, masl
#We can remove DO_mg_L, as it has the highest VIF, and it is still represented in DO %
#we can remove water_temp_C as it is represented in the iButton water temp
#we can remove dist_ds_dam_km as it is represented in width and depth and discharge and masl

vif_df_d <- vif_df_d[order(vif_df_d$VIF, decreasing = TRUE), ]
#variables with VIF >10 are dist_ds_dam_km, Distance_to_sourc_km, DO_mg_L, water_temp_C, masl, DO %, DailyMeanWaterTemp_C, OM_flux_g_m2_s, temp_C, fine_substrate, sed_OM_percent
#We can remove dist_ds_dam_km and dist_to_source_km as they are represented in masl
#We can remove DO_mg_L, as it is still represented in DO %
#we can remove water_temp_C as it is represented in the iButton DailyMeanWaterTemp_C
#can remove air temperature (temp_C) as water temperature is more representative

vif_df_CH4 <- vif_df_CH4[order(vif_df_CH4$VIF, decreasing = TRUE), ]
#variables with VIF >10 are DO_mg_L, water_temp_C, DO_., dist_ds_dam_km, masl, dist_to_source_km, #we can remove distance to dam and masl as they are represented in distance to source
#We can remove water_temp_C as it's represented in DailyMeanWaterTemp_C

vif_df_d_CH4 <- vif_df_d_CH4[order(vif_df_d_CH4$VIF, decreasing = TRUE), ]
#variables with VIF >10 are dist_ds_dam_km, dist_to_source_km, DO_mg_L, water_temp_C, masl, DO_., DailyMeanWaterTemp_C, OM_flux_g_m2_s, temp_C, fine_substrate, sed_OM_percent 
#Can remove dist ds dam and dist to source because they are represented in masl
#We can remove DO_mg_L, as it is still represented in DO %
#we can remove water_temp_C as it is represented in the iButton DailyMeanWaterTemp_C
#can remove temp_C, as it is represented in water temperature

vif_df_OMstock <- vif_df_OMstock[order(vif_df_OMstock$VIF, decreasing = TRUE), ]
# Variables with VIF >10 are DO_mg_L, water_temp_C, DO_., dist_ds_dam_km, masl, dist_to_source_km, water_conductivity_us_cm 
# Can remove DO_mg_L as it's represented in DO%
# Can remove water temp C, is it's represented in DailyMeanWaterTemp
#can remoce dist ds dam and masl is it's represented in distance to source

vif_df_OMstock_d <- vif_df_OMstock_d[order(vif_df_OMstock_d$VIF, decreasing = TRUE), ]
# Variables with VIF >10 are dist_ds_dam_km, dist_to_source_km, masl, DO_mg_L, water_temp_C, DO_., DailyMeanWaterTemp_C, OM_flux_g_m2_s, fine_substrate, temp_C, sed_OM_percent 
#Can remove dist ds dam and dist to source, as it's represented in masl
#Can remove DO_mg_L as it's represented in DO%
# Can remove water temp C as it's represented in Daily mean water temp

vif_df_OMflux <- vif_df_OMflux[order(vif_df_OMflux$VIF, decreasing = TRUE), ]
# Variables with VIF >10 are DO_mg_L, water_temp_C, DO_., dist_ds_dam_km, masl, dist_to_source_km, water_conductivity_us_cm 
# Can remove DO_mg_L as it's represented in DO%
# Can remove water temp C as it's represented in Daily mean water temp
#Can remove dist ds dam and masl, as it's represented in dist to source

vif_df_OMflux_d <- vif_df_OMflux_d[order(vif_df_OMflux_d$VIF, decreasing = TRUE), ]
# Variables with VIF >10 are dist_ds_dam_km dist_to_source_km masl DO_mg_L water_temp_C DailyMeanWaterTemp_C DO_. temp_C fine_substrate 
#Can remove dist_ds_dam_km and dist_to_source_km as they are represented in masl
# Can remove DO_mg_L as it's represented in DO%
# Can remove water temp C as it's represented in Daily mean water temp

vif_df_kcoarse <- vif_df_kcoarse[order(vif_df_kcoarse$VIF, decreasing = TRUE), ]
# Variables with VIF >10 are DO_mg_L water_temp_C DO_. dist_ds_dam_km dist_to_source_km masl discharge_m3_s wetted_width_m water_conductivity_us_cm mean_velocity_m_s 
# Can remove DO_mg_L as it's represented in DO%
# Can remove water temp C as it's represented in Daily mean water temp
#Can remove dist_ds_dam_km and distance to source as it's represented in masl

vif_df_kcoarse_d <- vif_df_kcoarse_d[order(vif_df_kcoarse_d$VIF, decreasing = TRUE), ]
#Many vars are >10 VIF
#remove distance to source
#remove fine substrate
#remove DailyMeanWaterTemp_C 

#### Rerun the ols_vif_tol analysis with these vars removed ####

dat_colin <- dat %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -position)

dat_colin_d <- dat %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -bedrock, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -temp_C, -fine_substrate)

dat_colinCH4 <- dat %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -masl, -position)

dat_colin_d_CH4 <- dat %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -bedrock, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -temp_C)

dat_colin_OMstock <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -site, -campaign, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -masl)

dat_colin_OMstock_d <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date, -site, -campaign,  -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -bedrock, -DailyMeanWaterTemp_C)

dat_colin_OMflux <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -site,   -campaign,  -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -masl)

dat_colin_OMflux_d <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date, -site,  -campaign, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -bedrock, -DailyMeanWaterTemp_C)

dat_colin_OMflux <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -site,  -campaign, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -masl)

dat_colin_kcoarse <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-site, -season, -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -discharge_m3_s)

dat_colin_kcoarse_d <- dat_means %>% 
    filter(site %in% downstream_sites)  %>%
    select(-site, -season, -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -DailyMeanWaterTemp_C, -fine_substrate, -bedrock, -mean_velocity_m_s, -mean_depth_cm)


#Try again

#create a linear model with all of the data. #try with dat_colin to check if removal improved VIFs
CO2_model_collinearity2 <- lm(CO2_C_mg_m2_h~., data = dat_colin) 
summary(CO2_model_collinearity2)

CO2_model_collinearity2_d <- lm(CO2_C_mg_m2_h~., data = dat_colin_d) 
summary(CO2_model_collinearity2_d)

CH4_model_collinearity2 <- lm(CH4_C_mg_m2_h~., data = dat_colinCH4) 
summary(CH4_model_collinearity2)

CH4_model_collinearity2_d_CH4 <- lm(CH4_C_mg_m2_h~., data = dat_colin_d_CH4) 
summary(CH4_model_collinearity2_d_CH4)

OMstock_model_collinearity2 <- lm(OM_stock_g_m2~., data = dat_colin_OMstock) 
summary(OMstock_model_collinearity2)

OMstock_d_model_collinearity2 <- lm(OM_stock_g_m2~., data = dat_colin_OMstock_d) 
summary(OMstock_d_model_collinearity2)

OMflux_model_collinearity2 <- lm(OM_stock_g_m2~., data = dat_colin_OMflux) 
summary(OMflux_model_collinearity2)

OMflux_d_model_collinearity2 <- lm(OM_stock_g_m2~., data = dat_colin_OMflux_d) 
summary(OMflux_d_model_collinearity2)

kcoarse_model_collinearity2 <- lm(k_dday_coarse~., data = dat_colin_kcoarse) 
summary(kcoarse_model_collinearity2)

kcoarse_d_model_collinearity2 <- lm(k_dday_coarse~., data = dat_colin_kcoarse_d) 
summary(kcoarse_d_model_collinearity2)

#Test tolerance and VIF. Tolerance measures the percent of the variance in the independent variable that cannot be accounted for by the other independent variables. The Variance Inflation Factor (VIF) measures the inflation in the coefficient of the independent variable due to the collinearities among the other independent variables (VIF>10 indicates multicollinearity).

vif_df <-as.data.frame(ols_vif_tol(CO2_model_collinearity2))
vif_df_d <-as.data.frame(ols_vif_tol(CO2_model_collinearity2_d))
vif_dfCH4 <-as.data.frame(ols_vif_tol(CH4_model_collinearity2))
vif_df_d_CH4 <-as.data.frame(ols_vif_tol(CH4_model_collinearity2_d_CH4))
vif_df_OMstock <-as.data.frame(ols_vif_tol(OMstock_model_collinearity2))
vif_df_OMstock_d <-as.data.frame(ols_vif_tol(OMstock_d_model_collinearity2))
vif_df_OMflux <-as.data.frame(ols_vif_tol(OMflux_model_collinearity2))
vif_df_OMflux_d <-as.data.frame(ols_vif_tol(OMflux_d_model_collinearity2))
vif_df_kcoarse <-as.data.frame(ols_vif_tol(kcoarse_model_collinearity2))
vif_df_kcoarse_d <-as.data.frame(ols_vif_tol(kcoarse_d_model_collinearity2))

# Order the data frame by VIF values in decreasing order
vif_df <- vif_df[order(vif_df$VIF, decreasing = TRUE), ] #good! 
vif_df_d <- vif_df_d[order(vif_df_d$VIF, decreasing = TRUE), ] #good!
vif_dfCH4 <- vif_dfCH4[order(vif_dfCH4$VIF, decreasing = TRUE), ] #good!
vif_df_d_CH4 <- vif_df_d_CH4[order(vif_df_d_CH4$VIF, decreasing = TRUE), ] #good!
vif_df_OMstock <- vif_df_OMstock[order(vif_df_OMstock$VIF, decreasing = TRUE), ] #good!
vif_df_OMstock_d <- vif_df_OMstock_d[order(vif_df_OMstock_d$VIF, decreasing = TRUE), ] #remove DailyMeanWaterTemp_C, good!
vif_df_OMflux <- vif_df_OMflux[order(vif_df_OMflux$VIF, decreasing = TRUE), ] #good!
vif_df_OMflux_d <- vif_df_OMflux_d[order(vif_df_OMflux_d$VIF, decreasing = TRUE), ] #remove DailyMeanWaterTemp_C, good!
vif_df_kcoarse <-  vif_df_kcoarse[order(vif_df_kcoarse$VIF, decreasing = TRUE), ] #remove discharge
vif_df_kcoarse_d <- vif_df_kcoarse_d[order(vif_df_kcoarse_d$VIF, decreasing = TRUE), ] #remove mean_velocity_m_s and mean_depth_cm 

####################################################################################

#### Run LMM for CO2 upstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
dat_lm_CO2_up <- dat %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C)


#scale the datam except for CO2
dat_lm_CO2_up_scaled <- dat_lm_CO2_up %>%
  select(-CO2_C_mg_m2_h) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the CO2 column back
dat_scaled <- cbind(dat_lm_CO2_up_scaled, CO2_C_mg_m2_h = dat_lm_CO2_up$CO2_C_mg_m2_h)

#start with a saturated model 
CO2_lmer <- lmer( log(CO2_C_mg_m2_h+16.3819) ~  water_pH + water_conductivity_us_cm + DO_. + canopy_cover_. +  wetted_width_m+ discharge_m3_s + mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + season*OM_stock_g_m2+ season*OM_flux_g_m2_s + sed_OM_percent + temp_C + DailyMeanWaterTemp_C + masl + bedrock + fine_substrate +  (1 + season | site), data = dat_scaled)

summary(CO2_lmer) 

plot(CO2_lmer) #log transforming gets rid of the fan shape
qqnorm(resid(CO2_lmer)) #maybe 2 outliers?

#Check outliers with Cook's Distance
cooksD <- cooks.distance(CO2_lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #60, 67

n <- nrow(dat_scaled_up)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_scaled[names_of_influential,]
dat_scaled_up_out <- dat_scaled_up %>% anti_join(outliers)

#Run again with outliers removed

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_scaled_up_out)[sapply(dat_scaled_up_out, is.numeric)]

#Saturated model 
CO2_lmer_sat <- lmer( log(CO2_C_mg_m2_h +16.3819) ~ 
                   water_pH + water_conductivity_us_cm + DO_. + canopy_cover_. +  wetted_width_m+ discharge_m3_s + mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + OM_stock_g_m2+ OM_flux_g_m2_s + sed_OM_percent + temp_C + DailyMeanWaterTemp_C + masl + bedrock + fine_substrate +
                   (1 + season | site), data = dat_scaled_up_out)

summary(CO2_lmer_sat) 

car::vif(CO2_lmer_sat) #all under 10

plot(CO2_lmer_sat) #Looks absolutely perfecto with log transformation and 2 outliers removed. log transforming gets rid of the fan shape
qqnorm(resid(CO2_lmer_sat)) #2 outliers removed

#Dredge is not working, too many variables I think, 
## try stepwise model selection ##

lmerTest::step(CO2_lmer_sat)

#Model found:
#log(CO2_C_mg_m2_h + 16.3819) ~ DO_. + mean_velocity_m_s + mean_depth_cm + OM_flux_g_m2_s + (1 + season | site)

# run the model 

CO2_lmer_up <- lmer( log(CO2_C_mg_m2_h +16.3819)~ DO_. + mean_velocity_m_s + mean_depth_cm + OM_flux_g_m2_s + (1 + season | site), data = dat_scaled_up_out)

summary(CO2_lmer_up) 
plot(CO2_lmer_up)
qqnorm(resid(CO2_lmer_up))
vif(CO2_lmer_up) #no issues, all around 1
r.squaredGLMM(CO2_lmer_up) #0.12 marginal, 0.79 conditional

#Run model based on predicitons

CO2_lmer_pred <- lmer( log(CO2_C_mg_m2_h +16.3819)~ poly(dist_to_source_km, 2)  +  mean_velocity_m_s + mean_depth_cm + OM_flux_g_m2_s + sed_OM_percent + fine_substrate + (1 + season | site), data = dat_scaled_up_out)

summary(CO2_lmer_pred)


#Maybe we also need to run a dredge on this bc quite a lot of variables
#Run dredge
options(na.action = "na.fail") # Required for dredge to run

dredge_CO2_lmer_pred<- dredge(CO2_lmer_pred, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

head(dredge_CO2_lmer_pred)

top.lmer_dredge_CO2_lmer_pred  <- get.models(dredge_CO2_lmer_pred , subset=delta <= 2)
mavg <- model.avg(top.lmer_dredge_CO2_lmer_pred)


#Run the dredge model averaged model 

CO2_lmer_pred <- lmer( log(CO2_C_mg_m2_h +16.3819)~ poly(dist_to_source_km, 2)  +  mean_velocity_m_s + fine_substrate + (1 + season | site), data = dat_scaled_up_out)

summary(CO2_lmer_pred)

r.squaredGLMM(CO2_lmer_pred) .32


#Simple LMM with dist to source

CO2up_lm <- lmer(log(CO2_C_mg_m2_h)~ dist_to_source_km + (1 | site), data=subset(dat, season=="Fall" & position_d=="downstream"))
summary(CO2up_lm)
plot(CO2up_lm)

x <- subset(dat, season=="Fall" & position_d=="downstream")
min(x$CO2_C_mg_m2_h)


####################################################################################

#### Run LMM for CO2 downstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
dat_lm_CO2_d <- dat %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -bedrock, -pebble_gravel, -DO_mg_L)


#scale the data, except for CO2
dat_lm_CO2_d_scaled <- dat_lm_CO2_d %>%
  select(-CO2_C_mg_m2_h) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the CO2 column back
dat_scaled_d <- cbind(dat_lm_CO2_d_scaled, CO2_C_mg_m2_h = dat_lm_CO2_d$CO2_C_mg_m2_h)

#start with a saturated model 
min(dat_scaled_d$CO2_C_mg_m2_h) #-38.37704

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_scaled_d)[sapply(dat_scaled_d, is.numeric)]

CO2_lmer_d_sat <- lmer(CO2_C_mg_m2_h ~ 
    water_pH +  water_conductivity_us_cm + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s +      mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s +     sed_OM_percent + DailyMeanWaterTemp_C + masl   
              +  (1 + season | site), data = dat_scaled_d)
summary(CO2_lmer_d_sat)

plot(CO2_lmer_d_sat) #look good, no need to log transform
qqnorm(resid(CO2_lmer_d_sat)) #maybe one outlier?
vif(CO2_lmer_d_sat) #remove fine substrate, at 18

#Check outliers with Cook's Distance
cooksD <- cooks.distance(CO2_lmer_d_sat) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #74

n <- nrow(dat_scaled_d)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_scaled_d[names_of_influential,]
dat_scaled_d_out <- dat_scaled_d %>% anti_join(outliers)

#Run again with outliers removed

CO2_lmer_d_sat2 <- lmer(CO2_C_mg_m2_h ~ 
                         water_pH +  water_conductivity_us_cm + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s +  mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + DailyMeanWaterTemp_C + masl        
                       +  (1 + season | site), data = dat_scaled_d_out)
summary(CO2_lmer_d_sat2)
plot(CO2_lmer_d_sat2) 
qqnorm(resid(CO2_lmer_d_sat2)) #improved with outlier removed


#run stepwise model selection
lmerTest::step(CO2_lmer_d_sat2)

#Model found:
#CO2_C_mg_m2_h ~ canopy_cover_. + discharge_m3_s + mean_velocity_m_s + k_dday_coarse + k_dday_fine + (1 + season | site)

#run the model

CO2_lmer_d <- lmer(CO2_C_mg_m2_h ~ canopy_cover_. + discharge_m3_s + mean_velocity_m_s + k_dday_coarse + k_dday_fine +  (1 + season | site) + (1 + season | site), data = dat_scaled_d)

summary(CO2_lmer_d)
plot(CO2_lmer_d) 
qqnorm(resid(CO2_lmer_d))

#Run model based on predicitons

CO2_lmer_pred <- lmer( log(CO2_C_mg_m2_h +16.3819)~ poly(dist_to_source_km, 2) + (1 + season | site), data = dat_scaled_d)

summary(CO2_lmer_pred)

##############################################################################

#### LMM CO2 all sites ####

#Use the dataset with the covarying variables removed, but add site and season back in 
dat_lm_CO2 <- dat %>% 
  select(-date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -bedrock, -pebble_gravel)


#scale the data, except for CO2
dat_lm_CO2_scaled <- dat_lm_CO2 %>%
  select(-CO2_C_mg_m2_h) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the CO2 column back
dat_scaled <- cbind(dat_lm_CO2_scaled, CO2_C_mg_m2_h = dat_lm_CO2$CO2_C_mg_m2_h)

#start with a saturated model 
min(dat_scaled$CO2_C_mg_m2_h) #-34.92515


#try model with percent openness
CO2_open <- lmer(log(CO2_C_mg_m2_h+35.92515) ~ dist_to_dam_km + dist_ds_dam_km + dist_to_source_km +percent_open_upstream + percent_open_downstream + (1 + season | site), data = dat_scaled) 
summary(CO2_open)
plot(CO2_open)
qqnorm(resid(CO2_open)) #outliers

#Check outliers with Cook's Distance
cooksD <- cooks.distance(CO2_open) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #134, 309, 310

n <- nrow(dat_scaled)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_scaled[names_of_influential,]
dat_scaled_out <- dat_scaled %>% anti_join(outliers)

#Run again with outliers removed

#try model with percent openness
CO2_open <- lmer(log(CO2_C_mg_m2_h+35.92515) ~ dist_to_dam_km  + dist_to_source_km + season*percent_open_upstream + season*percent_open_downstream + water_pH + DO_. +  (1 + season | site), data = dat_scaled_out) 
summary(CO2_open) #removed  dist_ds_dam
plot(CO2_open)
qqnorm(resid(CO2_open)) #improved
vif(CO2_open)

lmerTest::step(CO2_open)

options(na.action = "na.fail") # Required for dredge to run

dredge_CO2_open<- dredge(CO2_open, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

head(dredge_CO2_open)

top.lmerdredge_CO2_open  <- get.models(dredge_CO2_open , subset=delta <= 2)
mavg <- model.avg(top.lmerdredge_CO2_open)

CO2_open2 <- lmer(log(CO2_C_mg_m2_h + 35.92515) ~ DO_. + dist_to_source_km +  (1 + season | site), data = dat_scaled_out)
summary(CO2_open2)

#### Run LMMs by season ####

CO2_pos_lmm <- lmer(log(CO2_C_mg_m2_h) ~ dist_to_dam_km + (1 | site), data=subset(dat, season=="Fall"))

CO2_pos_lmm <- lmer(log(CO2_C_mg_m2_h) ~ position_d + (1 +season | site), data=dat)
summary(CO2_pos_lmm)
plot(CO2_pos_lmm)

CO2_pos_lm <- lm(log(CO2_C_mg_m2_h+14.22301) ~ dist_to_source_km, data=subset(dat, season=="Fall" & position_d=="upstream"))
summary(CO2_pos_lm)
plot(CO2_pos_lm)


x <- subset(dat, season=="Fall")
min(x$CO2_C_mg_m2_h)

emmeans(CO2_pos_lmm, pairwise ~ position_d) #post hoc test


CO2_season_lmm <- lmer(log(CO2_C_mg_m2_h+35.92515) ~ season + (1 +season | site), data=dat)
summary(CO2_season_lmm)
plot(CO2_season_lmm)

emmeans(CO2_season_lmm, pairwise ~ season) 



##############################################################################

#### Run LMM for CH4 upstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
dat_lm_CH4_up <- dat %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -masl, -position)


#scale the data except for CH4
dat_lm_CH4_up_scaled <- dat_lm_CH4_up %>%
  select(-CH4_C_mg_m2_h) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the CH4 column back
dat_scaled_up_CH4 <- cbind(dat_lm_CH4_up_scaled, CH4_C_mg_m2_h = dat_lm_CH4_up$CH4_C_mg_m2_h)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_scaled_up_CH4)[sapply(dat_scaled_up_CH4, is.numeric)]

#start with a saturated model 
min(dat_scaled_up_CH4$CH4_C_mg_m2_h)

CH4_lmer <- lmer((CH4_C_mg_m2_h^(1/3)) ~  #log(CH4_C_mg_m2_h + 1.02248368)
water_pH + water_conductivity_us_cm + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s +     sed_OM_percent + temp_C + DailyMeanWaterTemp_C + dist_to_source_km + bedrock + fine_substrate  +  (1 + season | site), data = dat_scaled_up_CH4)

summary(CH4_lmer) 

plot(CH4_lmer) #cube root transformation improves the fan shape the most (from log and square root)
qqnorm(resid(CH4_lmer)) #1 or 2 outliers

#Check outliers with Cook's Distance
cooksD <- cooks.distance(CH4_lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #173, 194

n <- nrow(dat_scaled_up_CH4)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_scaled_up_CH4[names_of_influential,]
dat_scaled_up_CH4_out <- dat_scaled_up_CH4 %>% anti_join(outliers)

#Run again with outliers removed

#Saturated model 
CH4_lmer_sat <- lmer((CH4_C_mg_m2_h^(1/3)) ~  #log(CH4_C_mg_m2_h + 1.02248368)
                       water_pH + water_conductivity_us_cm + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s +     sed_OM_percent + temp_C + DailyMeanWaterTemp_C + dist_to_source_km + bedrock + fine_substrate  +  (1 + season | site), data = dat_scaled_up_CH4_out)

summary(CH4_lmer_sat) 

plot(CH4_lmer_sat) #Improves the fan shape
qqnorm(resid(CH4_lmer_sat)) #2 outliers removed
car::vif(CH4_lmer_sat) #all under 10

#run stepwise model selection
lmerTest::step(CH4_lmer_sat)

#Model found:
#(CH4_C_mg_m2_h^(1/3)) ~ canopy_cover_. + mean_depth_cm + OM_stock_g_m2 + sed_OM_percent + temp_C + dist_to_source_km + (1 + season | site)

#run the model

CH4_lmer_up <- lmer( (CH4_C_mg_m2_h^(1/3)) ~ canopy_cover_. + mean_depth_cm + OM_stock_g_m2 + sed_OM_percent + temp_C + dist_to_source_km +  (1 + season | site) + (1 + season | site), data = dat_scaled_up_CH4_out)

summary(CH4_lmer_up)
plot(CH4_lmer_up) 
qqnorm(resid(CH4_lmer_up))


#Run model based on predictions
#Run model based on predicitons

CH4_lmer_pred <- lmer(  (CH4_C_mg_m2_h^(1/3)) ~ poly(dist_to_source_km, 2)  +  mean_velocity_m_s + mean_depth_cm + OM_flux_g_m2_s + sed_OM_percent + fine_substrate + (1 + season | site), data = dat_scaled_up_CH4_out)

summary(CH4_lmer_pred)


#Simple LMM with dist to source

CH4up_lm <- lmer((CH4_C_mg_m2_h^(1/3))~ dist_to_source_km + (1 | site), data=subset(dat, season=="Fall" & position_d=="downstream"))
summary(CH4up_lm)
plot(CH4up_lm)


####################################################################################

#### Run LMM for CH4 downstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
dat_lm_CH4_d <- dat %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -cobble, -boulder, -bedrock, -pebble_gravel, -DO_mg_L, -water_temp_C, -temp_C)

#scale the data, except for CH4
dat_lm_CH4_d_scaled <- dat_lm_CH4_d %>%
  select(-CH4_C_mg_m2_h) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the CH4 column back
dat_scaled_d_CH4 <- cbind(dat_lm_CH4_d_scaled, CH4_C_mg_m2_h = dat_lm_CH4_d$CH4_C_mg_m2_h)

#start with a saturated model 
min(dat_scaled_d_CH4$CH4_C_mg_m2_h) #-0.003544475

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_scaled_d_CH4)[sapply(dat_scaled_d_CH4, is.numeric)]

min(dat_scaled_d_CH4$CH4_C_mg_m2_h)

CH4_lmer_d_sat <- lmer((CH4_C_mg_m2_h^(1/3))  ~ water_pH + water_conductivity_us_cm + DO_. +  canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + DailyMeanWaterTemp_C + masl + (1 + season | site), data = dat_scaled_d_CH4)

summary(CH4_lmer_d_sat)

plot(CH4_lmer_d_sat) #log transforming improves residuals, but still not perfect, bu square root transformation is better 
qqnorm(resid(CH4_lmer_d_sat)) #maybe 3 outliers
vif(CH4_lmer_d_sat) #remove fine substrate *=(VIF of 12)


#Check outliers with Cook's Distance
cooksD <- cooks.distance(CH4_lmer_d_sat) #run the model with dat_sub first
influential <- cooksD[(cooksD > (8* mean(cooksD, na.rm = TRUE)))]
influential #22 , 74, 66

n <- nrow(dat_scaled_d_CH4)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 15/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_scaled_d_CH4[names_of_influential,]
dat_scaled_d_CH4_out <- dat_scaled_d_CH4 %>% anti_join(outliers)

#rerun with outliers removed
CH4_lmer_d <- lmer((CH4_C_mg_m2_h^(1/3))  ~ water_pH + water_conductivity_us_cm + DO_. +  canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm + k_dday_coarse + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + DailyMeanWaterTemp_C + masl + (1 + season | site), data = dat_scaled_d_CH4_out)

summary(CH4_lmer_d)
plot(CH4_lmer_d)
qqnorm(resid(CH4_lmer_d)) 

#Run model selection

#run stepwise model selection
lmerTest::step(CH4_lmer_d)

#Model found:
#(CH4_C_mg_m2_h^(1/3)) ~ DO_. + discharge_m3_s + mean_velocity_m_s + OM_stock_g_m2 + (1 + season | site)

#Run the model

CH4lmer_down <- lmer( (CH4_C_mg_m2_h^(1/3)) ~ DO_. + discharge_m3_s + mean_velocity_m_s + OM_stock_g_m2 + (1 + season | site),
data = dat_scaled_d_CH4_out)

summary(CH4lmer_down)
plot(CH4lmer_down)
qqnorm(resid(CH4lmer_down))
vif(CH4lmer_down)

#Run model based on predictions

CH4_lmer_pred <- lmer(  (CH4_C_mg_m2_h^(1/3)) ~ poly(dist_to_source_km, 2)  +  (1 + season | site), data = dat_scaled_d_CH4_out)

summary(CH4_lmer_pred)


##############################################################################

#### LMM CH4 all sites ####

#Use the dataset with the covarying variables removed, but add site and season back in 
dat_lm_CH4 <- dat %>% 
  select(-date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -flow_state, -cobble, -boulder, -bedrock, -pebble_gravel)


#scale the data, except for CH4
dat_lm_CH4_scaled <- dat_lm_CH4 %>%
  select(-CH4_C_mg_m2_h) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the CO2 column back
dat_scaled_CH4 <- cbind(dat_lm_CH4_scaled, CH4_C_mg_m2_h = dat_lm_CH4$CH4_C_mg_m2_h)

#start with a saturated model 
min(dat_scaled_CH4$CH4_C_mg_m2_h) #[1] -0.0198468


#try model with percent openness
CH4_open <- lmer(CH4_C_mg_m2_h ~ season*position_d + dist_to_dam_km + percent_open_upstream + percent_open_downstream  + fine_substrate + DailyMeanWaterTemp_C + (1 + season | site), data = dat_scaled_CH4) 
summary(CH4_open) 
plot(CH4_open) #improved a lot but still a bit fanned
qqnorm(resid(CH4_open))
vif(CH4_open)

#model selection
lmerTest::step(CH4_open)

options(na.action = "na.fail") # Required for dredge to run

dredge_CH4_open<- dredge(CH4_open, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

head(dredge_CH4_open)

top.lmer_dredge_CH4_open  <- get.models(dredge_CH4_open , subset=delta <= 2)
mavg <- model.avg(top.lmer_dredge_CH4_open)

#Run with variables from the top dredge model
CH4_open2 <- lmer( (CH4_C_mg_m2_h^(1/3))  ~ fine_substrate + position_d*season + percent_open_upstream +(1 + season | site), data = dat_scaled_CH4) 
summary(CH4_open2) 
plot(CH4_open2)

## Run with env vars

CH4_env <- lmer( (CH4_C_mg_m2_h^(1/3))  ~ water_pH + water_conductivity_us_cm + DO_mg_L + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + DailyMeanWaterTemp_C + temp_C + substrate_complexity + fine_substrate + mean_depth_cm +  masl + dist_to_source_km + OM_stock_g_m2 + (1 + season | site), data = dat_scaled_CH4) 

summary(CH4_env)
plot(CH4_env)
qqnorm(resid(CH4_env))
vif(CH4_env) #all good

#model selection
lmerTest::step(CH4_env)

#Model found:
#(CH4_C_mg_m2_h^(1/3)) ~ water_pH + DO_. + canopy_cover_. + mean_velocity_m_s + mean_depth_cm + OM_stock_g_m2 + (1 + season | site)

#run the top model


CH4_env2 <- lmer(  (CH4_C_mg_m2_h^(1/3)) ~ water_pH + DO_. + canopy_cover_. + mean_velocity_m_s + mean_depth_cm + OM_stock_g_m2 + (1 + season | site), data = dat_scaled_CH4) 
summary(CH4_env2)

#### Run LMMs by season ####


CH4_pos_lmm <- lm(CH4_C_mg_m2_h ~ dist_to_dam_km, data=dat)
summary(CH4_pos_lmm)
plot(CH4_pos_lmm)


CH4_pos_lmm <- lmer(CH4_C_mg_m2_h^(1/3) ~ dist_to_dam_km + (1 | site), data=subset(dat, season=="Summer"))
summary(CH4_pos_lmm)
plot(CH4_pos_lmm)
#to interpret coefficient need to untransform:
-0.05986^3 #-0.0002144915
(-0.05986*10)^3 #-0.2144915

CH4_position_lmm <- lmer(CH4_C_mg_m2_h^(1/3) ~ position_d + (1 | site), data=subset(dat, season=="Summer"))
summary(CH4_position_lmm)
plot(CH4_position_lmm)

emmeans(CH4_position_lmm, pairwise ~ position_d) #post hoc test


CH4_season_lmm <- lmer(CH4_C_mg_m2_h^(1/3) ~ season + (1 +season | site), data=dat)
summary(CH4_season_lmm)
plot(CH4_season_lmm)

emmeans(CH4_season_lmm, pairwise ~ season) 

################################################################################

#### Run LMM for CO2 equivalents ############

min(dat$CO2e[dat$season == "Fall"]) #Spring -13.25686, Summer -31.12941
min(dat$CO2e[dat$CO2e >= 0]) #0.1447593 
0.1447593+31.12941 #31.27417

dat_CO2e_Fall <- subset(dat, season=="Fall") #remove 1 outlier in Fall
dat_CO2e_Fall_out <- subset(dat_CO2e_Fall, CO2e<=1000)

CO2e_source_lmm <- lmer(log(CO2e) ~ dist_to_source_km + (1 | site), data=subset(dat, season=="Fall"))
summary(CO2e_source_lmm) #no sig relationship, p = 0.09 in the fall...
plot(CO2e_source_lmm)

CO2e_source_lmm_all <- lmer(log(CO2e) ~ dist_to_source_km + (1 +season | site), data=dat)
summary(CO2e_source_lmm_all) #not sig
plot(CO2e_source_lmm_all)

CO2e_source_lmm_season <- lmer(log(CO2e) ~ season + (1 +season | site), data=dat)
summary(CO2e_source_lmm_season) 
plot(CO2e_source_lmm_season)
emmeans(CO2e_source_lmm_season, pairwise ~ season) #no effect of season

CO2e_dam_lmm <- lmer(log(CO2e+32.12941) ~ dist_to_dam_km + (1 | site), data=subset(dat_CO2e_Summer_out, season=="Summer"))
summary(CO2e_dam_lmm) #no sig relationships
plot(CO2e_dam_lmm)

CO2e_pos_lmm <- lmer(log(CO2e) ~ position_d + (1 | site), data=subset(dat, season=="Fall")) #dat_CO2e_Summer_out
summary(CO2e_pos_lmm) #in the fall, reservoir sig lower than upstream
plot(CO2e_pos_lmm)
emmeans(CO2e_pos_lmm, pairwise ~ position_d) 

#1 outlier in summer - check  with Cook's Distance
dat_CO2e_Summer <- subset(dat, season=="Summer")
cooksD <- cooks.distance(CO2e_pos_lmm) #run the model with dat_sub first
influential <- cooksD[(cooksD > (8* mean(cooksD, na.rm = TRUE)))]
influential #15

n <- nrow(dat_CO2e_Summer)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 15/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_CO2e_Summer[names_of_influential,]
dat_CO2e_Summer_out <- dat_CO2e_Summer %>% anti_join(outliers)






emmeans(CH4_pos_lmm, pairwise ~ position_d) #post hoc test

CH4_season_lmm <- lmer(CH4_C_mg_m2_h^(1/3) ~ season + (1 +season | site), data=dat)
summary(CH4_season_lmm)
plot(CH4_season_lmm)




####################################################################################

#### Run LMM for OM stock upstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
#consider the ecology, water physicochemistry are not likely drivers, select only plausible drivers
#note we can use the site averages here 
dat_lm_OMstock <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -campaign, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -masl, -water_pH, -water_conductivity_us_cm, -sed_OM_percent, -DailyMeanWaterTemp_C, -DO_., -bedrock, -fine_substrate, -k_dday_fine, -k_dday_coarse)

#scale the data, except for CH4
dat_lm_OMstock_scaled <- dat_lm_OMstock %>%
  select(-OM_stock_g_m2) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_OMstock_scaled <- cbind(dat_lm_OMstock_scaled, OM_stock_g_m2 = dat_lm_OMstock$OM_stock_g_m2)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_OMstock_scaled)[sapply(dat_lm_OMstock_scaled, is.numeric)]


OMstock_lmer_sat <- lmer( OM_stock_g_m2 ~ season + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm +   OM_flux_g_m2_s  + temp_C + dist_to_source_km +  dist_ds_dam_km +   (1 | site), data = dat_lm_OMstock_scaled)

summary(OMstock_lmer_sat)

plot(CH4_lmer_d_sat) # log and cube root transformation does not improve residuals
qqnorm(resid(CH4_lmer_d_sat)) # 2 or 3 outliers
vif(OMstock_lmer_sat) #all <10

#Check outliers with Cook's Distance
cooksD <- cooks.distance(OMstock_lmer_sat) #run the model with dat_sub first
influential <- cooksD[(cooksD > (8* mean(cooksD, na.rm = TRUE)))]
influential #60, 77

n <- nrow(dat_scaled_d_CH4)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 15/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_lm_OMstock_scaled[names_of_influential,]
dat_lm_OMstock_scaled_out <- dat_lm_OMstock_scaled %>% anti_join(outliers)

#Rerun with outliers removed

OMstock_lmer <- lmer( OM_stock_g_m2 ~ season + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm +   OM_flux_g_m2_s  + temp_C + dist_to_source_km +  dist_ds_dam_km +   (1 | site), data = dat_lm_OMstock_scaled_out)

summary(OMstock_lmer)

plot(OMstock_lmer) # better
qqnorm(resid(OMstock_lmer)) # better
vif(OMstock_lmer) #all <10

#run stepwise model selection
lmerTest::step(OMstock_lmer)

#Model found:
#OM_stock_g_m2 ~ wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm + OM_flux_g_m2_s + (1 | site)

#Run the model
OMstock_lmer2 <- lmer( OM_stock_g_m2 ~ discharge_m3_s + mean_velocity_m_s + OM_flux_g_m2_s + dist_to_source_km + (1 | site), data = dat_lm_OMstock_scaled_out)

summary(OMstock_lmer2)
plot(OMstock_lmer2)
qqnorm(resid(OMstock_lmer2)) 
vif(OMstock_lmer2)


#Fit several regressions to determine highest R2, then see it it matches predicted patterns from Serial Disc. Concept

#fit linear model 

OMstock_lm <- lm(log(OM_stock_g_m2) ~ dist_to_source_km, data = dat_lm_OMstock_scaled)
summary(OMstock_lm) #8.4% variance 
plot(OMstock_lm)
qqnorm(resid(OMstock_lm)) 

#fit linear mixed effects model 
OMstock_lmer<- lmer(log(OM_stock_g_m2) ~ dist_to_source_km  + (1 | site), data = dat_lm_OMstock_scaled)
summary(OMstock_lmer)
plot(OMstock_lmer)
qqnorm(resid(OMstock_lmer)) 
r.squaredGLMM(OMstock_lmer) #.5% marginal, 47 conditional
AIC(OMstock_lmer) #320.8983

OMstock_lmer_poly<- lmer(log(OM_stock_g_m2) ~ poly(dist_to_source_km,2)*season  + (1 | site), data = dat_lm_OMstock_scaled)
summary(OMstock_lmer_poly)
plot(OMstock_lmer_poly)
qqnorm(resid(OMstock_lmer_poly)) 
r.squaredGLMM(OMstock_lmer_poly) #5.4% marginal, 48 conditional
AIC(OMstock_lmer_poly) #132.4304 thus it is a better fit to the data 


#try GAM model (then GLMM)

gam_model <- gam(log(OM_stock_g_m2) ~  s(dist_to_source_km) , data = dat_lm_OMstock_scaled, method="REML")
summary(gam_model)
anova.Gam(gam_model)
par(mfrow = c(1, 2))
plot(gam_model, all.terms = TRUE)
r.squaredGLMM(gam_model) 
AIC(OMstock_lm, gam_model)

data_plot <- ggplot(data = dat_lm_OMstock_scaled, aes(y = log(OM_stock_g_m2), x = dist_to_source_km)) +
  geom_point() + geom_line(aes(y = fitted(OMstock_lm)), colour = "red",
                           size = 1.2) + theme_bw()
data_plot

data_plot <- data_plot + geom_line(aes(y = fitted(gam_model)),
                                   colour = "blue", size = 1.2)
data_plot

#Try GAMM model where you can specify random effects: Random intercepts adjust the height of other model terms with a constant value: s(fac, bs = "re") 

gamm_stock <- gam(log(OM_stock_g_m2) ~ s(dist_to_source_km, by=season) + s(site, bs = "re"), data = dat_lm_OMstock_scaled, method = "REML")
summary(gamm_stock)
par(mfrow = c(2, 2))
gam.check(gamm_stock) #check assumptions for normality, good!
plot(gamm_stock)
edf(gamm_stock)
plot(ggeffects::ggpredict(gamm_stock), facets = TRUE)



####################################################################################

#### Run LMM for OM stock downstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
#consider the ecology, water physicochemistry are not likely drivers, select only plausible drivers
dat_lm_OMstock_d <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date,   -campaign,  -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -masl, -water_pH, -water_conductivity_us_cm, -sed_OM_percent, -DailyMeanWaterTemp_C, -DO_., -bedrock, -fine_substrate, -bedrock, -k_dday_coarse, -k_dday_fine)

#scale the data, except for CH4
dat_lm_OMstock_d_scaled <- dat_lm_OMstock_d %>%
  select(-OM_stock_g_m2) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_OMstock_d_scaled <- cbind(dat_lm_OMstock_d_scaled, OM_stock_g_m2 = dat_lm_OMstock_d$OM_stock_g_m2)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_OMstock_d_scaled)[sapply(dat_lm_OMstock_d_scaled, is.numeric)]


OMstock_d_lmer_sat <- lmer( log(OM_stock_g_m2) ~ canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm +   OM_flux_g_m2_s  + temp_C +  dist_ds_dam_km +   (1 | site), data = dat_lm_OMstock_d_scaled)

summary(OMstock_d_lmer_sat)

plot(OMstock_d_lmer_sat) # OK
qqnorm(resid(OMstock_d_lmer_sat)) # OK
vif(OMstock_d_lmer_sat) # remove distance to source and season


#run stepwise model selection
lmerTest::step(OMstock_d_lmer_sat)

#Model found:
 # log(OM_stock_g_m2) ~ mean_velocity_m_s

#Run the selected model

OMstock_d_lmer <- lmer( log(OM_stock_g_m2) ~  mean_velocity_m_s +  (1 | site), data = dat_lm_OMstock_d_scaled)

summary(OMstock_d_lmer)
plot(OMstock_d_lmer) # OK
qqnorm(resid(OMstock_d_lmer)) #OK 

#Test which type of relationship best fits the data

#fit linear mixed effects model 
OMstock_lmer<- lmer(log(OM_stock_g_m2) ~ dist_to_source_km  + (1 | site), data = dat_lm_OMstock_d_scaled)
summary(OMstock_lmer)
plot(OMstock_lmer)
qqnorm(resid(OMstock_lmer)) 
r.squaredGLMM(OMstock_lmer) #.03
AIC(OMstock_lmer) #46.47102

#Fit quadratic polynomial relationship
OMstock_lmer_poly<- lmer(log(OM_stock_g_m2) ~ season*poly(dist_to_source_km,2)  + (1 | site), data = dat_lm_OMstock_d_scaled)
summary(OMstock_lmer_poly)
plot(OMstock_lmer_poly)
qqnorm(resid(OMstock_lmer_poly)) 
r.squaredGLMM(OMstock_lmer_poly) #.06
AIC(OMstock_lmer_poly) #43.63861 thus it is a better fit to the data 

##Try GAMM model 

gamm_stock_d <- gam(log(OM_stock_g_m2) ~ s(dist_to_source_km, by=season) + s(site, bs = "re"), data = dat_lm_OMstock_d_scaled, method = "REML") #model has too few data points, so it won't run




###############################################################################
#### Run LMM for OM stock all sites ####

#include all variables for now, remove via VIF later
#exclude ecologically implausible predictors, like water chemistry
dat_lm_OMstock_all <- dat_means %>% 
  select(-date,   -campaign, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -water_pH, -water_conductivity_us_cm, -sed_OM_percent, -fine_substrate, -k_dday_coarse, -k_dday_fine, -bedrock, -mean_velocity_m_s)

#scale the data, except for CH4
dat_lm_OMstock_all_scaled <- dat_lm_OMstock_all %>%
  select(-OM_stock_g_m2) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_OMstock_all_scaled <- cbind(dat_lm_OMstock_all_scaled, OM_stock_g_m2 = dat_lm_OMstock_all$OM_stock_g_m2)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_OMstock_all_scaled)[sapply(dat_lm_OMstock_all_scaled, is.numeric)]

OMstock_lmer_sat <- lmer(log(OM_stock_g_m2) ~ season +  water_temp_C + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s +  mean_depth_cm + OM_flux_g_m2_s + temp_C + position*masl  +  (1 | site), data = dat_lm_OMstock_all_scaled)

summary(OMstock_lmer_sat)

plot(OMstock_lmer_sat) # log transforming improves residuals
qqnorm(resid(OMstock_lmer_sat)) # OK
vif(OMstock_lmer_sat) # remove DO_mg_L, dist_ds_dam_km, DailyMeanWaterTemp_C, dist_to_source_km

#Run stepwise model selection
lmerTest::step(OMstock_lmer_sat)

#Model found:
#log(OM_stock_g_m2) ~ season + wetted_width_m + mean_depth_cm + temp_C + (1 | site)

#Run top model

OMstock_lmer <- lmer(log(OM_stock_g_m2) ~ season + wetted_width_m + mean_depth_cm + temp_C +  (1 | site), data = dat_lm_OMstock_all_scaled)

summary(OMstock_lmer)
plot(OMstock_lmer) # 
qqnorm(resid(OMstock_lmer)) # 
vif(OMstock_lmer) #

# run model with just distance to dam

OMstock_dam <- lmer(log(OM_stock_g_m2) ~   position*poly(dist_to_dam_km, 2) + (1 | site), data = dat_lm_OMstock_all_scaled)
summary(OMstock_dam)
plot(OMstock_dam) # log transforming improves residuals
qqnorm(resid(OMstock_dam)) # 
vif(OMstock_dam)

lmerTest::step(OMstock_dam)

#Model based on predictions and top correlations
OMstock_open <- lmer(log(OM_stock_g_m2) ~  canopy_cover_. + season*percent_open_upstream + season*percent_open_downstream + (1 | site), data = dat_lm_OMstock_all_scaled) 
summary(OMstock_open)

#Maybe we also need to run a dredge on this bc quite a lot of variables
#Run dredge
options(na.action = "na.fail") # Required for dredge to run

dredge_OMstock_open<- dredge(OMstock_open, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

head(dredge_OMstock_open)

top.lmer_OMstock_open  <- get.models(dredge_OMstock_open , subset=delta <= 2)
mavg <- model.avg(top.lmer_OMstock_open)

#try model with percent openness
OMstock_open <- lmer(log(OM_stock_g_m2) ~   season*percent_open_upstream + season*percent_open_downstream + (1 | site), data = dat_lm_OMstock_all_scaled) #via dredge, remove canopy_cover_.
summary(OMstock_open)
plot(OMstock_open)


#Try GAMM model 

gamm_stock_all <- gam(log(OM_stock_g_m2) ~ s(dist_to_source_km, by=season) + s(site, bs = "re"), data = dat_lm_OMstock_all_scaled, method = "REML") 
summary(gamm_stock_all)
plot(gamm_stock_all)
plot(ggeffects::ggpredict(gamm_stock_all), facets = TRUE)


gam_stock <- gam(log(OM_stock_g_m2) ~ s(dist_to_source_km, bs = "cr"), data = dat_lm_OMstock_all_scaled)
summary(gam_stock)
plot(gam_stock)

#Spatial GAM
spatial_gam <- gam(OM_stock_g_m2 ~ s(lon, lat, bs = 'gp', k = 100, m = 2), data = dat_lm_OMstock_all_scaled)
summary(spatial_gam)


#### Run LMs by season ####

OMflux_lm <- lm(log(OM_stock_g_m2) ~ dist_to_source_km, data=subset(dat_means, season=="Summer" & position_d=="upstream"))
summary(OMflux_lm)
plot(OMflux_lm)

OMflux_lm <- lm(log(OM_g_m3) ~ dist_to_source_km, data=subset(dat_means, season=="Fall" & position_d=="downstream"))
summary(OMflux_lm)
plot(OMflux_lm)

OMflux_lm <- lm(log(OM_g_m3) ~ dist_to_source_km, data=subset(dat_means, season=="Fall" & position_d!="reservoir"))
summary(OMflux_lm)
plot(OMflux_lm)

OMflux_pos_lm <- lm(log(OM_g_m3) ~ position_d, data=subset(dat_means, season=="Fall"))
summary(OMflux_pos_lm)
plot(OMflux_pos_lm)

emmeans(OMflux_pos_lm, pairwise ~ position_d) #post hoc test

OMflux_season_lm <- lm(log(OM_stock_g_m2) ~ season, data=dat_means)
summary(OMflux_season_lm)
plot(OMflux_season_lm)

emmeans(OMflux_season_lm, pairwise ~ season) 

OMtransport_season_lm <- lm(log(OM_g_m3) ~ season, data=dat_means)
summary(OMtransport_season_lm)
plot(OMtransport_season_lm)

emmeans(OMtransport_season_lm, pairwise ~ season) 

###############################################################################
#### Run LMM for OM flux all sites ####

#include all variables for now, remove via VIF later
#exclude ecologically implausible predictors, like water chemistry
dat_lm_OMflux_all <- dat_means %>% 
  select(-date,   -campaign,  -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -water_pH, -water_conductivity_us_cm, -sed_OM_percent, -fine_substrate, -k_dday_coarse, -k_dday_fine, -bedrock, -mean_velocity_m_s)

#scale the data, except for CH4
dat_lm_OMflux_all_scaled <- dat_lm_OMflux_all %>%
  select(-OM_flux_g_m2_s) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_OMflux_all_scaled <- cbind(dat_lm_OMflux_all_scaled, OM_flux_g_m2_s = dat_lm_OMflux_all$OM_flux_g_m2_s)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_OMflux_all_scaled)[sapply(dat_lm_OMflux_all_scaled, is.numeric)]

OMflux_lmer_sat <- lmer(log(OM_flux_g_m2_s) ~ season + water_temp_C + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_depth_cm + OM_stock_g_m2 + temp_C +  position*masl + 
                          (1 | site), data = dat_lm_OMflux_all_scaled)

summary(OMflux_lmer_sat)

plot(OMflux_lmer_sat) # pattern improved by log transformation
qqnorm(resid(OMflux_lmer_sat)) # OK
vif(OMflux_lmer_sat) #remove DO_mg_L, dist_ds_dam_km, DailyMeanWaterTemp_C                       

#Run stepwise model selection
lmerTest::step(OMflux_lmer_sat)

#Model found:
#log(OM_flux_g_m2_s) ~ season + DO_. + mean_depth_cm + position

#run top model

OMflux_lmer <- lmer(log(OM_flux_g_m2_s) ~ season + DO_. + mean_depth_cm + position +
                          (1 | site), data = dat_lm_OMflux_all_scaled)

summary(OMflux_lmer)


#### Run LMM for OM flux upstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
#consider the ecology, water physicochemistry are not likely drivers, select only plausible drivers
dat_lm_OMflux<- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select(-date, -campaign, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -masl, -water_pH, -water_conductivity_us_cm, -sed_OM_percent, -DailyMeanWaterTemp_C, -DO_., -bedrock, -fine_substrate, -k_dday_fine, -k_dday_coarse)

#scale the data, except for CH4
dat_lm_OMflux_scaled <- dat_lm_OMflux %>%
  select(-OM_flux_g_m2_s) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_OMflux_scaled <- cbind(dat_lm_OMflux_scaled, OM_flux_g_m2_s = dat_lm_OMflux$OM_flux_g_m2_s)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_OMflux_scaled)[sapply(dat_lm_OMflux_scaled, is.numeric)]


OMflux_lmer_sat <- lmer(log(OM_flux_g_m2_s) ~ season + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm +   OM_stock_g_m2 + temp_C +   dist_ds_dam_km +   (1 | site), data = dat_lm_OMflux_scaled)

summary(OMflux_lmer_sat)

plot(OMflux_lmer_sat) # pattern improved by log transformation
qqnorm(resid(OMflux_lmer_sat)) # good
vif(OMflux_lmer_sat) # remove dist_to_source_km 

#run stepwise model selection
lmerTest::step(OMflux_lmer_sat)

#Model found:
#log(OM_flux_g_m2_s) ~ season

#Run top model
OMflux_lmer <- lmer(log(OM_flux_g_m2_s) ~  season +  (1 | site), data = dat_lm_OMflux_scaled)

summary(OMflux_lmer)
plot(OMflux_lmer) 
qqnorm(resid(OMflux_lmer)) 

OMflux_season_lm <- lm(log(OM_flux_g_m2_s) ~ season, data=dat_means)
summary(OMflux_season_lm)
plot(OMflux_season_lm)

emmeans(OMflux_season_lm, pairwise ~ season) 

## Run LMs by season ##

OMflux_lm <- lm(log(OM_flux_g_m2_s) ~ dist_to_source_km, data=subset(dat_means, season=="Fall" & position_d=="upstream"))
summary(OMflux_lm)
plot(OMflux_lm)


OMflux_lm <- lm(log(OM_flux_g_m2_s) ~ position_d, data=subset(dat_means, season=="Fall"))
summary(OMflux_lm)
plot(OMflux_lm)

emmeans(OMflux_lm, pairwise ~ position_d) #post hoc test




###############################################################################

#### Run LMM for OM flux downstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
#consider the ecology, water physicochemistry are not likely drivers, select only plausible drivers
dat_lm_OMflux_d <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select(-date, -campaign, -lat, -lon, -position, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -masl, -water_pH, -water_conductivity_us_cm, -sed_OM_percent, -DailyMeanWaterTemp_C, -DO_., -bedrock, -fine_substrate, -k_dday_fine, -k_dday_coarse)

#scale the data, except for CH4
dat_lm_OMflux_d_scaled <- dat_lm_OMflux_d %>%
  select(-OM_flux_g_m2_s) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_OMflux_d_scaled <- cbind(dat_lm_OMflux_d_scaled, OM_flux_g_m2_s = dat_lm_OMflux_d$OM_flux_g_m2_s)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_OMflux_d_scaled)[sapply(dat_lm_OMflux_d_scaled, is.numeric)]


OMflux_d_lmer_sat <- lmer(log(OM_flux_g_m2_s) ~ canopy_cover_. + wetted_width_m + discharge_m3_s +  mean_depth_cm +   OM_stock_g_m2 + temp_C +   dist_ds_dam_km +   (1 | site), data = dat_lm_OMflux_d_scaled)

summary(OMflux_d_lmer_sat)

plot(OMflux_d_lmer_sat) # pattern improved by log transformation
qqnorm(resid(OMflux_d_lmer_sat)) # good
vif(OMflux_d_lmer_sat) # remove season, remove velocity

#run stepwise model selection
lmerTest::step(OMflux_d_lmer_sat)

#Model found:
#log(OM_flux_g_m2_s) ~ canopy_cover_.

#Run model 

OMflux_d_lmer <- lmer(log(OM_flux_g_m2_s) ~ canopy_cover_. + (1 | site), data = dat_lm_OMflux_d_scaled)

summary(OMflux_d_lmer)

#### Run LMM for k coarse upstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
#consider the ecology, water physicochemistry are not likely drivers, select only plausible drivers
dat_lm_kcoarse<- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select( -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -DO_mg_L, -water_temp_C,  -discharge_m3_s, -bedrock, -DailyMeanWaterTemp_C)

#scale the data, except for CH4
dat_lm_kcoarse_scaled <- dat_lm_kcoarse %>%
  select(-k_dday_coarse) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_kcoarse_scaled <- cbind(dat_lm_kcoarse_scaled, k_dday_coarse = dat_lm_kcoarse$k_dday_coarse)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_kcoarse_scaled)[sapply(dat_lm_kcoarse_scaled, is.numeric)]


kcoarse_lmer_sat <- lmer(log(k_dday_coarse) ~  water_pH + water_conductivity_us_cm + DO_. + canopy_cover_. + wetted_width_m + mean_velocity_m_s + mean_depth_cm + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + temp_C + masl + fine_substrate +  (1 | site ), data = dat_lm_kcoarse_scaled)

summary(kcoarse_lmer_sat)

plot(kcoarse_lmer_sat) # slight fan shape, improved by log transformation
qqnorm(resid(kcoarse_lmer_sat)) # good
vif(kcoarse_lmer_sat) # 

lmerTest::step(kcoarse_lmer_sat)

#Model found:
#k_dday_coarse ~ water_conductivity_us_cm 

#Rerun model 

kcoarse_lmer <- lmer(log(k_dday_coarse) ~  water_conductivity_us_cm + (1 | site ), data = dat_lm_kcoarse_scaled)
summary(kcoarse_lmer)

plot(kcoarse_lmer) # slight fan shape, improved by log transformation
qqnorm(resid(kcoarse_lmer)) # OK

#Test model fit linear vs quadratic polynomial

#fit linear mixed effects model 
kcoarse_lmer<- lmer(log(k_dday_coarse)  ~ dist_to_source_km  + (1 | site), data = dat_lm_kcoarse_scaled)
summary(kcoarse_lmer)
plot(kcoarse_lmer)
qqnorm(resid(kcoarse_lmer)) 
r.squaredGLMM(kcoarse_lmer) #.07
AIC(kcoarse_lmer) #95.04504

kcoarse_lmer_poly<- lmer(log(k_dday_coarse) ~ poly(dist_to_source_km,2)  + (1 | site), data = dat_lm_kcoarse_scaled)
summary(kcoarse_lmer_poly)
plot(kcoarse_lmer_poly)
qqnorm(resid(kcoarse_lmer_poly)) 
r.squaredGLMM(kcoarse_lmer_poly) #.09
AIC(kcoarse_lmer_poly) #91.43



###############################################################################
#### Run LMM for k coarse downstream ####

#Use the dataset with the covarying variables removed, but add site and season back in 
#consider the ecology, water physicochemistry are not likely drivers, select only plausible drivers
dat_lm_kcoarse_d <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select( -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position, -DO_mg_L, -water_temp_C, -discharge_m3_s, -bedrock)

#scale the data, except for CH4
dat_lm_kcoarse_d_scaled <- dat_lm_kcoarse_d %>%
  select(-k_dday_coarse) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_kcoarse_d_scaled <- cbind(dat_lm_kcoarse_d_scaled, k_dday_coarse = dat_lm_kcoarse_d$k_dday_coarse)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_kcoarse_d_scaled)[sapply(dat_lm_kcoarse_d_scaled, is.numeric)]

kcoarse_d_lmer_sat <- lmer(log(k_dday_coarse) ~  water_pH + water_conductivity_us_cm + DO_. + canopy_cover_. + wetted_width_m + mean_velocity_m_s + k_dday_fine + OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + temp_C + masl +  (1 | site ), data = dat_lm_kcoarse_d_scaled)

summary(kcoarse_d_lmer_sat)

plot(kcoarse_d_lmer_sat) #ok
qqnorm(resid(kcoarse_d_lmer_sat)) #ok 
vif(kcoarse_d_lmer_sat) #remove DailyMeanWaterTemp_C and fine_substrate and mean_depth_cm

lmerTest::step(kcoarse_d_lmer_sat)
# no variables chosen

#Test model relationships linear vs polynomial

#fit linear mixed effects model 
kcoarse_lmer_d<- lmer(log(k_dday_coarse)  ~ dist_to_source_km  + (1 | site), data = dat_lm_kcoarse_d_scaled)
summary(kcoarse_lmer_d)
plot(kcoarse_lmer_d)
qqnorm(resid(kcoarse_lmer_d)) 
r.squaredGLMM(kcoarse_lmer_d) #.004
AIC(kcoarse_lmer_d) #29

kcoarse_lmer_poly<- lmer(log(k_dday_coarse) ~ poly(dist_to_source_km,2)  + (1 | site), data = dat_lm_kcoarse_d_scaled)
summary(kcoarse_lmer_poly)
plot(kcoarse_lmer_poly)
qqnorm(resid(kcoarse_lmer_poly)) 
r.squaredGLMM(kcoarse_lmer_poly) #.03
AIC(kcoarse_lmer_poly) #27.85


##############################################################################

#### Use all data with k coarse LMM ####

dat_lm_kcoarse_all <- dat_means %>% 
  select( -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel) #exclude DailyMeanWaterTemp_C because it's used in the calculation and water_temp_C point measure is not relevant

#scale the data, except for CH4
dat_lm_kcoarse_all_scaled <- dat_lm_kcoarse_all %>%
  select(-k_dday_coarse) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_kcoarse_all_scaled <- cbind(dat_lm_kcoarse_all_scaled, k_dday_coarse = dat_lm_kcoarse_all$k_dday_coarse)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_kcoarse_all_scaled)[sapply(dat_lm_kcoarse_all_scaled, is.numeric)]

kcoarse_all_lmer_sat <- lmer(log(k_dday_coarse) ~ season +  water_pH + water_conductivity_us_cm + DO_. + canopy_cover_. +  discharge_m3_s + mean_velocity_m_s + mean_depth_cm + position*OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + temp_C + masl + bedrock + fine_substrate + (1 | site ) , data = dat_lm_kcoarse_all_scaled)

summary(kcoarse_all_lmer_sat)
plot(kcoarse_all_lmer_sat) # 
qqnorm(resid(kcoarse_all_lmer_sat)) # 
vif(kcoarse_all_lmer_sat) # Remove dist_to_source_km and dist_ds_dam_km and wetted_width_m and DO_mg_L

lmerTest::step(kcoarse_all_lmer)

#Model found:
#log(k_dday_coarse) ~ season + water_conductivity_us_cm + DO_. + discharge_m3_s + OM_stock_g_m2 + OM_flux_g_m2_s + temp_C + position + dist_to_source_km + fine_substrate + position:dist_to_source_km

kcoarse_all_lmer <- lmer(log(k_dday_coarse) ~ season + water_conductivity_us_cm + DO_. + discharge_m3_s + OM_stock_g_m2 + OM_flux_g_m2_s + temp_C + position + dist_to_source_km + fine_substrate + position:dist_to_source_km + (1 | site ) , data = dat_lm_kcoarse_all_scaled)

summary(kcoarse_all_lmer)
plot(kcoarse_all_lmer) # OK
qqnorm(resid(kcoarse_all_lmer)) #  OK
vif(kcoarse_all_lmer) # OK


#Maybe we also need to run a dredge on this bc quite a lot of variables
#Run dredge
options(na.action = "na.fail") # Required for dredge to run

lmer_kcoarse_dredge<- dredge(kcoarse_all_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(lmer_kcoarse_dredge)  #how many models were run 1280
head(lmer_kcoarse_dredge)

top.lmer_kcoarse_dredge  <- get.models(lmer_kcoarse_dredge , subset=delta <= 4) #chose 4 because null model otherwise 

#Formula: log(k_dday_coarse) ~ OM_stock_g_m2 + season + water_conductivity_us_cm +      (1 | site)

#run the model 

kcoarse_dredged_lmer <- lmer(log(k_dday_coarse) ~ OM_stock_g_m2 + season + water_conductivity_us_cm + (1 | site ) , data = dat_lm_kcoarse_all_scaled)

summary(kcoarse_dredged_lmer)


## Run model from hypothesized variables

kcoarse__lmer <- lmer(log(k_dday_coarse) ~ poly(dist_to_source_km,2)*position + mean_velocity_m_s + substrate_complexity + OM_stock_g_m2 + (1 | site ) , data = dat_lm_kcoarse_all_scaled)

summary(kcoarse__lmer)
plot(kcoarse__lmer) # 1 outlier
qqnorm(resid(kcoarse__lmer)) #  OK
vif(kcoarse__lmer) # OK

#Check outliers with Cook's Distance
cooksD <- cooks.distance(kcoarse__lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (5* mean(cooksD, na.rm = TRUE)))]
influential #29

n <- nrow(dat_lm_kcoarse_all_scaled)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 15/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_lm_kcoarse_all_scaled[names_of_influential,]
dat_lm_kcoarse_all_scaled_out <- dat_lm_kcoarse_all_scaled %>% anti_join(outliers)

#Rerun model

kcoarse__lmer <- lmer(log(k_dday_coarse) ~ poly(dist_to_source_km,2)*position + position*percent_open_upstream + position_d*mean_velocity_m_s + position_d*substrate_complexity + position_d*OM_stock_g_m2 + position_d*DailyMeanWaterTemp_C + (1 | site ) , data = dat_lm_kcoarse_all_scaled_out)

summary(kcoarse__lmer)
plot(kcoarse__lmer) # improves after removing 1 outlier
qqnorm(resid(kcoarse__lmer)) #  OK
vif(kcoarse__lmer) # OK

#Run model based on predictions
kcoarse_dredge_lmer <- lmer(log(k_dday_coarse) ~ season*percent_open_downstream + season*percent_open_downstream +  mean_velocity_m_s + substrate_complexity + position_d*season + dist_to_source_km +  (1 | site ) , data = dat_lm_kcoarse_all_scaled) #try dropping percent_open_upstream 
summary(kcoarse_dredge_lmer)
plot(kcoarse_dredge_lmer) 
qqnorm(resid(kcoarse_dredge_lmer)) 
vif(kcoarse_dredge_lmer)


#Run dredge
options(na.action = "na.fail") # Required for dredge to run

kcoarse__lmer_dredge<- dredge(kcoarse_dredge_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default


nrow(kcoarse__lmer_dredge)  #how many models were run 152
head(kcoarse__lmer_dredge)

top.kcoarse__lmer_dredge  <- get.models(kcoarse__lmer_dredge , subset=delta <= 2) #chose 4 because null model otherwise 
mavg <-  model.avg(top.kcoarse__lmer_dredge)

step(kcoarse_dredge_lmer)

#Run the dredge model of the top averaged model
kcoarse_dredge_lmer <- lmer(log(k_dday_coarse) ~ season +  (1 | site ) , data = dat_lm_kcoarse_all_scaled_out)

summary(kcoarse_dredge_lmer)


#### Run LMs by season ####

kcoarse_spring_lm <- lm(log(k_dday_coarse) ~ dist_to_source_km, data=subset(dat_means, season=="Fall" & position_d=="downstream"))
summary(kcoarse_spring_lm)
plot(kcoarse_spring_lm)

kcoarse_lm <- lm(log(k_dday_coarse) ~ dist_to_source_km, data=subset(dat_means, season=="Summer" & k_dday_coarse <=0.013))
summary(kcoarse_lm)  #Note: with the highest value of kcoarse removed in the summer, the relationship is no longer significant.
plot(kcoarse_lm)

kcoarse_pos_lm <- lm(log(k_dday_coarse) ~ position_d, data=subset(dat_means, season=="Fall"))
summary(kcoarse_pos_lm)
plot(kcoarse_pos_lm)

emmeans(kcoarse_pos_lm, pairwise ~ position_d) #post hoc test

kcoarse_season_lm <- lm(log(k_dday_coarse) ~ season, data=dat_means)
summary(kcoarse_season_lm)
plot(kcoarse_season_lm)

emmeans(kcoarse_season_lm, pairwise ~ season) 

mean(dat_means$`k_dday_coarse`[dat_means$season == "Spring"])


############################################################################

#### Run LMM for k fine upstream  ####

dat_lm_kfine_up <- dat_means %>% 
  filter(site %in% upstream_sites)  %>%
  select( -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position)

#scale the data, except for CH4
dat_lm_kfine_up_scaled <- dat_lm_kfine_up %>%
  select(-k_dday_fine) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_kfine_up_scaled <- cbind(dat_lm_kfine_up_scaled, k_dday_fine = dat_lm_kfine_up$k_dday_fine)

#Test model relationships linear vs polynomial

#fit linear mixed effects model 
kfine_lmer<- lmer(log(k_dday_fine)  ~ dist_to_source_km  + (1 | site), data = dat_lm_kfine_up_scaled)
summary(kfine_lmer)
plot(kfine_lmer)
qqnorm(resid(kfine_lmer)) 
r.squaredGLMM(kfine_lmer) #.003
AIC(kfine_lmer) #20

kfine_lmer_poly<- lmer(log(k_dday_fine) ~ poly(dist_to_source_km,2)  + (1 | site), data = dat_lm_kfine_up_scaled)
summary(kfine_lmer_poly)
plot(kfine_lmer_poly)
qqnorm(resid(kfine_lmer_poly)) 
r.squaredGLMM(kfine_lmer_poly) #.04
AIC(kfine_lmer_poly) #18


############################################################################

#### Run LMM for k fine downstream  ####

dat_lm_kfine_d <- dat_means %>% 
  filter(site %in% downstream_sites)  %>%
  select( -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -position)

#scale the data, except for CH4
dat_lm_kfine_d_scaled <- dat_lm_kfine_d %>%
  select(-k_dday_fine) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_kfine_d_scaled <- cbind(dat_lm_kfine_d_scaled, k_dday_fine = dat_lm_kfine_d$k_dday_fine)

#Test model relationships linear vs polynomial

#fit linear mixed effects model 
kfine_lmer_d<- lmer(log(k_dday_fine)  ~ dist_to_source_km  + (1 | site), data = dat_lm_kfine_d_scaled)
summary(kfine_lmer_d)
plot(kfine_lmer_d)
qqnorm(resid(kfine_lmer_d)) 
r.squaredGLMM(kfine_lmer_d) #.0003
AIC(kfine_lmer_d) #7.9

kfine_lmer_poly<- lmer(log(k_dday_fine) ~ poly(dist_to_source_km,2)  + (1 | site), data = dat_lm_kfine_d_scaled)
summary(kfine_lmer_poly)
plot(kfine_lmer_poly)
qqnorm(resid(kfine_lmer_poly)) 
r.squaredGLMM(kfine_lmer_poly) #.04
AIC(kfine_lmer_poly) #8


##############################################################################

#### Use all data  with k fine LMM ####

dat_lm_kfine_all <- dat_means %>% 
  select( -date,  -campaign, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CO2_C_mg_m2_h, -CH4_C_mg_m2_h, -flow_state, -cobble, -boulder, -pebble_gravel, -water_temp_C)

#scale the data, except for CH4
dat_lm_kfine_all_scaled <- dat_lm_kfine_all %>%
  select(-k_dday_fine) %>%
  mutate(across(where(is.numeric), ~ scale(., center = T) %>% as.vector()))

# Add the OM stock column back
dat_lm_kfine_all_scaled <- cbind(dat_lm_kfine_all_scaled, k_dday_fine = dat_lm_kfine_all$k_dday_fine)

#dredge and step don't like the '.' in the formula, so put in all of the columns...
names(dat_lm_kfine_all_scaled)[sapply(dat_lm_kfine_all_scaled, is.numeric)]


kfine_all_lmer_sat <- lmer(k_dday_fine ~ season +  + water_pH + water_conductivity_us_cm + DO_. +  canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm + position*OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + temp_C + masl + position + bedrock + fine_substrate + (1 | site ) , data = dat_lm_kfine_all_scaled)

summary(kfine_all_lmer_sat)
plot(kfine_all_lmer_sat) # 
qqnorm(resid(kfine_all_lmer_sat)) # 1 outlier
vif(kfine_all_lmer_sat) # Remove DO_mg_L and dist_to_source_km

#Check outliers with Cook's Distance
cooksD <- cooks.distance(kfine_all_lmer_sat) #run the model with dat_sub first
influential <- cooksD[(cooksD > (5* mean(cooksD, na.rm = TRUE)))]
influential #6, 29

n <- nrow(dat_lm_kfine_all_scaled)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 15/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_lm_kfine_all_scaled[names_of_influential,]
dat_lm_kfine_all_scaled_out <- dat_lm_kfine_all_scaled %>% anti_join(outliers)

#Re run without outliers
kfine_all_lmer_sat2 <- lmer(k_dday_fine ~ season +  + water_pH + water_conductivity_us_cm + DO_. +  canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + mean_depth_cm + position*OM_stock_g_m2 + OM_flux_g_m2_s + sed_OM_percent + temp_C + masl + position + bedrock + fine_substrate +  (1 | site ) , data = dat_lm_kfine_all_scaled_out)

summary(kfine_all_lmer_sat2)
plot(kfine_all_lmer_sat2) # Much better
qqnorm(resid(kfine_all_lmer_sat2)) # still 1 outlier
vif(kfine_all_lmer_sat2) # 


lmerTest::step(kfine_all_lmer_sat2)

#Top model found:
#k_dday_fine ~ position + sed_OM_percent

#Run model
kfine_all_lmer <- lmer(k_dday_fine ~ position + sed_OM_percent + (1 | site ) , data = dat_lm_kfine_all_scaled_out)

summary(kfine_all_lmer)
plot(kfine_all_lmer) # 
qqnorm(resid(kfine_all_lmer)) # 
vif(kfine_all_lmer) # 


#Run model based on predictions

kfine__lmer <- lmer(k_dday_fine ~ position*poly(dist_to_source_km,2) + position*percent_open_upstream + position*percent_open_downstream + position_d*DailyMeanWaterTemp_C + (1 | site ) , data = dat_lm_kfine_all_scaled_out)

summary(kfine__lmer)
plot(kfine__lmer) # OK
qqnorm(resid(kfine__lmer)) #  OK, maybe 1 outlier
 
options(na.action = "na.fail") # Required for dredge to run

kfine__lmer_dredge<- dredge(kfine_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(kfine__lmer_dredge)  #35
head(kfine__lmer_dredge)

top.kfine__lmer_dredge <- get.models(kfine__lmer_dredge , subset=delta <= 2) #chose 4 because null model otherwise 
mavg <-  model.avg(top.kfine__lmer_dredge)

#Run model based on predictions

kfine_lmer <- lmer(log(k_dday_fine) ~ season*percent_open_downstream + season*percent_open_upstream +  DailyMeanWaterTemp_C + dist_to_source_km + position_d + (1 | site ) , data = dat_lm_kfine_all_scaled)
summary(kfine_lmer)
plot(kfine_lmer)
vif(kfine_lmer)

#### Run LMs by season ###

kfine_spring_lm <- lm(log(k_dday_fine) ~ dist_to_source_km, data=subset(dat_means, season=="Fall" ))
summary(kfine_spring_lm)
plot(kfine_spring_lm)

kfine_pos_lm <- lm(log(k_dday_fine) ~ position_d, data=subset(dat_means, season=="Fall"))
summary(kfine_pos_lm)
plot(kfine_pos_lm)

emmeans(kfine_pos_lm, pairwise ~ position_d) #post hoc test

kfine_season_lm <- lm(log(k_dday_fine) ~ season, data=dat_means)
summary(kfine_season_lm)
plot(kfine_season_lm)

emmeans(kfine_season_lm, pairwise ~ season) 

##################################################################################
##### For Fanny: LMMs with k per day (not dday) ####

# k coarse

#All data
kcoarse_source_lm <- lm(log(k_day_coarse) ~ dist_to_source_km, data=dat_means)
summary(kcoarse_source_lm) #not sig
plot(kcoarse_source_lm)

kcoarse_dam_lm <- lm(log(k_day_coarse) ~ dist_to_dam_km, data=dat_means)
summary(kcoarse_dam_lm) #not sig
plot(kcoarse_dam_lm)

kcoarse_season_lm <- lm(log(k_dday_coarse) ~ season, data=dat_means)
summary(kcoarse_season_lm)
plot(kcoarse_season_lm)

emmeans(kcoarse_season_lm, pairwise ~ season) #Spring - Summer and Spring - Fall sig diff


#by season
kcoarse_spring_lm <- lm(log(k_day_coarse) ~ dist_to_dam_km, data=subset(dat_means, season=="Summer" & k_day_coarse <= 0.22))
summary(kcoarse_spring_lm)
plot(kcoarse_spring_lm)

#by position
kcoarse_pos_lm <- lm(log(k_day_coarse) ~ position_d, data=subset(dat_means, season=="Fall"))
summary(kcoarse_pos_lm)
plot(kcoarse_pos_lm)

emmeans(kcoarse_pos_lm, pairwise ~ position_d)

#kfine

kfine_source_lm <- lm(log(k_day_fine) ~ dist_to_dam_km, data=dat_means)
summary(kfine_source_lm) #not sig
plot(kfine_source_lm)

kfine_season_lm <- lm(log(k_day_fine) ~ season, data=dat_means)
summary(kfine_season_lm) 
emmeans(kfine_season_lm, pairwise ~ season) #all seasons sig  different

#by season

kfine_season_lm <- lm(log(k_day_fine) ~ dist_to_source_km, data=subset(dat_means, season=="Fall" & position_d=="downstream"))
summary(kfine_season_lm)
plot(kfine_season_lm)

#by position
kfine_pos_lm <- lm(log(k_day_fine) ~ position_d, data=subset(dat_means, season=="Fall"))
summary(kfine_pos_lm)
plot(kfine_pos_lm)

emmeans(kfine_pos_lm, pairwise ~ position_d) #

##############################################################################

#### predict CO2 values ####

#subset downstream sites
dat_predict <- dat %>% 
  filter(site %in% downstream_sites)  %>%
           select(-date, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state,  -cobble, -boulder, -pebble_gravel, -DO_mg_L, -water_temp_C, -dist_ds_dam_km, -dist_to_source_km, -position)

names(dat_predict)[sapply(dat_predict, is.numeric)]

lme4:::predict.merMod(CO2_lmer_up, newdata = dat_predict, allow.new.levels = TRUE)

dat_predict$predicted_CO2_log <- predict(CO2_lmer_up, newdata = dat_predict, allow.new.levels = TRUE)
dat_predict$predicted_CO2_original <- exp(dat_predict$predicted_CO2 - 16.3819)


#############################################################################

#### Run a PCA ####


dat_pca <- dat %>% 
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine,-flow_state, -season, -cobble, -boulder, -pebble_gravel, -water_temp_C)

dat_pca2 <- dat %>% 
 select(water_temp_C, water_pH, water_conductivity_us_cm, DO_mg_L, DO_., canopy_cover_., wetted_width_m, discharge_m3_s, mean_velocity_m_s, DailyMeanWaterTemp_C, temp_C, substrate_complexity, mean_depth_cm, percent_open_upstream, percent_open_downstream, masl, OM_stock_g_m2, OM_flux_g_m2_s, dist_to_source_km) #CO2_C_mg_m2_h, CH4_C_mg_m2_h, k_dday_coarse, k_dday_fine)

dat_pca3 <- dat %>% 
select(water_pH ,water_conductivity_us_cm , DO_mg_L, DO_. , canopy_cover_. , wetted_width_m , discharge_m3_s , mean_velocity_m_s , DailyMeanWaterTemp_C , temp_C , substrate_complexity , fine_substrate , mean_depth_cm ,  masl , dist_to_source_km , OM_stock_g_m2)  %>% 
rename(
  pH = water_pH,
  Conductivity = water_conductivity_us_cm,
  DO = DO_mg_L,
  DO_percent = DO_.,  
  Canopy_Cover = canopy_cover_.,
  Wetted_Width = wetted_width_m,
  Discharge = discharge_m3_s,
  Velocity = mean_velocity_m_s,
  Water_Temp = DailyMeanWaterTemp_C,
  Air_Temp = temp_C,
  Substrate_Complexity = substrate_complexity,
  Fine_Substrate = fine_substrate,
  Mean_Depth = mean_depth_cm,
  masl = masl,
  Dist_to_Source = dist_to_source_km,
  OM_Stock = OM_stock_g_m2
)

dat_PCA <- dat_pca3[, sapply(dat_pca3, is.numeric)]

PCA_CO2 <- PCA(dat_PCA, scale.unit = TRUE, ncp = 5, graph = TRUE)
PCA_CO2

my.col.var <- c("#91278e", "#3d9970", "#f7921e")



tiff("PCA_by_position", units="in", width=8, height=6, res=300)
PCA_fig <- fviz_pca_biplot(PCA_CO2, 
                                   col.ind = dat$position_d,
                                   addEllipses = TRUE, label = "var",
                                   pointsize=3,
                                   alpha.ind=0.25,
                                   mean.point=F,
                                   palette=my.col.var,
                                   col.var = "black", repel = TRUE,
                                   legend.title = " ") + ggtitle(NULL) 
PCA_fig
dev.off()

#There are a lot of varialbes so maybe we can drop the less influential variables based on their eigenvalues

eig.val <- get_eigenvalue(PCA_CO2)
eig.val #to view egien values. An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data. This is commonly used as a cutoff point for which PCs are retained

fviz_eig(PCA_CO2, addlabels = TRUE, ylim = c(0, 25)) #visualize the percentage of explained variance per dimension

#test groupwise differences using MANOVA

# Extract the PCA scores from PCA_CO2
PCA_scores <- PCA_CO2$ind$coord

# Create a data frame with PCA scores and position_d
dat_with_PCA <- data.frame(PCA1 = PCA_scores[, 1], PCA2 = PCA_scores[, 2], position_d = dat$position_d)

# Perform MANOVA
fit <- manova(cbind(PCA1, PCA2) ~ position_d, data = dat_with_PCA)
summary(fit) #significant effect of position_d

#permanova is non-parametric alternative

result <- adonis2(dat_with_PCA[, c("PCA1", "PCA2")] ~ dat_with_PCA$position_d, permutations = 999, method = "euclidean")

#pairwise differences

pairwise_result <- pairwise.perm.manova(resp = dat_with_PCA[, c("PCA1", "PCA2")],
                                        fact = dat_with_PCA$position_d, nperm = 999)

#variable percent contributions

var <- get_pca_var(PCA_CO2)

head(var$contrib, 100)

# Contributions of variables to PC1
fviz_contrib(PCA_CO2, choice = "var", axes = 1, top = 10)
fviz_contrib(PCA_CO2, choice = "var", axes = 2, top = 10)

#The variables most strongly associated with standing water are: fine substrate, mean depth, water_conductivity, wetted width, sed_OM percent, DO mg L, OM stock,  and OM flux.

#The variables most strongly positively associated with upstream are: masl, canopy cover, k_dday_coarse, and k_dday_fine

#run a PCA with the substrate variables

dat_substrate <- dat %>%
select(bedrock, boulder, cobble, pebble_gravel, fine_substrate)

dat_substrate <- lapply(dat_substrate, function(x) asin(sqrt(x/100)))


PCA_substrate <- PCA(dat_substrate, scale.unit = TRUE, ncp = 5, graph = TRUE)
PCA_substrate

PCA_substrate_vis <- fviz_pca_biplot(PCA_substrate, 
                               col.ind = dat$position_d, #change for flow_state
                               addEllipses = TRUE, label = "var",
                               pointsize=3,
                               alpha.ind=0.4,
                               mean.point=F,
                               ellipse.method = "covariance",
                               col.var = "black", repel = TRUE,
                               legend.title = " ") + ggtitle(NULL) 
PCA_substrate_vis

#####################################################################
#### Try structural equation modeling SEM ####

library(lavaan)
library(semPlot)

dat <-as.data.frame(dat)

row.names(dat) <- NULL  # Remove existing row names
row.names(dat) <- 1:nrow(dat)  # Assign new row names

#we want to create two latent variables which are "dam effect" and "network effect", to see how these impact CO2 fluxes

#Fit a test model
mCO2 <- '
# measurement model
dam =~ fine_substrate + mean_depth_cm + water_conductivity_us_cm + wetted_width_m + OM_stock_g_m2 
network =~ discharge_m3_s + masl

# regressions
CO2_C_mg_m2_h ~ dam + network

# variances and covariances 
               discharge_m3_s ~~ discharge_m3_s
    
               dam ~~ network 
'
fitmCO2 <- sem(mCO2, data=dat)

summary(fitmCO2, standardized=TRUE, fit.measures=TRUE)


#get p values of coefficients
table2<-parameterEstimates(fitmCO2,standardized=TRUE)  %>%  head(9)

#turning the chosen parameters into text
b<-gettextf('%.3f \n p=%.3f', table2$est, digits=table2$pvalue)

node_labels <- c("Fine\nsubstrate", "Depth", "Conduct\nivity", "Width", "OM\nstock", "masl", "Discharge", "CO2", "Dam", "Network")

# Create path diagram with adjusted settings
semPaths(fitmCO2, intercept = FALSE, 
         whatLabel = "est",
         what="std",
         edgeLabels = b,
        # edge.color= ifelse(table2$est > 0, "green", "red"), 
         residuals = FALSE, exoCov = FALSE, 
         nodeLabels = node_labels,
         style='ram', edge.label.cex = 1.3,
         layout = "tree", sizeMan = 7)


semPaths(fitmCO2, "std", edge.label.cex = 1.0, curvePivot = TRUE) #nice layout

################################################################################

#### Run Mantel tests ####

#response variable matrices
CO2flux <- dat %>%
  select(CO2_C_mg_m2_h)

CH4flux <- dat %>%
  select(CH4_C_mg_m2_h)

stock_OM <- dat_means %>%
  select(OM_stock_g_m2)

k_coarse <- dat_means %>%
  select(k_dday_coarse)

k_fine <- dat_means %>%
  select(k_dday_fine)


#Make a matrix of euclidean distances
geo <- data.frame(dat$lon, dat$lat)
geo_means <- data.frame(dat_means$lon, dat_means$lat)

#Load in network distances (calculated by Herve Pella)
net_dist <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Network distances/network_dist_tarentaine.csv")

#Need to make it into a matrix
site <- sort(unique(c(net_dist$nom1, net_dist$nom2)))
matrix_data <- matrix(NA, nrow = length(site), ncol = length(site), dimnames = list(site, site))

# Fill the matrix with distances
for (i in 1:nrow(net_dist)) {
  site1 <- net_dist$nom1[i]
  site2 <- net_dist$nom2[i]
  distance <- net_dist$distance.m[i]
  
  matrix_data[site1, site2] <- distance
  matrix_data[site2, site1] <- distance
}

diag(matrix_data) <- 0

#rename
net_dist <- matrix_data

#write as csv

write.csv(matrix_data, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Network distances/network_dist_tarentaine_matrix.csv")

#make a vector for damming metrics
dam.metrics <-  dat %>%
  select(percent_open_downstream, percent_open_downstream, dist_ds_dam_km)

#make a vector for damming metrics
dam.metrics.means <-  dat_means %>%
  select(percent_open_downstream, percent_open_downstream, dist_ds_dam_km)

#environmental vector
env <- dat  %>%
  select(water_temp_C, water_pH, water_conductivity_us_cm, DO_mg_L, DO_., canopy_cover_., wetted_width_m, discharge_m3_s, mean_velocity_m_s, DailyMeanWaterTemp_C, temp_C, substrate_complexity)

env_means <- dat_means  %>%
  select(water_temp_C, water_pH, water_conductivity_us_cm, DO_mg_L, DO_., canopy_cover_., wetted_width_m, discharge_m3_s, mean_velocity_m_s, DailyMeanWaterTemp_C, temp_C, substrate_complexity)

#scale environmental data
scale.env <- scale(env, center = TRUE, scale = TRUE)
scale.env.means <- scale(env_means, center = TRUE, scale = TRUE)

#can also make a vector for dam matrices

#Now we have to convert these subsets into distance matrices.

#for the response variable
dist.CO2flux <- dist(CO2flux, method = "euclidean")
dist.CH4flux <- dist(CH4flux, method = "euclidean")
dist.OMstock <- dist(stock_OM, method = "euclidean")
dist.kcoarse <- dist(k_coarse, method = "euclidean")
dist.kfine <- dist(k_fine, method = "euclidean")
dist.dam.metrics <- dist(dam.metrics, method = "euclidean")
dist.dam.metrics.means <- dist(dam.metrics.means, method = "euclidean")

#environmental vector - euclidean distance - or bray? Tutorial says euclidean
dist.env <- dist(scale.env, method = "euclidean")
dist.env.means <- dist(scale.env.means, method = "euclidean")

#geographic distance - Haversine distance

#geographic data frame - haversine distance 
d.geo = distm(geo, fun = distHaversine)
dist.geo.hav = as.dist(d.geo)

d.geo.means = distm(geo_means, fun = distHaversine)
dist.geo.hav.means = as.dist(d.geo.means)

#now we can run the Mantel command

#do environmental variables vary with geographic distance?
envgeo <- mantel(dist.env, dist.geo.hav, method = "spearman", permutations = 999, na.rm = TRUE)
envgeo #Mantel statistic r: 0.1762 Significance: 0.001 


#does the response variable vary with environmental distance?
CO2_env <- mantel(dist.CO2flux, dist.env, method = "spearman", permutations = 999, na.rm = TRUE)
CO2_env #if taking a super long time, try reducing permutations 
#Mantel statistic r: 0.1031 Significance: 0.001 

CH4_env <- mantel(dist.CH4flux, dist.env, method = "spearman", permutations = 999, na.rm = TRUE)
CH4_env #Mantel statistic r: 0.1166   Significance: 0.003

OMstock_env <- mantel(dist.OMstock, dist.env.means, method = "spearman", permutations = 999, na.rm = TRUE)
OMstock_env  #Mantel statistic r: 0.3276 Significance: 0.001 

kcoarse_env <- mantel(dist.kcoarse, dist.env.means, method = "spearman", permutations = 999, na.rm = TRUE)
kcoarse_env #Mantel statistic r: -0.008908  Significance: 0.508

kfine_env <- mantel(dist.kfine, dist.env.means, method = "spearman", permutations = 999, na.rm = TRUE)
kfine_env  #Mantel statistic r: -0.06784  Significance: 0.819 

# Run a partial Mantel test
#Which is fluxes vs geographical distance, accounting for env distance

CO2partial <- mantel.partial(dist.CO2flux, dist.env, dist.geo.hav, method="pearson", permutations=999)
CO2partial #this is removing any spatial autocorrelation if the env matrix explains the flux matrix
#Mantel statistic r: 0.2208     Significance: 0.004

CH4partial <- mantel.partial(dist.CH4flux, dist.env, dist.geo.hav, method="pearson", permutations=999)
CH4partial # Mantel statistic r: 0.08764 Significance: 0.043

OMstockpartial <- mantel.partial(dist.OMstock, dist.env.means, dist.geo.hav.means, method="pearson", permutations=999)
OMstockpartial #Mantel statistic r: 0.3048 Significance: 0.022 

kcoarsepartial <- mantel.partial(dist.kcoarse, dist.env.means, dist.geo.hav.means, method="pearson", permutations=999)
kcoarsepartial  #Mantel statistic r: -0.02291   Significance: 0.561

kfinepartial <- mantel.partial(dist.kfine, dist.env.means, dist.geo.hav.means, method="pearson", permutations=999)
kfinepartial  #  Mantel statistic r: -0.1001 Significance: 0.891

CO2partial2 <- mantel.partial(dist.CO2flux, dist.geo.hav, dist.env, method="pearson", permutations=999)
CO2partial2 ##This is removing the effect of env variables if euclidean distance explains fluxes
#Mantel statistic r: -0.03577 Significance: 0.871

CH4partial2 <- mantel.partial(dist.CH4flux, dist.geo.hav, dist.env, method="pearson", permutations=999)
CH4partial2 #Mantel statistic r: -0.126  Significance: 1

OMstockpartial2 <- mantel.partial(dist.OMstock, dist.geo.hav.means, dist.env.means, method="pearson", permutations=999)
OMstockpartial2 #Mantel statistic r: -0.1148 Significance: 0.996 

kcoarsepartial2 <- mantel.partial(dist.kcoarse, dist.geo.hav.means, dist.env.means, method="pearson", permutations=999)
kcoarsepartial2  #Mantel statistic r: 0.05975  Significance: 0.227

kfinepartial2 <- mantel.partial(dist.kfine, dist.geo.hav.means, dist.env.means, method="pearson", permutations=999)
kfinepartial2  #  Mantel statistic r: 0.008369  Significance: 0.409 

CO2partial3 <- mantel.partial(dist.CO2flux, dist.dam.metrics, dist.geo.hav, method="pearson", permutations=999) #Mantel statistic r: -0.06446   Significance: 0.998
CO2partial3 ##This is removing the effect of env distance if damming variables explain fluxes

CH4partial3 <- mantel.partial(dist.CH4flux, dist.dam.metrics, dist.geo.hav, method="pearson", permutations=999)
CH4partial3  #Mantel statistic r: 0.1569  Significance: 0.001

OMstockpartial3 <- mantel.partial(dist.OMstock, dist.dam.metrics.means, dist.geo.hav.means, method="pearson", permutations=999)
OMstockpartial3 #Mantel statistic r: 0.008879 Significance: 0.34

kcoarsepartial3 <- mantel.partial(dist.kcoarse, dist.dam.metrics.means, dist.geo.hav.means, method="pearson", permutations=999)
kcoarsepartial3 #Mantel statistic r: -0.1215  Significance: 0.99 

kfinepartial3 <- mantel.partial(dist.kfine, dist.dam.metrics.means, dist.geo.hav.means, method="pearson", permutations=999)
kfinepartial3 #Mantel statistic r: -0.09975  Significance: 0.937


######################################################################
#### Mantel test: effects of network distance? must do this by campaign ####

#Make response variable data frame per campaign, with site averages to run against the network distances

#response variable matrices
CO2flux_1 <- dat_means %>%
  filter(season == "Spring") %>%
  select(CO2_C_mg_m2_h) # 20 obs

CO2flux_2 <- dat_means %>%
  filter(season == "Summer") %>%
  select(CO2_C_mg_m2_h)  #19 obs #No TA12 for C2

CO2flux_3 <- dat_means %>%
  filter(season == "Fall") %>%
  select(CO2_C_mg_m2_h)  #19 obs #no TA02 for C3, instead we have TA02, and no TA12 for C3

CH4flux_1 <- dat_means %>%
  filter(season == "Spring") %>%
  select(CH4_C_mg_m2_h)  #20 obs

CH4flux_2 <- dat_means %>%
  filter(season == "Summer") %>%
  select(CH4_C_mg_m2_h)  #19 obs

CH4flux_3 <- dat_means %>%
  filter(season == "Fall") %>%
  select(CH4_C_mg_m2_h)   #19 obs

stock_OM_1 <- dat_means %>%
  filter(season == "Spring") %>%
  select(OM_stock_g_m2) #20 obs

stock_OM_2 <- dat_means %>%
  filter(season == "Summer") %>%
  select(OM_stock_g_m2) #19 obs

stock_OM_3 <- dat_means %>%
  filter(season == "Fall") %>%
  select(OM_stock_g_m2) #19 obs

k_coarse_1 <- dat_means %>%
  filter(season == "Spring") %>%
  select(k_dday_coarse)   #20 obs

k_coarse_2 <- dat_means %>%
  filter(season == "Summer") %>%
  select(k_dday_coarse)   #19 obs

k_coarse_3 <- dat_means %>%
  filter(season == "Fall") %>%
  select(k_dday_coarse)   #19 obs

k_fine_1 <- dat_means %>%
  filter(season == "Spring") %>%
  select(k_dday_fine)   #20 obs

k_fine_2 <- dat_means %>%
  filter(season == "Summer") %>%
  select(k_dday_fine)   # 19 obs
 
k_fine_3 <- dat_means %>%
  filter(season == "Fall") %>%
  select(k_dday_fine)   #19 obs

#Make distance matrix

#for the response variable
dist.CO2flux_1 <- dist(CO2flux_1, method = "euclidean")
dist.CO2flux_2 <- dist(CO2flux_2, method = "euclidean")
dist.CO2flux_3 <- dist(CO2flux_3, method = "euclidean")
dist.CH4flux_1 <- dist(CH4flux_1, method = "euclidean")
dist.CH4flux_2 <- dist(CH4flux_2, method = "euclidean")
dist.CH4flux_3 <- dist(CH4flux_3, method = "euclidean")
dist.stock_OM_1 <- dist(stock_OM_1, method = "euclidean")
dist.stock_OM_2 <- dist(stock_OM_2, method = "euclidean")
dist.stock_OM_3 <- dist(stock_OM_3, method = "euclidean")
dist.k_coarse_1 <- dist(k_coarse_1, method = "euclidean")
dist.k_coarse_2 <- dist(k_coarse_2, method = "euclidean")
dist.k_coarse_3 <- dist(k_coarse_3, method = "euclidean")
dist.k_fine_1 <- dist(k_fine_1, method = "euclidean")
dist.k_fine_2 <- dist(k_fine_2, method = "euclidean")
dist.k_fine_3 <- dist(k_fine_3, method = "euclidean")

#env var matrix by campaign
env_means_1 <- dat_means  %>%
  filter(season == "Spring") %>%
  select(water_temp_C, water_pH, water_conductivity_us_cm, DO_mg_L, DO_., canopy_cover_., wetted_width_m, discharge_m3_s, mean_velocity_m_s, DailyMeanWaterTemp_C, temp_C, substrate_complexity)

env_means_2 <- dat_means  %>%
  filter(season == "Summer") %>%
  select(water_temp_C, water_pH, water_conductivity_us_cm, DO_mg_L, DO_., canopy_cover_., wetted_width_m, discharge_m3_s, mean_velocity_m_s, DailyMeanWaterTemp_C, temp_C, substrate_complexity)

env_means_3 <- dat_means  %>%
  filter(season == "Fall") %>%
  select(water_temp_C, water_pH, water_conductivity_us_cm, DO_mg_L, DO_., canopy_cover_., wetted_width_m, discharge_m3_s, mean_velocity_m_s, DailyMeanWaterTemp_C, temp_C, substrate_complexity)

#scale environmental data
scale.env.means1 <- scale(env_means_1, center = TRUE, scale = TRUE)
scale.env.means2 <- scale(env_means_2, center = TRUE, scale = TRUE)
scale.env.means3 <- scale(env_means_3, center = TRUE, scale = TRUE)

dist.env.means1 <- dist(scale.env.means1, method = "euclidean")
dist.env.means2 <- dist(scale.env.means2, method = "euclidean")
dist.env.means3 <- dist(scale.env.means3, method = "euclidean")

#damming metrics by campaign
#make a vector for damming metrics
dam.metrics1 <-  dat_means %>%
  filter(season == "Spring") %>%
  select(percent_open_downstream, percent_open_downstream, dist_ds_dam_km)

dam.metrics2 <-  dat_means %>%
  filter(season == "Summer") %>%
  select(percent_open_downstream, percent_open_downstream, dist_ds_dam_km)

dam.metrics3 <-  dat_means %>%
  filter(season == "Fall") %>%
  select(percent_open_downstream, percent_open_downstream, dist_ds_dam_km)

#scale
scale.dam.metrics1 <- scale(dam.metrics1, center = TRUE, scale = TRUE)
scale.dam.metrics2 <- scale(dam.metrics2, center = TRUE, scale = TRUE)
scale.dam.metrics3 <- scale(dam.metrics3, center = TRUE, scale = TRUE)

dist.dam.metrics1 <- dist(scale.dam.metrics1, method = "euclidean")
dist.dam.metrics2 <- dist(scale.dam.metrics2, method = "euclidean")
dist.dam.metrics3 <- dist(scale.dam.metrics3, method = "euclidean")

#Make a custom network distance matrix for each campaign depending on number of observations

#Spring, remove TA02R
net_dist1 <- subset(net_dist, select = -c(TA02R) )  #remove the site TA02R column
net_dist1 <-  net_dist1[-3,]   #and row

#Summer, remove TA12 and TA02R
net_dist2 <- subset(net_dist, select = -c(TA02R) )  #remove the site TA02R column
net_dist2 <-  net_dist2[-3,]   #and row
net_dist2 <- subset(net_dist2, select = -c(TA12) )   #remove site TA12 column
net_dist2 <-  net_dist2[-12,]  #... and row

#Fall, remove TA12 and TA02
net_dist3 <- subset(net_dist, select = -c(TA02) )  #remove the site TA02 column
net_dist3 <-  net_dist3[-2,]  
net_dist3 <- subset(net_dist3, select = -c(TA12) )  
net_dist3 <-  net_dist3[-12,] 

#Geographic distance
geo_means1 <- dat_means %>%
  filter(site != "TA02R" & season=="Spring") %>%
  select(lon, lat)

geo_means2 <- dat_means %>%
  filter(site != "TA12" & site != "TA02R" & season=="Summer") %>%
  select(lon, lat)

geo_means3 <- dat_means %>%
  filter(site != "TA12" & site != "TA02" & season=="Fall") %>%
  select(lon, lat)

#make matrix of geographic distances
d.geo.means1 = distm(geo_means1, fun = distHaversine)
dist.geo1 = as.dist(d.geo.means1)

d.geo.means2 = distm(geo_means2, fun = distHaversine)
dist.geo2 = as.dist(d.geo.means2)

d.geo.means3 = distm(geo_means3, fun = distHaversine)
dist.geo3 = as.dist(d.geo.means3)


## Run mantel tests

#effect of network distances

env_C1 <-  mantel(dist.env.means1, net_dist1, method = "spearman", permutations = 999)
env_C1 #not significant

env_C2 <-  mantel(dist.env.means2, net_dist2, method = "spearman", permutations = 999)
env_C2 #significant at 0.03

env_C3 <-  mantel(dist.env.means3, net_dist3, method = "spearman", permutations = 999)
env_C3 #not significant

# So we only need a partial Mantel test for C2...

geo_C1 <-  mantel(dist.geo1, net_dist1, method = "spearman", permutations = 999)
geo_C1 # p 0.001

geo_C2 <-  mantel(dist.geo2, net_dist2, method = "spearman", permutations = 999)
geo_C2 # p 0.001

geo_C3 <-  mantel(dist.geo3, net_dist3, method = "spearman", permutations = 999)
geo_C3 # p 0.001

#effect of dams on env vars

damsC1 <- mantel(dist.env.means1, dist.dam.metrics1, method = "spearman", permutations = 999)
damsC1

damsC2 <- mantel(dist.env.means2, dist.dam.metrics2, method = "spearman", permutations = 999)
damsC2

damsC3 <- mantel(dist.env.means3, dist.dam.metrics3, method = "spearman", permutations = 999)
damsC3

#effect of env vars on response vars

CO2_C1 <- mantel(dist.CO2flux_1, dist.env.means1, method = "spearman", permutations = 999)
CO2_C1 #not sig

CO2_C2 <- mantel(dist.CO2flux_2, dist.env.means2, method = "spearman", permutations = 999)
CO2_C2 #not sig

CO2_C3 <- mantel(dist.CO2flux_3, dist.env.means3, method = "spearman", permutations = 999)
CO2_C3  #sig

#effect of network distance on CO2

CO2net_C1 <- mantel(dist.CO2flux_1, net_dist1, method = "spearman", permutations = 999)
CO2net_C1 # not sig

CO2net_C2 <- mantel(dist.CO2flux_2, net_dist2, method = "spearman", permutations = 999)
CO2net_C2 # not sig

CO2net_C3 <- mantel(dist.CO2flux_3, net_dist3, method = "spearman", permutations = 999)
CO2net_C3 # not sig

#Effect of damming on CO2

CO2dam_C1 <- mantel(dist.CO2flux_1, dist.dam.metrics1, method = "spearman", permutations = 999)
CO2dam_C1 # not sig

CO2dam_C2 <- mantel(dist.CO2flux_2, dist.dam.metrics2, method = "spearman", permutations = 999)
CO2dam_C2 # not sig

CO2dam_C3 <- mantel(dist.CO2flux_3, dist.dam.metrics3, method = "spearman", permutations = 999)
CO2dam_C3 # not sig

## partial Mantel tests

CO2_partial1 <- mantel.partial(dist.CO2flux_1, dist.env.means1, net_dist1, method="pearson", permutations=999)
CO2_partial1 

CO2_partial2 <- mantel.partial(dist.CO2flux_2, dist.env.means2, net_dist2, method="pearson", permutations=999)
CO2_partial2 

CO2_partial3 <- mantel.partial(dist.CO2flux_3, dist.env.means3, net_dist3, method="pearson", permutations=999)
CO2_partial3 
#
#
CO2_part1 <- mantel.partial(dist.CO2flux_1,  net_dist1, dist.env.means1, method="pearson", permutations=999)
CO2_part1 

CO2_part2 <- mantel.partial(dist.CO2flux_2,  net_dist2, dist.env.means2,  method="pearson", permutations=999)
CO2_part2 

CO2_part3 <- mantel.partial(dist.CO2flux_3,  net_dist3, dist.env.means3, method="pearson", permutations=999)
CO2_part3 
#
#
CO2_par1 <- mantel.partial(dist.CO2flux_1, dist.dam.metrics1, net_dist1, method="pearson", permutations=999)
CO2_par1 

CO2_par2 <- mantel.partial(dist.CO2flux_2, dist.dam.metrics2,  net_dist2,  method="pearson", permutations=999)
CO2_par2 

CO2_par3 <- mantel.partial(dist.CO2flux_3, dist.dam.metrics3,  net_dist3,  method="pearson", permutations=999)
CO2_par3 

############################################################################
## CH4

#effect of env vars on response vars

CH4_C1 <- mantel(dist.CH4flux_1, dist.env.means1, method = "spearman", permutations = 999)
CH4_C1 #

CH4_C2 <- mantel(dist.CH4flux_2, dist.env.means2, method = "spearman", permutations = 999)
CH4_C2 #

CH4_C3 <- mantel(dist.CH4flux_3, dist.env.means3, method = "spearman", permutations = 999)
CH4_C3  #

#effect of network distance on CH4

CH4net_C1 <- mantel(dist.CH4flux_1, net_dist1, method = "spearman", permutations = 999)
CH4net_C1 # 

CH4net_C2 <- mantel(dist.CH4flux_2, net_dist2, method = "spearman", permutations = 999)
CH4net_C2 # 

CH4net_C3 <- mantel(dist.CH4flux_3, net_dist3, method = "spearman", permutations = 999)
CH4net_C3 # 

#Effect of damming on CH4

CH4dam_C1 <- mantel(dist.CH4flux_1, dist.dam.metrics1, method = "spearman", permutations = 999)
CH4dam_C1 # 

CH4dam_C2 <- mantel(dist.CH4flux_2, dist.dam.metrics2, method = "spearman", permutations = 999)
CH4dam_C2 # 

CH4dam_C3 <- mantel(dist.CH4flux_3, dist.dam.metrics3, method = "spearman", permutations = 999)
CH4dam_C3 # 

## partial Mantel tests

CH4_partial1 <- mantel.partial(dist.CH4flux_1, dist.env.means1, net_dist1, method="pearson", permutations=999)
CH4_partial1 

CH4_partial2 <- mantel.partial(dist.CH4flux_2, dist.env.means2, net_dist2, method="pearson", permutations=999)
CH4_partial2 

CH4_partial3 <- mantel.partial(dist.CH4flux_3, dist.env.means3, net_dist3, method="pearson", permutations=999)
CH4_partial3 
#
#
CH4_part1 <- mantel.partial(dist.CH4flux_1,  net_dist1, dist.env.means1, method="pearson", permutations=999)
CH4_part1 

CH4_part2 <- mantel.partial(dist.CH4flux_2,  net_dist2, dist.env.means2,  method="pearson", permutations=999)
CH4_part2 

CH4_part3 <- mantel.partial(dist.CH4flux_3,  net_dist3, dist.env.means3, method="pearson", permutations=999)
CH4_part3 
#
#
CH4_par1 <- mantel.partial(dist.CH4flux_1, dist.dam.metrics1, net_dist1, method="pearson", permutations=999)
CH4_par1 

CH4_par2 <- mantel.partial(dist.CH4flux_2, dist.dam.metrics2,  net_dist2,  method="pearson", permutations=999)
CH4_par2 

CH4_par3 <- mantel.partial(dist.CH4flux_3, dist.dam.metrics3,  net_dist3,  method="pearson", permutations=999)
CH4_par3 
#
#

####################################################
#OM

#effect of env vars on response vars

OM_C1 <- mantel(dist.stock_OM_1, dist.env.means1, method = "spearman", permutations = 999)
OM_C1 #

OM_C2 <- mantel(dist.stock_OM_2, dist.env.means2, method = "spearman", permutations = 999)
OM_C2 #

OM_C3 <- mantel(dist.stock_OM_3, dist.env.means3, method = "spearman", permutations = 999)
OM_C3  #

#effect of network distance on OM

OMnet_C1 <- mantel(dist.stock_OM_1, net_dist1, method = "spearman", permutations = 999)
OMnet_C1 # 

OMnet_C2 <- mantel(dist.stock_OM_2, net_dist2, method = "spearman", permutations = 999)
OMnet_C2 # 

OMnet_C3 <- mantel(dist.stock_OM_3, net_dist3, method = "spearman", permutations = 999)
OMnet_C3 # 

#Effect of damming on OM

OMdam_C1 <- mantel(dist.stock_OM_1, dist.dam.metrics1, method = "spearman", permutations = 999)
OMdam_C1 # 

OMdam_C2 <- mantel(dist.stock_OM_2, dist.dam.metrics2, method = "spearman", permutations = 999)
OMdam_C2 # 

OMdam_C3 <- mantel(dist.stock_OM_3, dist.dam.metrics3, method = "spearman", permutations = 999)
OMdam_C3 # 

## partial Mantel tests

OM_partial1 <- mantel.partial(dist.stock_OM_1, dist.env.means1, net_dist1, method="pearson", permutations=999)
OM_partial1 

OM_partial2 <- mantel.partial(dist.stock_OM_2, dist.env.means2, net_dist2, method="pearson", permutations=999)
OM_partial2 

OM_partial3 <- mantel.partial(dist.stock_OM_3, dist.env.means3, net_dist3, method="pearson", permutations=999)
OM_partial3 
#
#
OM_part1 <- mantel.partial(dist.stock_OM_1,  net_dist1, dist.env.means1, method="pearson", permutations=999)
OM_part1 

OM_part2 <- mantel.partial(dist.stock_OM_2,  net_dist2, dist.env.means2,  method="pearson", permutations=999)
OM_part2 

OM_part3 <- mantel.partial(dist.stock_OM_3,  net_dist3, dist.env.means3, method="pearson", permutations=999)
OM_part3 
#
#
OM_par1 <- mantel.partial(dist.stock_OM_1, dist.dam.metrics1, net_dist1, method="pearson", permutations=999)
OM_par1 

OM_par2 <- mantel.partial(dist.stock_OM_2, dist.dam.metrics2,  net_dist2,  method="pearson", permutations=999)
OM_par2 

OM_par3 <- mantel.partial(dist.stock_OM_3, dist.dam.metrics3,  net_dist3,  method="pearson", permutations=999)
OM_par3 

###############################################################

#k coarse

#effect of env vars on response vars

kcoarse_C1 <- mantel(dist.k_coarse_1, dist.env.means1, method = "spearman", permutations = 999)
kcoarse_C1 #

kcoarse_C2 <- mantel(dist.k_coarse_2, dist.env.means2, method = "spearman", permutations = 999)
kcoarse_C2 #

kcoarse_C3 <- mantel(dist.k_coarse_3, dist.env.means3, method = "spearman", permutations = 999)
kcoarse_C3  #

#effect of network distance on OM

kcoarsenet_C1 <- mantel(dist.k_coarse_1, net_dist1, method = "spearman", permutations = 999)
kcoarsenet_C1 # 

kcoarsenet_C2 <- mantel(dist.k_coarse_2, net_dist2, method = "spearman", permutations = 999)
kcoarsenet_C2 # 

kcoarsenet_C3 <- mantel(dist.k_coarse_3, net_dist3, method = "spearman", permutations = 999)
kcoarsenet_C3 # 

#Effect of damming on OM

kcoarsedam_C1 <- mantel(dist.k_coarse_1, dist.dam.metrics1, method = "spearman", permutations = 999)
kcoarsedam_C1 # 

kcoarsedam_C2 <- mantel(dist.k_coarse_2, dist.dam.metrics2, method = "spearman", permutations = 999)
OMdam_C2 # 

kcoarsedam_C3 <- mantel(dist.k_coarse_3, dist.dam.metrics3, method = "spearman", permutations = 999)
kcoarsedam_C3 # 

## partial Mantel tests

kcoarse_partial1 <- mantel.partial(dist.k_coarse_1, dist.env.means1, net_dist1, method="pearson", permutations=999)
kcoarse_partial1 

kcoarse_partial2 <- mantel.partial(dist.k_coarse_2, dist.env.means2, net_dist2, method="pearson", permutations=999)
kcoarse_partial2 

kcoarse_partial3 <- mantel.partial(dist.k_coarse_3, dist.env.means3, net_dist3, method="pearson", permutations=999)
kcoarse_partial3 
#
#
kcoarse_part1 <- mantel.partial(dist.k_coarse_1,  net_dist1, dist.env.means1, method="pearson", permutations=999)
kcoarse_part1 

kcoarse_part2 <- mantel.partial(dist.k_coarse_2,  net_dist2, dist.env.means2,  method="pearson", permutations=999)
kcoarse_part2 

kcoarse_part3 <- mantel.partial(dist.k_coarse_3,  net_dist3, dist.env.means3, method="pearson", permutations=999)
kcoarse_part3 
#
#
kcoarse_par1 <- mantel.partial(dist.k_coarse_1, dist.dam.metrics1, net_dist1, method="pearson", permutations=999)
kcoarse_par1 

kcoarse_par2 <- mantel.partial(dist.k_coarse_2, dist.dam.metrics2,  net_dist2,  method="pearson", permutations=999)
kcoarse_par2 

kcoarse_par3 <- mantel.partial(dist.k_coarse_3, dist.dam.metrics3,  net_dist3,  method="pearson", permutations=999)
kcoarse_par3 

############################################################
##kfine

#effect of env vars on response vars

kfine_C1 <- mantel(dist.k_fine_1, dist.env.means1, method = "spearman", permutations = 999)
kfine_C1 #

kfine_C2 <- mantel(dist.k_fine_2, dist.env.means2, method = "spearman", permutations = 999)
kfine_C2 #

kfine_C3 <- mantel(dist.k_fine_3, dist.env.means3, method = "spearman", permutations = 999)
kfine_C3  #

#effect of network distance on OM

kfinenet_C1 <- mantel(dist.k_fine_1, net_dist1, method = "spearman", permutations = 999)
kfinenet_C1 # 

kfinenet_C2 <- mantel(dist.k_fine_2, net_dist2, method = "spearman", permutations = 999)
kfinenet_C2 # 

kfinenet_C3 <- mantel(dist.k_fine_3, net_dist3, method = "spearman", permutations = 999)
kfinenet_C3 # 

#Effect of damming on OM

kfinedam_C1 <- mantel(dist.k_fine_1, dist.dam.metrics1, method = "spearman", permutations = 999)
kfinedam_C1 # 

kfinedam_C2 <- mantel(dist.k_fine_2, dist.dam.metrics2, method = "spearman", permutations = 999)
OMdam_C2 # 

kfinedam_C3 <- mantel(dist.k_fine_3, dist.dam.metrics3, method = "spearman", permutations = 999)
kfinedam_C3 # 

## partial Mantel tests

kfine_partial1 <- mantel.partial(dist.k_fine_1, dist.env.means1, net_dist1, method="pearson", permutations=999)
kfine_partial1 

kfine_partial2 <- mantel.partial(dist.k_fine_2, dist.env.means2, net_dist2, method="pearson", permutations=999)
kfine_partial2 

kfine_partial3 <- mantel.partial(dist.k_fine_3, dist.env.means3, net_dist3, method="pearson", permutations=999)
kfine_partial3 
#
#
kfine_part1 <- mantel.partial(dist.k_fine_1,  net_dist1, dist.env.means1, method="pearson", permutations=999)
kfine_part1 

kfine_part2 <- mantel.partial(dist.k_fine_2,  net_dist2, dist.env.means2,  method="pearson", permutations=999)
kfine_part2 

kfine_part3 <- mantel.partial(dist.k_fine_3,  net_dist3, dist.env.means3, method="pearson", permutations=999)
kfine_part3 
#
#
kfine_par1 <- mantel.partial(dist.k_fine_1, dist.dam.metrics1, net_dist1, method="pearson", permutations=999)
kfine_par1 

kfine_par2 <- mantel.partial(dist.k_fine_2, dist.dam.metrics2,  net_dist2,  method="pearson", permutations=999)
kfine_par2 

kfine_par3 <- mantel.partial(dist.k_fine_3, dist.dam.metrics3,  net_dist3,  method="pearson", permutations=999)
kfine_par3 

#### Run PLS ####

#Fit PCR model

#no need to centre the variables as this is intrinsic to the pls algorithm 

#subset by season
pls_OMstock_spring <- plsr(log(CO2_C_mg_m2_h+35.92515)~water_pH + water_conductivity_us_cm + DO_mg_L + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + DailyMeanWaterTemp_C + temp_C + substrate_complexity + fine_substrate + mean_depth_cm + percent_open_upstream + percent_open_downstream + masl + dist_to_source_km + dist_to_dam_km + dist_ds_weir_km + percent_open_upstream_weir + No_weirs_dams_up + OM_stock_g_m2, data=dat, scale=TRUE, validation="LOO", method = "oscorespls") #or can use LOOCV

pls_OMstock <- plsr(log(k_day_fine)~ water_pH + water_conductivity_us_cm + DO_mg_L + DO_. + canopy_cover_. + wetted_width_m + discharge_m3_s + mean_velocity_m_s + DailyMeanWaterTemp_C + temp_C + substrate_complexity + fine_substrate + mean_depth_cm +  masl + dist_to_source_km + OM_stock_g_m2, data=dat, scale=TRUE, validation="LOO", method = "oscorespls") 
summary(pls_OMstock)
plot(pls_OMstock)
coef(pls_OMstock) ## Get the PLS coefficients
coef(pls_OMstock, ncomp = 2) #for the 1st two components


validationplot(pls_OMstock, val.type="MSEP")
validationplot(pls_OMstock, val.type="R2")
#Now choose the number of PLS components worth keeping.
#See when RMSE start to increase, 
#Training indicates the percentage of the variance in the response variable explained by the PLS components. 

#Calculate the variable importance 
## VIP returns all VIP values for all variables and all number of components,
## as a ncomp x nvars matrix.
VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

VIP_OMstock_spring <- VIP(pls_OMstock)

#Subset just the first two components
VIP_OMstock_spring <- VIP_OMstock_spring[1:2, ]


# Get the variable names and component names
variables <- rownames(VIP_OMstock_spring)
components <- colnames(VIP_OMstock_spring)

# Create a data frame
VIP_OMstock_spring_df <- data.frame(
  Variables = rep(variables, each = length(components)),
  Components = rep(components, times = length(variables)),
  VIP_Values = as.vector(t(VIP_OMstock_spring))
)

VIP_OMstock_spring_df <- VIP_OMstock_spring_df %>%
  group_by(Components) %>%
  summarize(Mean_VIP = mean(VIP_Values))

# Create the data frame for the coefficients
pls_coef <- coef(pls_OMstock, ncomp = 2)
coefficients <- pls_coef[, , 1]
variables <- rownames(pls_coef) # Get variable names

pls_coef_df <- data.frame(Coefficients = coefficients, Components=variables)

VIP1 <- merge(VIP_OMstock_spring_df, pls_coef_df, by="Components")

#write.csv(VIP1, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/VIP_kcoarse_fall.csv")

#Variables with high influence are identified as VIP > 1, variables with moderate influence as VIP between 1 and 8, and low influence as VIP < 0.8 (Wold et al. 2001).

#plotting predictions
plot(pls_OMstock, ncomp = 5, asp = 1, line = TRUE)

#Plot PLS biplot

plot(VIP_OMstock[, 1], VIP_OMstock[, 2], 
     xlab = "VIP on PC1", ylab = "VIP on PC2", 
     main = "Variable Importance Plot")


#############################################################################
#### ROE obstacles d'ecoulements ####

obst <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Site selection/ObstEcoul.txt")

str(obst)

#SUbset only the validated impoundments

obst_v <- subset(obst, StObstEcou== "Validé" & LbEtOuvrag == "Existant") #all are validated


summary(as.factor(obst_v_NAr$LbHautChut)) #3057 uncategorized, remove these
obst_v_NAr$LbHautChut <- as.factor(obst_v_NAr$LbHautChut)
levels(obst_v_NAr$LbHautChut)

obst_v_NAr <- obst_v[obst_v$LbHautChut != " ", ]

summary(obst_v_NAr$HautChutEt)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.000   0.400   0.920   1.488   1.680 260.000   15620

#How many under 15m? 
sum(obst_v_NAr$HautChutEt < 15, na.rm=T) #38637
sum(obst_v_NAr$HautChutEt > 15, na.rm=T)  #205
(205/38637)*100 #0.5305795% of dams are >15m
1-0.005305795

sum(obst_v_NAr$LbHautChut != "Supérieure ou égale à 10m") #50954
sum(obst_v_NAr$LbHautChut == "Supérieure ou égale à 10m") #594
