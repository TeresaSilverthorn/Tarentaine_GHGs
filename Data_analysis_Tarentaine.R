#### Data analysis and visualization of the OM and GHG dynamics in a river network (Tarentaine) impacted by damming ####

#Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(RColorBrewer)

#Load and merge the 1. ancillary data, 2. GHG flux data, 3. OM stock and flux data, 4. leaf litter decomposition data, 5. daily temperature data, 6. Sediment OM content. Check for quality control 

#### 1. Load ancillary data ####

ancil_dat <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Ancillary data/Data_entry_Tarentaine_2022_2023-08-07.csv", header=T) #Make sure it's the latest version from the Google Drive

str(ancil_dat) #327 obs of 42 vars

#Drop the useless columns, get the depth and velocity average per transect
ancil_dat <- ancil_dat %>%
  select(-CO2_start_offset_mins, -CO2_end_offset_mins, -CH4_start_offset_mins, -CH4_end_offset_mins, -ibutton_air, -ibutton_water, -X, -chamber_area, -chamber_volume, -OM_flux_net_time_mins) %>%
  rowwise() %>%
  select(-starts_with("depth"), -starts_with("velocity"))

#Add a "ID_unique" column with this format : "2022-05-30_TA01_1"




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

str(decomp) #323 obs. of  20 variables

#reshape the data so coarse and fine mesh LML become their own columns per site/campaign
decomp <- decomp %>%
  pivot_wider(names_from = Mesh_size, values_from = LML_perday) %>%
  rename(LMLperday_coarse = Large, LMLperday_fine = Small) %>%
  group_by(Campaign, Site) %>%   #to solve for the NA in alternating rows, group
  summarize(LMLperday_coarse = mean(LMLperday_coarse, na.rm = TRUE),   #then take the mean
            LMLperday_fine = mean(LMLperday_fine, na.rm = TRUE)) %>%
  mutate(across(everything(), ~replace(., is.nan(.), NA)))    #replace NaN with NA

str(decomp) #58 obs of 4 vars

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

#Remove the weird outlying value of sed_OM_percent "14.5095801"

sed_OM <- sed_OM %>%
  filter(sed_OM_percent <= 14)

#Rename columns to match
sed_OM <- sed_OM %>%
  rename(site = Site,
         campaign = Campaign)

#### Merge the dataframes together ####

#Rename columns to match
decomp <- decomp %>%
  rename(site = Site,
       campaign = Campaign)

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


#remove unnecessary columns

str(dat)

dat <- dat %>%
  select(-X.x, -X.y, -Notes) 

#########################################################################

#### Fill any missing data ####

## add the substrate type to all of the campaigns

columns_to_fill <- c("X._bedrock", "X._boulder", "X._cobble", "X._pebble_gravel", "X._fine_substrate")

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
  X._bedrock= 0,
  X._boulder = 70,
  X._cobble = 20,
  X._pebble_gravel = 10, 
  X._fine_substrate = 0 )

fill_values_TA06 <- data.frame(
  site = "TA06",
  X._bedrock = 0,
  X._boulder = 15,
  X._cobble = 65,
  X._pebble_gravel = 15,
  X._fine_substrate = 5 )

# Combine fill values for both sites
all_fill_values <- bind_rows(fill_values_TA02, fill_values_TA06)

# Replace "TA02" and "TA06" rows with fill values
dat <- dat %>%
  left_join(all_fill_values, by = "site") %>%
  mutate(
    X._bedrock = if_else(is.na(X._bedrock.y), X._bedrock.x, X._bedrock.y),
    X._boulder = if_else(is.na(X._boulder.y), X._boulder.x, X._boulder.y),
    X._cobble = if_else(is.na(X._cobble.y), X._cobble.x, X._cobble.y),
    X._pebble_gravel = if_else(is.na(X._pebble_gravel.y), X._pebble_gravel.x, X._pebble_gravel.y),
    X._fine_substrate = if_else(is.na(X._fine_substrate.y), X._fine_substrate.x, X._fine_substrate.y)
  ) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>%
  arrange(campaign, site, transect)

###############################################################################

#### Save as .csv ####

write.csv(dat, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/dat.csv")

###############################################################################

#### Make preliminary plots ####

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

#### Make a correlation plot ####

#subset just the numeric columns
dat_numeric <- dat %>%
  select_if(is.numeric)

dat_numeric <- na.omit(dat_numeric) #remove NAs

M <-cor(dat_numeric)

corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))


## GHG and width/depth

CO2_width <- ggplot(dat, aes(wetted_width_m, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) +   facet_wrap(~ season, ncol = 1)  # Facet by campaign
CO2_width

CO2_depth <- ggplot(dat, aes(mean_depth_cm, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_depth

CH4_width <- ggplot(dat, aes(wetted_width_m, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) +   facet_wrap(~ season, ncol = 1)  # Facet by campaign
CH4_width

CH4_depth <- ggplot(dat, aes(mean_depth_cm, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_depth

## GHG temperature

CO2_airtemp <- ggplot(dat, aes(temp_C, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_airtemp

CH4_airtemp <- ggplot(dat, aes(temp_C, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_airtemp

CO2_watertemp <- ggplot(dat, aes(DailyMeanWaterTemp_C, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_watertemp

CH4_watertemp <- ggplot(dat, aes(DailyMeanWaterTemp_C, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_watertemp

## GHG and Velocity /discharge

CO2_veloc <- ggplot(dat, aes(mean_velocity_m_s, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_veloc

CH4_veloc <- ggplot(dat, aes(mean_velocity_m_s, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_veloc

CO2_discharge <- ggplot(dat, aes(discharge_m3_s, CO2_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_discharge

CH4_discharge <- ggplot(dat, aes(discharge_m3_s, CH4_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_discharge

## GHG water chemistry

CO2_DO <- ggplot(dat, aes(DO_mg_L, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_DO

CH4_DO <- ggplot(dat, aes(DO_mg_L, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_DO

CO2_pH <- ggplot(dat, aes(water_pH, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_pH

CH4_pH <- ggplot(dat, aes(water_pH, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_pH

CO2_cond <- ggplot(dat, aes(water_conductivity_us_cm, CO2_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_cond

CH4_cond <- ggplot(dat, aes(water_conductivity_us_cm, CH4_C_mg_m2_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_cond

## GHG and canopy cover

CO2_canopy <- ggplot(dat, aes(canopy_cover_., CO2_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) + facet_wrap(~ season, ncol = 1) 
CO2_canopy

CH4_canopy <- ggplot(dat, aes(canopy_cover_., CH4_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black"))  +facet_wrap(~ season, ncol = 1) 
CH4_canopy

## Decomp 

LML_fine <- ggplot(dat, aes(site, LMLperday_fine)) +
  geom_bar(stat = "identity", aes(fill = season)) +
  theme_bw() +
  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 12, colour = "black"))
LML_fine

LML_coarse<- ggplot(dat, aes(site, LMLperday_coarse)) +
  geom_bar(stat = "identity", aes(fill = season)) +
  theme_bw() +
  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 12, colour = "black"))
LML_coarse

## OM 

OM_stock <-  ggplot(dat, aes(wetted_width_m, OM_stock_g_m2)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) +   facet_wrap(~ season, ncol = 1)  # Facet by campaign
OM_stock

OM_flux <-  ggplot(dat, aes(wetted_width_m, OM_flux_g_m3_h)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) +   facet_wrap(~ season, ncol = 1)  # Facet by campaign
OM_flux



