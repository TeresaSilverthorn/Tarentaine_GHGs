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

dist <- dist  %>%
  rename(site=Site, , 
         dist_ds_dam_km = "Distance_downstream_dam_5m._km")

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

#### Save as .csv ####

str(dat) #327 × 32 variables

write.csv(dat, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/dat.csv")

###############################################################################

##### Calculate site averages ####

dat_means <- dat %>% 
  dplyr::select(-ID_unique, -transect, -start_time, -end_time) %>% # use select with -var_name to eliminate columns 
  dplyr::group_by(site, campaign, date, flow_state, season) %>% # we group by the two values
  dplyr::summarise(across(.cols = everything(), .fns = mean, .names = "{.col}"))
   
dat_means <- as.data.frame(dat_means)

str(dat_means) #58 obs of 38 vars


###############################################################################
#### Make preliminary plots ####


#### Make a correlation plot ####

#subset just the numeric columns
dat_numeric <- dat %>%
  select_if(is.numeric)

dat_numeric <- na.omit(dat_numeric) #remove NAs

M <-cor(dat_numeric)

corrplot(M, type="upper", order="hclust",
         tl.cex = 0.7, tl.srt = 45,
         col=brewer.pal(n=8, name="RdYlBu"))


## CO2 by site and campaign 
CO2_site<- ggplot(data=dat, aes(x=as.factor(site), y=CO2_C_mg_m2_h, fill=as.factor(season))) +
  geom_bar(stat="identity", position=position_dodge())+  theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) 
CO2_site

CH4_site<- ggplot(data=dat, aes(x=as.factor(site), y=CH4_C_mg_m2_h, fill=as.factor(season))) +
  geom_bar(stat="identity", position=position_dodge())+  theme_bw() +   scale_fill_brewer(palette="Paired")+   theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) 
CH4_site

## GHG by altitude (masl)
CO2_masl <- ggplot(dat, aes(masl, CO2_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_masl

CH4_masl <- ggplot(dat, aes(masl, CH4_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_masl


## GHG and distance to source
CO2_distance <- ggplot(dat, aes(Distance_to_source_km, CO2_C_mg_m2_h)) + geom_point(aes(colour=dist_ds_dam_km), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_distance

dat <- dat %>%
  mutate(position = case_when(
    dist_ds_dam_km < 0 ~ "upstream",
    dist_ds_dam_km > 0 ~ "downstream",
  ))

CO2_distancebyposition <- ggplot(subset(dat_means, CO2_C_mg_m2_h <= 300), aes(Distance_to_source_km, CO2_C_mg_m2_h)) + geom_point(aes(colour=position), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) + geom_smooth(data = subset(dat_means, position %in% c("upstream", "downstream")), method = "lm", se = FALSE, aes(group = position, colour = position), size = 1)
CO2_distancebyposition

CO2_distancebyposition <- ggplot(dat, aes(dist_ds_dam_km, CO2_C_mg_m2_h)) + geom_point(aes(colour=position), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) + geom_smooth(data = subset(dat_means, position %in% c("upstream", "downstream")), method = "lm", se = F, aes(group = position, colour = position), size = 1) + ylim(-10, 500)  + scale_colour_manual(values=c("#f7921e", "#91278e"))
CO2_distancebyposition

CH4_distance <- ggplot(dat, aes(Distance_to_source_km, CH4_C_mg_m2_h)) + geom_point(aes(colour=site), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_distance

CH4_distancebyposition <- ggplot(dat, aes(Distance_to_source_km, CH4_C_mg_m2_h)) + geom_point(aes(colour=position), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) + geom_smooth(data = subset(dat_means, position %in% c("upstream", "downstream")), method = "lm", se = F, aes(group = position, colour = position), size = 1) + scale_colour_manual(values=c("#f7921e", "#91278e"))
CH4_distancebyposition

## GHG and distance to dam
CO2_distancedam <- ggplot(dat, aes(Distance_downstream_dam_5m._km, CO2_C_mg_m2_h)) + geom_point(aes(colour=site), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_distancedam


CH4_distancedam <- ggplot(dat, aes(Distance_downstream_dam_5m._km, CH4_C_mg_m2_h)) + geom_point(aes(colour=site), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_distancedam

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

## GHG and sed OM

CO2_sedOM <- ggplot(dat, aes(sed_OM_percent, CO2_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CO2_sedOM

CH4_sedOM <- ggplot(dat, aes(sed_OM_percent, CH4_C_mg_m2_h)) + geom_point(size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) 
CH4_sedOM

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

LML_fine <- ggplot(dat, aes(site, LML_fine)) +
  geom_bar(stat = "identity", aes(fill = season)) +
  theme_bw() +
  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 12, colour = "black"))
LML_fine

LML_coarse<- ggplot(dat, aes(site, LML_coarse)) +
  geom_bar(stat = "identity", aes(fill = season)) +
  theme_bw() +
  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 12, colour = "black"))
LML_coarse

## OM 

OM_stock <-  ggplot(dat, aes(wetted_width_m, OM_stock_g_m2)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) +   facet_wrap(~ season, ncol = 1)  # Facet by campaign
OM_stock

OM_flux <-  ggplot(dat, aes(wetted_width_m, OM_flux_g_m2_s)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black")) +   facet_wrap(~ season, ncol = 1)  # Facet by campaign
OM_flux


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
OMspring <-ggplot() + 
  geom_sf(data = tarentaine_catch, fill = "lightblue", alpha = 0.3)+
  geom_sf(data = tarentaine) + 
  geom_sf(data = sites_C1, aes(fill = OM_stock_g_m2, shape=flow_state),  size = 4, alpha = 0.3, colour="black") +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_c(option="plasma", na.value="grey40", oob=scales::squish, direction = -1) +
  theme_classic()+ ggtitle("OM stock Fall") + xlab("Longitude") + ylab("Latitude") +  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 8), legend.position ="bottom",legend.text = element_text(size=8)) 


##############################################################################

##### Run Random Forest Model on the data ####

# Goal: to use the sites upstream of the 2 main dams, find the drivers of GHG fluxes and OM, and then using that relationship predict the downstream values, then compare those to the actual impacted values to see what effect the dams have. 


#Should we also try with site means? We need to get rid of ID, Site, campaign, etc. 

#Subset just the upstream sites, TA06, TA11, TA07, TA08, TA12 and higher

upstream_sites <- c("TA06", "TA11", "TA07", "TA08", "TA12", "TA13", "TA14", "TA15", "TA17", "TA20", "TA21", "TA22", "TA24", "TA02R")
#15 sites
downstream_sites <- c("TA04", "TA05", "TA09", "TA10", "TA01", "TA02", "TA03")
#7 sites

# Subset 
dat_rf <- dat %>% 
 # filter(site %in% upstream_sites)  %>%
  select(-date, -site, -start_time, -end_time, -ID_unique,  -campaign, -transect, -lat, -lon, -LML_coarse, -LML_fine, -k_day_coarse,  -k_day_fine, -CH4_C_mg_m2_h, -flow_state, -season, -cobble, -boulder, -pebble_gravel)

str(dat_upstream)
str(dat_rf)


# Convert the tibble to a data frame
dat_upstream <- as.data.frame(dat_upstream)
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
bestmtry <- tuneRF(data_train,data_train$CO2_C_mg_m2_h,stepFactor = 1.2, improve = 0.01, trace=T, plot= T)  #7

#create RF model: 
model <- randomForest(CO2_C_mg_m2_h~.,data= data_train, ntree = 1000, mtry = 7, importance = TRUE) #need to determine what is the best value for mtry and ntree
model 

plot(model)

imp <- importance(model)  # returns the importance of the variables

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

### Visualize variable importance ----------------------------------------------

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


#### Try structural equation modeling SEM ####

library(lavaan)
library(semPlot)

dat <-as.data.frame(dat)

row.names(dat) <- NULL  # Remove existing row names
row.names(dat) <- 1:nrow(dat)  # Assign new row names

#Fit a test model
model <-'CO2_C_mg_m2_h ~ sed_OM_percent + DO_. + water_pH 
CH4_C_mg_m2_h ~ sed_OM_percent + DO_. + water_pH 
CO2_C_mg_m2_h ~ dist_ds_dam_km 
CH4_C_mg_m2_h ~ dist_ds_dam_km 
sed_OM_percent  ~ dist_ds_dam_km 
DO_. ~ dist_ds_dam_km 
water_pH ~ dist_ds_dam_km 
'

fit <- cfa(model, data = dat)

summary(fit)

node_labels <- c("CO2", "CH4",  "sed_OM", "DO", "water_pH", "dist_dam")

# Create path diagram with adjusted settings
semPaths(fit, intercept = FALSE, whatLabel = "est",
         residuals = FALSE, exoCov = FALSE, 
         nodeLabels = node_labels,
         style='lisrel', edge.label.cex = 1.3,
         layout = "tree", sizeMan = 7)

#### Find best LMM and make predictions ####

# Can make a model explaining GHG fluxes using all of the sites upstream of the two dams, then predict the CO2/CH4 values for downstream of the dams, and see if they differ from the actual measured values

#Maybe first you can use RF to determine top drivers with a subset of the data set. 
# sedOM%, DO%, water pH, water conductivity, wetted width, distance ds dam



