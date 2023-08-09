#### R script to calculate the organic matter stock and flux across the Tarentaine River, over 3 campaigns ####

#Load necessary packages
library(dplyr)
library(readxl)
library(ggplot2)

#Load the ancillary data

ancil_dat <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Ancillary data/Data_entry_Tarentaine_2022_2023-08-07.csv", header=T)

str(ancil_dat)  #327 obs. of  43 variables


#Load in the OM flux and stock data

OM <-read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/OM_flux_stock/OMflux_stock_weights_Tarentaine2022.csv")

str(OM)

#### Calculate ash free dry mass (AFDM) of the stock and flux samples ####

# AFDM = DM - ash

#But because we used a subset, we can calculate the percent of OM in the burnt sample, and multiply that by the dry leaves
#First calculate the % of organic matter in the subsamples: 
# %OM = (subsample dry mass - subsample ash dry mass ) / subsample dry mass * 100

OM$OM_percent <-( ((OM$DM_to_burn.pan_g - OM$ceramic_pan_g) - (OM$ash.pan_g - OM$ceramic_pan_g)) /  (OM$DM_to_burn.pan_g - OM$ceramic_pan_g) ) *100

#dealing with the wet weight subsets, calculate percent moisture of original sample
#subtract the dry weight from the wet weight, divide the result by the wet weight, then multiply by 100.
OM$percent_moist <- ( (OM$wet_weight_g_sub - (OM$oven_dry_weight_g...pan-OM$pan_weight_g) ) / OM$wet_weight_g_sub) *100

#Subtract the percent moisture from the original wet sample
OM$OM_dry_tot <- OM$wet_weight_g* (1 - (OM$percent_moist/100) )

#Multiple the oven dry mass by the OM percentage, taking into accounth the cases where we subset the wet volume
OM <- OM %>%
  mutate(AFDM = ifelse(!is.na(OM_dry_tot),
                       (OM_dry_tot) * (OM_percent / 100),
                       (oven_dry_weight_g...pan - pan_weight_g) * (OM_percent / 100)))


# Calculate the AFDM of OM stock on a per area basis

#0.105625 m2 quadrat for stock (32.5 cm by 32.5 cm)

#OM$OM_stock_g_m2 <- OM$AFDM / (OM$No_samples_in_reach*0.11)

OM_stock <- subset(OM, type=="stock")

OM_stock <- OM_stock %>%
  mutate(OM_stock_g_m2 = if_else(type == "stock", AFDM / (No_samples_in_reach * 0.105625), NA_real_))

str(OM_stock)
#59 obs. of  3 variables 

#why are there two values for TA05 in C2 and TA10 C1? These were two sub samples from the ziplocks 
#no TA11 for C3 is normal, the flow was too high
OM_stock <- OM_stock %>%
group_by(campaign, site) %>%
  summarize(OM_stock_g_m2 = mean(OM_stock_g_m2, na.rm = TRUE))

str(OM_stock) #57 obs of 3 vars

#Calculate AFDM of flux

# can chose to calculate flux as g/m2/h or g/h depending on whether you use velocity or discharge

#0.086125 m2 net opening (32.5 cm by 26.5 cm)

#Flux should be associated with T6, or whatever is the last transect 

#Need to account for the time the net was out

#Need to combine ancil dat and OM
ancil_dat_sub <- ancil_dat %>%
  select(site, campaign, transect, mean_velocity_m_s, discharge_m3_s, OM_flux_net_time_mins)

# Group the data by 'site' and select the row with the highest 'transect' number, because the nets were placed at the most upstream of the reach, and thus this transect is most representative in terms of velocity and discharge

ancil_dat_sub <- ancil_dat_sub %>%
  group_by(site, campaign) %>%
  filter(transect == max(transect)) %>%
  ungroup()

OM_flux <- subset(OM, type=="flux")

OM_flux <- merge(ancil_dat_sub, OM_flux, by = c("site", "campaign"))


# Calculate OM flux 

#With the help of Jerome leCoz, this equation is: 
# Flux (g/m2/s) = AFDM (g) / (area (m2) * time (s)  ) 

OM_flux <- OM_flux %>%
  mutate(OM_flux_g_m2_s = AFDM/ ((No_samples_in_reach*0.086125) * (OM_flux_net_time_mins*60))    )

# Calculate OM water column concentration, accounting for the discharge 
#F (g/s) = OM_flux_g_m2_s * the surface area of the transect m2 (or alternatively by the width of the transect)
#Then divide that value by the discharge, C (g/m3) = F/Q

#####################################################


#### plot ####

flux <- ggplot(data=OM_flux, aes(x=as.factor(site), y=OM_flux_g_m2_s, fill=as.factor(campaign))) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+   theme_minimal()
flux


stock <- ggplot(data=OM_stock, aes(x=as.factor(site), y=OM_stock_g_m2, fill=as.factor(campaign))) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+  scale_y_log10()  +
  theme_minimal()
stock


#Merge the data into one file

OM_stock_sub <- OM_stock %>%
  select(campaign, site, OM_stock_g_m2) %>% #subset relevant columns
  na.omit()
  
str(OM_stock_sub) #56 obs of 3 vars

OM_flux_sub <- OM_flux %>%
  select(campaign, site, OM_flux_g_m2_s)

str(OM_flux_sub) #56 obs. of  3 variables

OM_stock_flux <- merge(OM_stock_sub, OM_flux_sub, by = c("campaign", "site"))

str(OM_stock_flux) #56 obs of 4 vars

#Save as csv

write.csv(OM_stock_flux, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/OM_flux_stock/OM_stock_flux.csv")
