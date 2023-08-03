#### Script to calculate leaf litter decomposition of fine and coarse mesh leaf packs along the Tarentaine river network across 2 campaigns ####

# Load the necessary packages: 
library(dplyr)


##########################

# Input the leaf pack mass data

lp <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Leaf packs/Leaf_pack_weights_Tarentaine_all_campaigns_one_file.csv")

str(lp) #363 obs. of  15 variables

# Exclude any leaf pack that were found open or out of the water

lp <- lp %>% 
  mutate(LP_open_or_dry = na_if(LP_open_or_dry, ""))   %>%     #First make blank cells NA
  filter(is.na(LP_open_or_dry))  %>% 
  select(-Pan_ID, -Al_pan_g ) #Drop non-useful columns

str(lp) #338 obs. of  15 variables, so we lost 25 leaf packs


#### Calculate ash free dry mass remaining ####

#First calculate the % of organic matter in the subsamples: 

# %OM = (sub-sample dry mass - sub-sample ash dry mass ) / sub-sample dry mass * 100
lp$OM_percent <-( ((lp$DM_to_burn.pan_g - lp$Ceramic_pan_g) - (lp$Ash.pan_g - lp$Ceramic_pan_g)) /  (lp$DM_to_burn.pan_g - lp$Ceramic_pan_g) ) *100


#Calculate the AFDM by multiplying the dry sample by the %OM content
lp$final_AFDM <- lp$Final_dry_weight_g * (lp$OM_percent/100)

#Need the OM percentage of the initial leaves (can use the handling loss bags)

Initial <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Leaf packs/Initial_leaves_ashed.csv")

str(Initial)

#Calculate teh AFDM of the initial leaves

Initial$OM_percent <-( ((Initial$DM_to_burn.pan_g - Initial$Ceramic_pan_g) - (Initial$Ash.pan_g - Initial$Ceramic_pan_g)) /  (Initial$DM_to_burn.pan_g - Initial$Ceramic_pan_g) ) *100

#The average OM content of the initial leaves: 
mean(Initial$OM_percent, na.rm=T) #94.53203

#Calculate the decomposition rate (k)

#Decomposition rate (k) is the proportion of litter mass loss (LML) per degree day (dd), to account for differences in temperature across sites; 

#LML = [initial AFDM (g) – final AFDM (g)]/initial AFDM (g), where initial AFDM was previously corrected by leaching, drying and ash content, which were estimated in the laboratory
lp$LML <- ( (lp$Initial_dry_weight_g * 0.9453203) - lp$final_AFDM) / lp$Initial_dry_weight

#Get rid of the LML percentages which are negatives, super weird, and not sure where this comes form, type or a rock in the bag?

lp_site_means <- lp %>%
  mutate(LML = ifelse(LML < 0, NA, LML))  %>%
  group_by(Campaign, Site) %>%
  mutate(mean_LML = mean(LML, na.rm = TRUE))

#plot

lp_plot <- ggplot(lp_site_means, aes(x=Site, y=LML, fill=as.factor(Campaign))) + 
  geom_bar(stat="identity", position=position_dodge())
lp_plot






