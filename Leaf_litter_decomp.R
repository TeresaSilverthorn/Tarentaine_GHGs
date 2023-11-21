#### Script to calculate leaf litter decomposition of fine and coarse mesh leaf packs along the Tarentaine river network across 2 campaigns ####

# Load the necessary packages: 
library(dplyr)


##########################

# Input the leaf pack mass data

lp <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Leaf_packs/Leaf_pack_weights_Tarentaine_all_campaigns_one_file.csv")

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

#Determine the average handling loss and subtract from final dry weight
#Calculate the percent handling loss

#Subset handling loss 
handling_loss <- subset(lp, Site=="Handling loss")

str(handling_loss)

handling_loss$percent_loss <-( ((handling_loss$Initial_dry_weight_g) - (handling_loss$Final_dry_weight_g)) /  (handling_loss$Initial_dry_weight_g) ) * 100
#In all of the cases we actually gained weight somehow.... so can assume a negligible handling loss


#Calculate the AFDM by multiplying the dry sample by the %OM content
lp$final_AFDM <- lp$Final_dry_weight_g * (lp$OM_percent/100)

#Need the OM percentage of the initial leaves (can use the handling loss bags)

Initial <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Leaf_packs/Initial_leaves_ashed.csv")

str(Initial)

#Calculate the AFDM of the initial leaves
Initial$OM_percent <-( ((Initial$DM_to_burn.pan_g - Initial$Ceramic_pan_g) - (Initial$Ash.pan_g - Initial$Ceramic_pan_g)) /  (Initial$DM_to_burn.pan_g - Initial$Ceramic_pan_g) ) *100

#The average OM content of the initial leaves: 
mean(Initial$OM_percent, na.rm=T) #94.53203

#Calculate the decomposition rate (k)

#Decomposition rate (k) is the proportion of litter mass loss (LML) per degree day (dd), to account for differences in temperature across sites; 

#LML = [initial AFDM (g) – final AFDM (g)]/initial AFDM (g), where initial AFDM was previously corrected by leaching, drying and ash content, which were estimated in the laboratory per Boyero et al 2021
#this is really a percent of the LML
lp$LML <- ( (lp$Initial_dry_weight_g * 0.9453203) - lp$final_AFDM) / lp$Initial_dry_weight
 
##### Divide by degree days using water temperature #####

water_temp <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Daily_mean_water_temp.csv")

# Convert Date to Date type
water_temp$Date <- as.Date(water_temp$Date)

# Add the campaign column based on date ranges
water_temp <- water_temp %>%
  mutate(
    campaign = case_when(
      Date >= as.Date("2022-05-30") & Date <= as.Date("2022-06-22") ~ 1,
      Date >= as.Date("2022-07-18") & Date <= as.Date("2022-08-10") ~ 2,
      Date >= as.Date("2022-10-10") & Date <= as.Date("2022-11-04") ~ 3,
      TRUE ~ NA_integer_    )  ) %>%
  filter(!is.na(campaign))

#Make lp Date as.Date too 
lp$Date_in <- as.Date(lp$Date_in)
lp$Date_out <- as.Date(lp$Date_out)

# Calculate the number of days between Date_out and Date_in
lp$Duration_days <- as.numeric(lp$Date_out - lp$Date_in)

mean(lp$Duration_days, na.rm=T)
sd(lp$Duration_days, na.rm=T)

lp <- lp %>%
rename(site = Site, campaign = Campaign) #rename to harmonize

lp_summarized <- lp %>%
  group_by(site, campaign, Date_in, Date_out) %>%
  filter(!is.na(Date_in) & !is.na(Date_out)) %>%
  summarize()

unique_campaigns <- unique(lp_summarized$campaign)
result_list <- list()

for (campaign in unique_campaigns) {
  campaign_lp <- lp_summarized %>% filter(campaign == campaign)
  
  campaign_degree_days <- campaign_lp %>%
    left_join(water_temp, by = c("site", "campaign")) %>%
    filter(Date >= Date_in & Date <= Date_out) %>%
    group_by(site, campaign, Date_in, Date_out) %>%
    summarize(
      Degree_Days = sum(ifelse(DailyMeanWaterTemp_C > 0, DailyMeanWaterTemp_C, 0))
    )
  
  result_list[[as.character(campaign)]] <- campaign_degree_days
}

# Combine the results from each campaign
degree_days <- bind_rows(result_list)
head(degree_days)

#Merge degree days with lp 

#select useful columns
degree_days <- degree_days %>%
  ungroup%>%
  select(site, campaign, Degree_Days)

# merge DD to the original lp dataframe
lp <- merge(lp, degree_days, by = c("site", "campaign")) 


#Check

sum(water_temp$DailyMeanWaterTemp_C[water_temp$site == "TA24" & water_temp$Date >= as.Date("2022-10-14") & water_temp$Date <= as.Date("2022-11-03")]) #209.2917

sum(water_temp$DailyMeanWaterTemp_C[water_temp$site == "TA01" & water_temp$Date >= as.Date("2022-07-19") & water_temp$Date <= as.Date("2022-08-08")]) #396.3333


#Calculate decomposition rate corrected by degree days 
#k = -ln (M final / M initial ) / degree days

lp$k_day <- -log( lp$final_AFDM / lp$Initial_dry_weight_g * 0.9453203) / lp$Duration_days  #Initial weight multiplied by OM content, final AFDM is already multiplied by OM content

lp$k_dday <- -log( lp$final_AFDM / lp$Initial_dry_weight_g * 0.9453203) / lp$Degree_Days

#Calculate site means
  
#Get rid of the LML percentages which are negatives, super weird, and not sure where this comes form, type or a rock in the bag?
lp_site_means <- lp %>%
  mutate(LML = ifelse(LML < 0, NA, LML))  %>%  #remove the 3 negatives, data entry error
  filter(site != "Handling loss" & site != "") %>%   #remove the handling loss and blank site rows
  group_by(campaign, site, Mesh_size) %>%   #group
  summarize(LML = mean(LML, na.rm = TRUE),
            k_day = mean(k_day, na.rm = TRUE),
            k_dday = mean(k_dday, na.rm = TRUE),
            )   #get mean per the 3 replicates

#plot
lp_plot <- ggplot(lp_site_means, aes(x=site, y=k_day, fill=as.factor(campaign))) + 
  geom_bar(stat="identity", position=position_dodge())
lp_plot


#### Save as .csv ####


write.csv(lp_site_means, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Leaf_packs/Leaf_pack_decomposition.csv")


