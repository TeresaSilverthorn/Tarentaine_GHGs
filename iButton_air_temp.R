#### Script for iButton data for the Tarentaine catchment ####

# 1. Extract air temperature data from the iButton dataloggers to get an average air temperature per GHG measurement to correct for in the gasflux calculation

# 2. Calculate degree days for decomposition correction

# 3. Calculate daily mean water temperature




## load necessary packages ##
library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(survival)
library(scales)
library(tidyr)
library(purrr)
library(ggpubr)
library(broom)



#### 1. Air temperature for GHG measurements ####

#iButton data can be found here
# C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons 

setwd("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons")
#
# Note that for C1 we had an air iButton at each site (17 datasets, so missing 3), and for C2 and C3 we had a single iButton on the Picarro
#
#
# Note when missing data, the datalogger from the next closest site was used. NA in the iButton ID column denotes missing datalogger, the third "missing" one is becaus ethe same ibutton data is used for TA10 and TA09

#missing:   76_5200000031CD6D21_072922.csv  (TA12)        81_4200000031FB2121_072922.csv (TA04)

# The closest sites are 
# TA12: TA13 is super close by
# TA04: TA05 is super close by

#
#
# ibutton information (file name, type, site, ibutton ID) can be found here: 
ibutton_master <- read.csv("iButton_programming_Tarentaine_2022 - Campaign1.csv", header=T)

head(ibutton_master)

#Add a new column for the file name
ibutton_master <- ibutton_master %>%
  mutate(
    file_name = paste0(ID, "_", Device_address, "_", "072922.csv" ) 
    ) %>%
  select(file_name, everything()) %>%    #reorder
  filter(habitat == "Air")   #Keep only the air iButtons


str(ibutton_master)  #20 obs. of  14 variables

# replace the missing sites with the next closest site data: 

#missing:   76_5200000031CD6D21_072922.csv  (TA12 to TA13)        81_4200000031FB2121_072922.csv (TA04 to TA05)

ibutton_master <- ibutton_master %>%
  mutate(file_name = replace(file_name, file_name == "76_5200000031CD6D21_072922.csv", "87_1D00000031FE6A21_072922.csv"))%>%
  mutate(file_name = replace(file_name, file_name == "81_4200000031FB2121_072922.csv", "84_A900000031FA3D21_072922.csv"))

head(ibutton_master)

########################################################################

#Load in the ancillary GHG data with the start and end times of the GHG measurements
ancil_dat <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Ancillary data/Data_entry_Tarentaine_2022_2023-07-25.csv", header=T)

str(ancil_dat) #327obs

#add an ID column  
ancil_dat <- ancil_dat %>%
  mutate(ID_unique = paste(date, site, transect, sep = "_")) %>%
  select(ID_unique, everything()) #make this ID_unique column come first

#Combine the date and time #Note that here we are using the actual noted watch time, and not the time corrected for the Picarro drift, as the latter will should correspond better to the iButton temperature
ancil_dat <- ancil_dat %>%
  mutate(
    start_time = ifelse(nchar(start_time) == 5, paste0(start_time, ":00"), start_time),
    end_time = ifelse(nchar(end_time) == 5, paste0(end_time, ":00"), end_time),
    datetime_start = as.POSIXct(paste(date, start_time), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris"),
    datetime_end = as.POSIXct(paste(date, end_time), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris")
  ) %>%
  select(ID_unique, campaign, date, site, flow_state, datetime_start, datetime_end, everything()) %>% #reorder
  mutate(date = as.POSIXct(date, format = "%Y-%m-%d", tz = "Europe/Paris")) #make date as.POSIXct too

#Merge the data logger master sheet  with the ancillary data
#subset just the first campaign
ancil_dat_C1 <- subset(ancil_dat, campaign=="1")

air_temp <- merge(ancil_dat_C1,ibutton_master,by="site")

str(air_temp)   # 	119 obs. of  56 variables

#################################################################################################
### Use a loop to select the correct datalogger file and choose the closest air temperature #####

#try with just a subset
#air_temp <- subset(air_temp, file_name=="11_0100000031F65021_072922.csv" | file_name=="20_5F00000031EEA221_072922.csv")


#create a matrix with one column for the temperature associated with each unique IF
matTemp<-matrix(rep(NA,nrow(air_temp)*1),nrow=nrow(air_temp),ncol=1,dimnames=(list(air_temp$ID_unique ,c("Temp")))) 


## Loop for finding nearest temperature to each start time (ibuttons) 

for(i in 1:nrow(air_temp)){
  I<-air_temp[i,]
  ID<-I$file_name
  Dat <- fread(paste("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Air_iButtons/",ID,sep=""), header=TRUE)
  colnames(Dat) <- c("DateTime", "Unit", "Temp_C")
  Dat$date <- as.POSIXct(Dat$DateTime, format = "%d/%m/%y %I:%M:%S %p", tz="Europe/Paris", origin = "1970-01-01")
  Dat <- Dat %>% drop_na()
  I$datetime_start<- as.POSIXct(I$datetime_start, format="%Y-%m-%d %H:%M", tz="Europe/Paris", origin = "1970-01-01")
  
  setDT(Dat)            ## convert to data.table by reference
  setDT(I)            ## same
  
  I[, date := datetime_start] 
  setkey(I, datetime_start)    ## set the column to perform the join on
  setkey(Dat, date)    ## same as above
  
  ans = Dat[I, roll=Inf] ## perform rolling join
  
  matTemp[i,"Temp"]<-ans$Temp_C
  
}

write.table(matTemp,file="C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/ibutton_air_temp_C1.csv") 


### Read in the data and check for any outliers

#first you will need to read in the CSV's using fread because of the weird formatting of write.table

air.temp <- fread("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/ibutton_air_temp_C1.csv") #For C1 

colnames(air.temp) <- c("ID_unique", "Air_temp")


#Merge with ancil_dat

ancil_dat_temp <- merge(ancil_dat_C1, air.temp, by = "ID_unique")

#plot

May30 <- ggplot(data = subset(ancil_dat_temp, month(datetime_start) == 5 & day(datetime_start) == 30), aes(x = datetime_start, y = Air_temp)) + geom_point() + scale_x_datetime(date_labels = "%H", breaks = "1 hour")
May30 

May31 <- ggplot(data = subset(ancil_dat_temp, month(datetime_start) == 5 & day(datetime_start) == 31), aes(x = datetime_start, y = Air_temp)) + geom_point() + scale_x_datetime(date_labels = "%H", breaks = "1 hour")
May31 

Jun01 <- ggplot(data = subset(ancil_dat_temp, month(datetime_start) == 6 & day(datetime_start) == 1), aes(x = datetime_start, y = Air_temp)) + geom_point() + scale_x_datetime(date_labels = "%H", breaks = "1 hour")
Jun01 

Jun02 <- ggplot(data = subset(ancil_dat_temp, month(datetime_start) == 6 & day(datetime_start) == 2), aes(x = datetime_start, y = Air_temp)) + geom_point() + scale_x_datetime(date_labels = "%H", breaks = "1 hour")
Jun02 

Jun03 <- ggplot(data = subset(ancil_dat_temp, month(datetime_start) == 6 & day(datetime_start) == 3), aes(x = datetime_start, y = Air_temp)) + geom_point() + scale_x_datetime(date_labels = "%H", breaks = "1 hour")
Jun03 


###############################################################################

#### For campaign 2 and 3 ####

#We used a single iButton attached to the Picarro in these cases, so need to watch out for artefacts from storing the material in the car!

#Read in the iButton for C2
temp_C2 <- fread("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/C2_Picarro_iButton/55_460000003200B421_072722.csv", header=T)

#rename columns
colnames(temp_C2) <- c("date_time", "unit", "temp_C")

#change time format
temp_C2$date_time <- as.POSIXct(temp_C2$date_time, format = "%d/%m/%y %I:%M:%S %p")


#Read in the iButton for C3
temp_C3 <- fread("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/C3_Picarro_iButton/64_9A00000031D48F21_101822.csv", header=T)

#rename columns
colnames(temp_C3) <- c("date_time", "unit", "temp_C")

#change time format
temp_C3$date_time <- as.POSIXct(temp_C3$date_time, format = "%d/%m/%y %I:%M:%S %p")

#################

#Merge 

temp_C2_C3 <- bind_rows(temp_C2, temp_C3)

#################

#################################################################################################
### Use a loop to select the correct datalogger file and choose the closest air temperature #####

#Subset campaigns 2 and 3 
ancil_dat_C2_C3 <- subset(ancil_dat, campaign=="2" | campaign=="3")

# Convert both dataframes to data.table for efficient rolling join
setDT(ancil_dat_C2_C3)
setDT(temp_C2_C3)

# Perform the rolling join to find the closest datetime
test <- temp_C2_C3[ancil_dat_C2_C3, roll = "nearest", on = .(date_time = datetime_start)]

#Make some plots for quality control to see how this worked
plotC1 <- ggplot(data = subset(temp_C2_C3, month(date_time) == 7), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%d %H", breaks = "12 hour")
plotC1 

plot_testC1 <- ggplot(data = subset(test, month(date_time) == 7), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%d %H", breaks = "12 hour")
plot_testC1 

plotC2 <- ggplot(data = subset(temp_C2_C3, month(date_time) == 10), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%d %H", breaks = "12 hour")
plotC2 

plot_testC2 <- ggplot(data = subset(test, month(date_time) == 10), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%d %H", breaks = "12 hour")
plot_testC2 

Jul18 <- ggplot(data = subset(temp_C2_C3, month(date_time) == 7 & day(date_time) == 18), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%H:%m", breaks = "1 hour")
Jul18 

Jul18 <- ggplot(data = subset(temp_C2_C3, month(date_time) == 7 & day(date_time) == 18), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%H:%m", breaks = "1 hour")
Jul18 

###############################################

#Combine the temperature data from all campaigns into one dataframe 

C2C3 <-  test[, c("ID_unique", "temp_C")]

colnames(air.temp) <- c("ID_unique", "temp_C")

C1 <- air.temp

#merge 

air_temp_all <- bind_rows(C1, C2C3)

#save as csv and load into GHG script

write.csv(air_temp_all,"C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/ibutton_air_temp.csv") 

##########################################################################
#### 2. Calculate degree days for decomposition correction ####

#### Run a loop for the air temp data but get the daily mean, starting with C1 where we had an iButton per site ####

# Initialize an empty data frame
all_data <- data.frame(DateTime = character(),
                         Unit = character(),
                         Temp_C = numeric(),
                         site = character(),
                         stringsAsFactors = FALSE)

for (i in 1:nrow(ibutton_master)) {
  site <- ibutton_master$site[i]
  file_name <- ibutton_master$file_name[i]
  file_path <- paste("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Air_iButtons/", file_name, sep = "")
  
  # Load the data from the file
  Dat <- fread(file_path, header = TRUE)
  colnames(Dat) <- c("DateTime", "Unit", "Temp_C")
  Dat$DateTime <- as.POSIXct(Dat$DateTime, format = "%d/%m/%y %I:%M:%S %p", tz = "Europe/Paris", origin = "1970-01-01")
  Dat$file_name <- file_name
  
  # Add the 'site' column and fill it with the site value
  Dat$site <- site
  
  Dat <- Dat %>% drop_na()
  
  # Append the data to the all_data data frame
  all_data <- rbind(all_data, Dat)
}

str(all_data) #29411 obs. of  6 variables

levels(as.factor(all_data$site)) #All sites except TA02R

# Now calculate the daily average temperature for each Site

daily_ave_air_temp_C1 <- all_data %>%
  group_by(site, Date = as.Date(DateTime)) %>%
  summarise(DailyMeanTemp_C = mean(Temp_C, na.rm = TRUE))

# Plot to check

# Convert the 'Date' column to POSIXct format for plotting
daily_ave_air_temp_C1$Date <- as.POSIXct(daily_ave_air_temp_C1$Date, format = "%Y-%m-%d")

TA01 <- ggplot(data = subset(daily_ave_air_temp_C1, site=="TA01"), aes(x = Date, y = DailyMeanTemp_C)) + geom_point() + scale_x_datetime(date_labels = "%m/%d", breaks = "5 days")
TA01

TA05 <- ggplot(data = subset(daily_ave_air_temp_C1, site=="TA05"), aes(x = Date, y = DailyMeanTemp_C)) + geom_point() + scale_x_datetime(date_labels = "%m/%d", breaks = "5 days")
TA05

TA24 <- ggplot(data = subset(daily_ave_air_temp_C1, site=="TA24"), aes(x = Date, y = DailyMeanTemp_C)) + geom_point() + scale_x_datetime(date_labels = "%m/%d", breaks = "5 days")
TA24

#Save as .csv
write.csv(daily_ave_air_temp_C1, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/daily_ave_air_temp_C1.csv")


#Looks good!

# Now calculate the daily average air temp for C2 and C3, where we only had 1 air temperature iButton so we can't really see site differences, check weather stations maybe



################################################################################

#### 3. Calculate the average daily water temperature ####

#### Check in C3, the ID 20 iButton was noted for TA24 and TA13, check the data, as 24 should be colder ####

Ib20 <- fread("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Leaf_packs/Campaign 3/20_5F00000031EEA221_120122.csv")

#rename columns
colnames(Ib20) <- c("date_time", "unit", "temp_C")

#change time format
Ib20$date_time <- as.POSIXct(Ib20$date_time, format = "%d/%m/%y %I:%M:%S %p")

#plot
Ib20_C3 <- ggplot(data = subset(Ib20, month(date_time) == 10  &  day(date_time) >=15), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%m/%d", breaks = "1 day")
Ib20_C3

#load in ID # 02
Ib02 <- fread("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Leaf_packs/Campaign 3/02_E900000031CAED21_120122.csv")

#rename columns
colnames(Ib02) <- c("date_time", "unit", "temp_C")

#change time format
Ib02$date_time <- as.POSIXct(Ib02$date_time, format = "%d/%m/%y %I:%M:%S %p")

#plot
Ib02_C3 <- ggplot(data = subset(Ib02, month(date_time) == 10 &  day(date_time) >=15), aes(x = date_time, y = temp_C)) + geom_point() + scale_x_datetime(date_labels = "%m/%d", breaks = "1 day")
Ib02_C3

#What is the average temp
mean_temp02 <- mean(subset(Ib02, month(date_time) == 10 & day(date_time) >= 15)$temp_C, na.rm = TRUE)
mean_temp20 <- mean(subset(Ib20, month(date_time) == 10 & day(date_time) >= 15)$temp_C, na.rm = TRUE)

max_temp02 <- max(subset(Ib02, month(date_time) == 10 & day(date_time) >= 15)$temp_C, na.rm = TRUE)
memax_temp20 <- max(subset(Ib20, month(date_time) == 10 & day(date_time) >= 15)$temp_C, na.rm = TRUE)

#20 is on average a bit higher thus TA13, and 02 is more likely TA24


#####################################################


#### Load in the iButton master file for the leaf pack water iButtons ####

# ibutton information (file name, type, site, ibutton ID) can be found here: 
ibutton_master2 <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Leaf_packs/Leaf_pack_iButton_master_sheet.csv", header=T)

head(ibutton_master2)

#Add a new column for the file name
ibutton_master2 <- ibutton_master2 %>%
  rename(Device_address = File_name)   %>%
  mutate(
    iButton_ID = ifelse(nchar(iButton_ID) == 1, paste0("0", iButton_ID), iButton_ID),  #if single digit add a 0 in front of it
    file_ending = case_when(
      Campaign == 1 ~ "071222.csv",
      Campaign == 2 ~ "090222.csv",
      Campaign == 3 ~ "120122.csv",
      TRUE ~ ".csv"   # Change the ending depending on the campaign
    ),
    file_name = paste0(iButton_ID, "_", Device_address, "_", file_ending)
  ) %>%
  select(file_name, everything()) %>%
  rename(site = Site) #rename


str(ibutton_master2)  #58 obs. of  6 variables
head(ibutton_master2)

########################################################

##### Run loop to find the mean daily water temperature per site ####

#x <- subset(ibutton_master2, file_name=="05_ED00000031EC9621_071222.csv") #test with a subset first

all_data_water <- data.frame(DateTime = character(),
                       Unit = character(),
                       Temp_C = numeric(),
                       site = character(),
                       stringsAsFactors = FALSE)

for (i in 1:nrow(ibutton_master2)) {
  site <- ibutton_master2$site[i]
  file_name <- ibutton_master2$file_name[i]
  file_path <- paste("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/Leaf_packs/All/", file_name, sep = "")
  
  # Load the data from the file
  Dat <- fread(file_path, header = TRUE)
  colnames(Dat) <- c("DateTime", "Unit", "Temp_C")
  Dat$DateTime <- as.POSIXct(Dat$DateTime, format = "%d/%m/%y %I:%M:%S %p", tz = "Europe/Paris", origin = "1970-01-01")
  Dat$file_name <- file_name 
  
  # Add the 'site' column and fill it with the site value
  Dat$site <- site
  
  Dat <- Dat %>% drop_na()
  
  # Append the data to the all_data data frame
  all_data_water <- rbind(all_data_water, Dat)
}

str(all_data_water) #68746 obs. of  5 variables


levels(as.factor(all_data_water$site)) #all sites


# Now calculate the daily average temperature for each Site

daily_mean_water <- all_data_water %>%
  group_by(site, Date = as.Date(DateTime)) %>%
  summarise(DailyMeanWaterTemp_C = mean(Temp_C, na.rm = TRUE))

# Plot to check

# Convert the 'Date' column to POSIXct format for plotting
daily_mean_water$Date <- as.POSIXct(daily_mean_water$Date, format = "%Y-%m-%d")

TA02 <- ggplot(data = subset(daily_mean_water, site=="TA02"), aes(x = Date, y = DailyMeanWaterTemp_C)) + geom_point() + scale_x_datetime(date_labels = "%m/%d", breaks = "2 days")
TA02

TA42 <- ggplot(data = subset(daily_mean_water, site=="TA24"), aes(x = Date, y = DailyMeanWaterTemp_C)) + geom_point() + scale_x_datetime(date_labels = "%m/%d", breaks = "5 days")
TA42


##########################################################################

#### 4. Determine the relationship between daily mean air and water temperature for C1, by site. In order to extrapolate from the water temperature for C2 and C3 ####

#Merge and and water temperature by site and date

merged_data <- left_join(daily_ave_air_temp_C1, daily_mean_water, by = c("site", "Date"))  
   
merged_data_sub <- merged_data %>%
  ungroup() %>%
  filter(Date >= as.POSIXct("2022-05-30") & Date <= as.POSIXct("2022-06-19"))  #filter the date from the start of the campaign til the day before you collected the leaf packs

merged_data_sub <- as.data.frame(merged_data_sub) #make dataframe

str(merged_data_sub) #400 obs

#Drop certain rows where the iButton was not in the site yet
merged_data_sub <- merged_data_sub %>%
  filter(!(Date == as.Date("2022-05-30") & !(site %in% c("TA01", "TA02")))) %>%
  filter(!(Date == as.Date("2022-05-31") & !(site %in% c("TA01", "TA02", "TA03", "TA04", "TA05", "TA06")))) %>%
  filter(!(Date == as.Date("2022-06-01") & !(site %in% c("TA01", "TA02", "TA03", "TA04", "TA05", "TA06", "TA07", "TA08", "TA09", "TA10", "T011" )))) %>%
  filter(!(Date == as.Date("2022-06-02") & !(site %in% c("TA01", "TA02", "TA03", "TA04", "TA05", "TA06", "TA07", "TA08", "TA09", "TA10", "T011", "TA14", "TA15", "TA22", "TA24", "TA12" ))))

str(merged_data_sub) #353 obs

# Start by plotting the general relationship
all_temp <- ggplot(data=merged_data_sub, aes(x=DailyMeanWaterTemp_C, y=DailyMeanTemp_C) ) + geom_point() + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)  + stat_regline_equation(label.y = 29) + stat_cor(label.y = 27.5)
all_temp #R is 0.83 and significant so this is a strong relationship 

all_temp2 <- ggplot(data=merged_data_sub, aes(x=DailyMeanWaterTemp_C, y=DailyMeanTemp_C, colour=site) ) + geom_point() 
all_temp2

# Fit site-specific linear regressions and store coefficients in a new dataframe
site_coefficients <- merged_data_sub %>%
  group_by(site) %>%
  do(tidy(lm(DailyMeanTemp_C ~ DailyMeanWaterTemp_C, data = .))) %>%
  select(site, term, estimate) %>%
  spread(term, estimate) %>%
  rename(site = site, intercept = `(Intercept)`, slope = DailyMeanWaterTemp_C)

# Use site-specific coefficients to predict new air temperature
air_temp_predicted <- merged_data %>%
  left_join(site_coefficients, by = "site") %>%
  mutate(predicted_air_temp = intercept + slope * DailyMeanWaterTemp_C)

#plot

pred_tempC2 <- ggplot(data = air_temp_predicted[which(air_temp_predicted$Date < "2022-07-22" & air_temp_predicted$Date > "2022-07-18"),], aes(x = Date, y = predicted_air_temp, colour=site)) +   geom_point() 
pred_tempC2
 # The predictions do not make sense, some as high as 60 degrees!





