#### Script to extract air temperature data from the iButton dataloggers to get an average air temperature per GHG measurement to correct for in the gasflux calculation ####

## load necessary packages ##
library(data.table)
library(dplyr)
library(lubridate)
library(ggplot2)
library(survival)
library(scales)
library(tidyr)
library(purrr)

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

write.table(matTemp,file="C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/ibutton_air_temp.csv") 


### Read in the data and check for any outliers

#first you will need to read in the CSV's using fread because of the weird formatting of write.table

air.temp <- fread("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/ibutton_air_temp.csv") #For C1 

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
