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
# Note when missing data, the datalogger from the next closest site was used. NA in the iButton ID column denotes missing datalogger
#
#
# ibutton information (file name, type, site, ibutton ID) can be found here: 
ibutton_master <- read.csv("iButton_programming_Tarentaine_2022 - Campaign1.csv", header=T)

head(ibutton_master)

#Add a new column for the file name
ibutton_master <- ibutton_master %>%
  mutate(
    file_name = paste0(substr(site, nchar(site) - 1, nchar(site)), "_", Device_address,   ".csv")  ) %>%
  select(file_name, everything()) %>%    #reorder
  filter(habitat != "Water")   #Keep only the air iButtons

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


