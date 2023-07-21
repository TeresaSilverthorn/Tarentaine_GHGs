#### Data and statistical analysis and visualization for GHG data collected along the Albarine River network, France in 2021 ####

# Associated with the publication: 

# Load necessary packages 
library(data.table) #for fread
library(lubridate)
library(ggplot2)
library(scales) #for date_breaks
library(dplyr)
library(ggplot2)

library(ggpmisc)
library(ggpubr)
library(plyr)
library(gasfluxes)
library(tidyverse)
library(stringr)
library(car)
library(dpseg)

#########################################

# Set working directory #

setwd("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork 2022")

#########################################
#NOTES#

# Used the loaner Picarro (so maybe it is not 2 minutes ahead, need to check)

#########################################
### load raw data and metadata
### put the date-time in the same format for all files
##########################################


#Since the Picarro measures in GMT time, need to check the actual time zone of the measurement period and adjust


## CAMPAIGN 1 ##
#LOAD in the raw Picarro output data from campaign 1 (May 30 to June 3 2022)

Picarro_Spring2022<-do.call(rbind, lapply(list.files("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork 2022/Data/Picarro_dat/Campaign_1", pattern='dat', full.names=T, recursive=TRUE), fread ,header=T))
#to include sub-directories, change the recursive T

str(Picarro_Spring2022) # 78862 obs. of  22 variables

# create a column merging date and time
time<-as.POSIXct(paste(Picarro_Spring2022$DATE, Picarro_Spring2022$TIME), format="%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")
#Including the timezone is important as it prevents ggplot from displaying GMT time

Picarro_Spring2022<-cbind(Picarro_Spring2022,time)

str(Picarro_Spring2022) #78862 obs. of  23 variables

#Since the Picarro data is in GMT time, we need to add two hours to correspond to the actual time (GMT+2 time) for campaign 1 
Picarro_Spring2022$time<- as.POSIXlt(Picarro_Spring2022$time) +7200


## CAMPAIGN 2 ##
#LOAD in the raw Picarro output data from campaign 2 (July  18 to 22)
Picarro_Summer2022<-do.call(rbind, lapply(list.files("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork 2022/Data/Picarro_dat/Campaign_2", pattern='dat', full.names=T, recursive=TRUE), fread ,header=T))

str(Picarro_Summer2022) #93194 obs. of  22 variables

# create a column merging date and time
time<-as.POSIXct(paste(Picarro_Summer2022$DATE, Picarro_Summer2022$TIME), format="%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin") #Including the timezone is important as it prevents ggplot from displaying GMT time

Picarro_Summer2022<-cbind(Picarro_Summer2022,time)

str(Picarro_Summer2022) #93194 obs. of  23 variables

#Since the Picarro data is in GMT time, we need to add two hours to correspond to the actual time (GMT+2 time) for campaign 2 
Picarro_Summer2022$time<- as.POSIXlt(Picarro_Summer2022$time) +7200



## CAMPAIGN 3 ##
#LOAD in the raw Picarro output data from campaign 3 (October 10 to 14)
Picarro_Fall2022<-do.call(rbind, lapply(list.files("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork 2022/Data/Picarro_dat/Campaign_3", pattern='dat', full.names=T, recursive=TRUE), fread ,header=T))

str(Picarro_Fall2022) #53393 obs. of  22 variables

# create a column merging date and time
time<-as.POSIXct(paste(Picarro_Fall2022$DATE, Picarro_Fall2022$TIME), format="%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin") #Including the timezone is important as it prevents ggplot from displaying GMT time

Picarro_Fall2022<-cbind(Picarro_Fall2022,time)

str(Picarro_Fall2022) #93194 obs. of  23 variables

#Since the Picarro data is in GMT time, we need to add two hours to correspond to the actual time (GMT+2 time) for campaign 3 (Note the clocks didn't go back 1h until October 30, 2022) 
Picarro_Fall2022$time<- as.POSIXlt(Picarro_Fall2022$time) +7200

#############################################################################

#Plot the raw data for a rough check

spring_plot<- ggplot(data=Picarro_Spring2022[which(Picarro_Spring2022$time<"2022-05-30 16:30" & Picarro_Spring2022$time>"2022-05-30 5:00"),],aes(time, CO2_dry))+ geom_point() + scale_x_datetime(breaks=date_breaks("1 hour"), date_labels = "%I:%M")
spring_plot #looks good, decent peaks


summer_plot<- ggplot(data=Picarro_Summer2022[which(Picarro_Summer2022$time<"2022-07-18 16:30" & Picarro_Summer2022$time>"2022-07-18 5:00"),],aes(time, CO2_dry))+ geom_point() + scale_x_datetime(breaks=date_breaks("1 hour"), date_labels = "%I:%M")
summer_plot #Looks good, nice peaks


fall_plot<- ggplot(data=Picarro_Fall2022[which(Picarro_Fall2022$time<"2022-10-11 14:30" & Picarro_Fall2022$time>"2022-10-11 12:00"),],aes(time, CO2_dry))+ geom_point() + scale_x_datetime(breaks=date_breaks("1 hour"), date_labels = "%I:%M") + ylim(0,1000)
fall_plot #Looks good, nice peaks

#############################################################################

### Combine the campaigns together vertically ###

Picarro_2022 <- rbind(Picarro_Spring2022, Picarro_Summer2022,  Picarro_Fall2022)
str(Picarro_2022) #225,449 obs. of  23 variables

write.csv(Picarro_2022,"C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork 2022/Data/Picarro_raw_2022.csv")

#################################################

############## #Plots to check the time lag ##############
#Make a few plots and check the plotted start and end times vs the written start and end times of the measurements

#CAMPAIGN 1  # 

May30<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-05-30 14:15" & Picarro_2022$time<"2022-05-30 15:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(415, 550) 
May30 #TA01 #Picarro on average 1.83 minutes ahead

May30<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-05-30 18:00" & Picarro_2022$time<"2022-05-30 19:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(415, 550) 
May30 #TA02 #Picarro on average 2.42 minutes ahead

May31<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-05-31 10:00" & Picarro_2022$time<"2022-05-31 11:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,650)
May31 #TA04 #Picarro on average 1.83 minutes ahead

May31<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-05-31 12:00" & Picarro_2022$time<"2022-05-31 13:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(415, 440) 
May31 #TA05 #Picarro on average 2.5 minutes ahead

May31<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-05-31 14:50" & Picarro_2022$time<"2022-05-31 15:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(415, 440) 
May31 #TA06 #Picarro on average 1.6 minutes ahead

May31<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-05-31 17:00" & Picarro_2022$time<"2022-05-31 18:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(410, 440) 
May31 #TA03  #Picarro on average 2.25 minutes ahead

Jun1<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-01 8:00" & Picarro_2022$time<"2022-06-01 11:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,800)
Jun1 #TA09 #Picarro on average 2 minutes ahead

Jun1<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-01 11:40" & Picarro_2022$time<"2022-06-01 12:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,500)
Jun1 #TA10  #Picarro on average 1.58 minutes ahead

Jun1<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-01 13:00" & Picarro_2022$time<"2022-06-01 13:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,800)
Jun1 #TA11 #Picarro on average 1.83 minutes ahead

Jun1<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-01 14:40" & Picarro_2022$time<"2022-06-01 16:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,500)
Jun1 #TA07 #Picarro on average 1.83 minutes ahead

Jun1<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-01 16:40" & Picarro_2022$time<"2022-06-01 17:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,600)
Jun1 #TA08 #Picarro on average 2 minutes ahead

Jun2<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-02 9:40" & Picarro_2022$time<"2022-06-02 10:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(400,1000)
Jun2 #TA14  #Picarro on average 1.75 minutes ahead

Jun2<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-02 11:40" & Picarro_2022$time<"2022-06-02 12:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(400,1000)
Jun2 #TA15 #Picarro on average 1.58 minutes ahead

Jun2<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-02 14:00" & Picarro_2022$time<"2022-06-02 15:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(400,550)
Jun2 #TA22   #Picarro on average 1.75 minutes ahead

Jun2<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-02 15:00" & Picarro_2022$time<"2022-06-02 16:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(400,550)
Jun2 #TA24  #Picarro on average 2 minutes ahead

Jun2<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-02 17:00" & Picarro_2022$time<"2022-06-02 18:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")+ ylim(400,1000)
Jun2 #TA12 #Picarro on average 1.7 minutes ahead

Jun3<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-03 10:00" & Picarro_2022$time<"2022-06-03 11:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,1500)
Jun3 #TA20 #Picarro on average 2 minutes ahead

Jun3<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-03 12:00" & Picarro_2022$time<"2022-06-03 13:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,1000)
Jun3 #TA21 #Picarro on average 2 minutes ahead

Jun3<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-03 13:50" & Picarro_2022$time<"2022-06-03 14:50"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,600)
Jun3 #TA13  #Picarro on average 1.92 minutes ahead

Jun3<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-03 15:00" & Picarro_2022$time<"2022-06-03 18:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,650)
Jun3  #TA17  #Picarro on average 1.25 minutes ahead




##########################################################
#load the ancillary data
##########################################################

#Import ancillary data 

#Update as necessary to most recent file in the Ancillary data folder in case you make any corrections to the Google Drive sheet

ancil_dat <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork 2022/Data/Ancillary data/Data_entry_Tarentaine_2022_2023-07-21.csv", header=T)

str(ancil_dat) # 331 obs. of  36 variables


#add an ID column  
ancil_dat <- ancil_dat %>%
  mutate(ID_unique = paste(date, site, transect, sep = "_")) %>%
  select(ID_unique, everything()) #make this ID_unique column come first

# Add date to the start and end time columns and format as as.POSITX
ancil_dat <- ancil_dat %>%
  mutate(datetime_start = as.POSIXct(paste(date, start_time), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris"),
         datetime_end = as.POSIXct(paste(date, end_time), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris"),
         datetime_start = if_else(is.na(datetime_start), as.POSIXct(paste(date, start_time), format = "%Y-%m-%d %H:%M", tz = "Europe/Paris"), datetime_start),
         datetime_end = if_else(is.na(datetime_end), as.POSIXct(paste(date, end_time), format = "%Y-%m-%d %H:%M", tz = "Europe/Paris"), datetime_end)) %>%
  select(ID_unique, campaign, date, site, flow_state, datetime_start, datetime_end, everything())  %>% #reorder
mutate(date = as.POSIXct(date, format = "%Y-%m-%d", tz = "Europe/Paris") ) #make date as.positx too


