#### Data and statistical analysis and visualization for GHG data collected along the Albarine River network, France in 2021 ####

# Associated with the publication: 

# Load necessary packages 
library(data.table) #for fread
library(lubridate)
library(ggplot2)
library(scales) #for date_breaks
library(dplyr)
library(ggplot2)
library(gasfluxes)

#library(ggpmisc)
#library(ggpubr)
#library(plyr)
#library(tidyverse)
#library(stringr)
#library(car)

#########################################

# Set working directory #

setwd("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022")

#########################################
#NOTES#

# Used the loaner Picarro (so maybe it is not 2 minutes ahead, need to check)

#########################################
### load raw Picarro data ####
### put the date-time in the same format for all files
##########################################


#Since the Picarro measures in GMT time, need to check the actual time zone of the measurement period and adjust


## CAMPAIGN 1 ##
#LOAD in the raw Picarro output data from campaign 1 (May 30 to June 3 2022)

Picarro_Spring2022<-do.call(rbind, lapply(list.files("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Picarro_dat/Campaign_1", pattern='dat', full.names=T, recursive=TRUE), fread ,header=T))
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
Picarro_Summer2022<-do.call(rbind, lapply(list.files("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Picarro_dat/Campaign_2", pattern='dat', full.names=T, recursive=TRUE), fread ,header=T))

str(Picarro_Summer2022) #93194 obs. of  22 variables

# create a column merging date and time
time<-as.POSIXct(paste(Picarro_Summer2022$DATE, Picarro_Summer2022$TIME), format="%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin") #Including the timezone is important as it prevents ggplot from displaying GMT time

Picarro_Summer2022<-cbind(Picarro_Summer2022,time)

str(Picarro_Summer2022) #93194 obs. of  23 variables

#Since the Picarro data is in GMT time, we need to add two hours to correspond to the actual time (GMT+2 time) for campaign 2 
Picarro_Summer2022$time<- as.POSIXlt(Picarro_Summer2022$time) +7200



## CAMPAIGN 3 ##
#LOAD in the raw Picarro output data from campaign 3 (October 10 to 14)
Picarro_Fall2022<-do.call(rbind, lapply(list.files("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Picarro_dat/Campaign_3", pattern='dat', full.names=T, recursive=TRUE), fread ,header=T))

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

write.csv(Picarro_2022,"C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Picarro_raw_2022.csv")

#################################################

############## #Plots to check the time lag ##############

#Make a few plots and check the plotted start and end times vs the written start and end times of the measurements

#### CAMPAIGN 1  #### 

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

Jun1<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-01 11:55" & Picarro_2022$time<"2022-06-01 12:10"),],aes(time, CH4_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") 
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

Jun2<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-06-02 15:00" & Picarro_2022$time<"2022-06-02 16:30"),],aes(time, CH4_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M")
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

#### CAMPAIGN 2  ##### 

Jul18<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-18 11:15" & Picarro_2022$time<"2022-07-18 13:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels = "%I:%M") + ylim(400,1100)
Jul18  #TA11 #Picarro on average 2.5 minutes ahead 

Jul18<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-18 14:15" & Picarro_2022$time<"2022-07-18 15:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,700)
Jul18  #TA10  #Picarro on average 2.25 minutes ahead 

Jul18<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-18 15:20" & Picarro_2022$time<"2022-07-18 16:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,1000)
Jul18  #TA09   #Picarro on average 2.3 minutes ahead 

Jul18<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-18 18:00" & Picarro_2022$time<"2022-07-18 19:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") 
Jul18 #TA08  #Picarro on average 2.42 minutes ahead 

Jul19<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-19 8:40" & Picarro_2022$time<"2022-07-19 10:10"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M")  + ylim(450,850)
Jul19 #TA01  #Picarro on average 2.2 minutes ahead 

Jul19<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-19 9:30" & Picarro_2022$time<"2022-07-19 10:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M")  + ylim(450,850)
Jul19 #TA01  #Picarro on average 2.2 minutes ahead 

Jul19<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-19 14:30" & Picarro_2022$time<"2022-07-19 15:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,500)
Jul19 #TA02  #Picarro on average 2.67 minutes ahead 

Jul19<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-19 17:30" & Picarro_2022$time<"2022-07-19 18:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(250,550)
Jul19 #TA03  #Picarro on average 2.25 minutes ahead 

Jul20<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-20 10:40" & Picarro_2022$time<"2022-07-20 11:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,650)
Jul20  #TA04  #Picarro on average 2.4 minutes ahead 

Jul20<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-20 11:30" & Picarro_2022$time<"2022-07-20 12:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(435,475)
Jul20  #TA05 #Picarro on average 2.2 minutes ahead 

Jul20<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-20 14:00" & Picarro_2022$time<"2022-07-20 15:20"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,480)
Jul20  #TA06 #Picarro on average 2.2 minutes ahead 

Jul20<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-20 16:30" & Picarro_2022$time<"2022-07-20 18:20"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(370,500)
Jul20  #TA07 #Picarro on average 2.6 minutes ahead 

Jul21<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-21 9:30" & Picarro_2022$time<"2022-07-21 10:50"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,700)
Jul21  #TA15 #Picarro on average 2.2 minutes ahead 

Jul21<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-21 11:00" & Picarro_2022$time<"2022-07-21 12:10"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,610)
Jul21  #TA14  #Picarro on average 2.75 minutes ahead 

Jul21<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-21 13:20" & Picarro_2022$time<"2022-07-21 14:20"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,710)
Jul21 #TA13  #Picarro on average 2.67 minutes ahead 

Jul21<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-21 15:20" & Picarro_2022$time<"2022-07-21 16:20"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,810)
Jul21  #TA20  #Picarro on average 2.25 minutes ahead 

Jul21<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-21 16:40" & Picarro_2022$time<"2022-07-21 17:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,800)
Jul21  #TA21   #Picarro on average 2.33 minutes ahead 

Jul22<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-22 10:20" & Picarro_2022$time<"2022-07-22 11:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,720)
Jul22  #TA17  #Picarro on average 2.5 minutes ahead 

Jul22<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-22 12:00" & Picarro_2022$time<"2022-07-22 13:15"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,550)
Jul22  #TA22  #Picarro on average 2.4 minutes ahead 

Jul22<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-07-22 14:40" & Picarro_2022$time<"2022-07-22 15:50"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,650)
Jul22 #TA24  #Picarro on average 2.67 minutes ahead 

##########################################################

#### CAMPAIGN 3  ##### 

Oct10<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-10 14:00" & Picarro_2022$time<"2022-10-10 15:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,650)
Oct10  # TA01    #Picarro on average 2.1 minutes ahead  #some weird flat lines, so 4 and 5 could be bad and need to be removed...

Oct10<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-10 16:20" & Picarro_2022$time<"2022-10-10 17:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,850)
Oct10  #TA02R  #Picarro on average 2.67 minutes ahead 

Oct11<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-11 9:00" & Picarro_2022$time<"2022-10-11 10:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(500,850)
Oct11  #TA03   #Picarro on average 2.3 minutes ahead 

Oct11<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-11 12:00" & Picarro_2022$time<"2022-10-11 13:50"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,600)
Oct11  #TA07 #Picarro on average 2.83 minutes ahead 

Oct11<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-11 15:00" & Picarro_2022$time<"2022-10-11 16:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(400,600)
Oct11 #TA15  #Picarro on average 2.83 minutes ahead 

Oct12<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-12 8:40" & Picarro_2022$time<"2022-10-12 10:20"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(450,600)
Oct12 #TA04 #Picarro on average 2.7 minutes ahead 

Oct12<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-12 10:40" & Picarro_2022$time<"2022-10-12 11:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(430,500)
Oct12 #TA05   #Picarro on average 2.9 minutes ahead 

Oct12<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-12 13:20" & Picarro_2022$time<"2022-10-12 14:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,700)
Oct12  #TA09  #Picarro on average 2.6 minutes ahead 

Oct12<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-12 14:30" & Picarro_2022$time<"2022-10-12 15:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,500)
Oct12  #TA10  #Picarro on average 2.7 minutes ahead 

Oct12<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-12 15:30" & Picarro_2022$time<"2022-10-12 16:10"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,1000)
Oct12  #TA11 #Picarro on average 2.67 minutes ahead 

Oct12<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-12 16:30" & Picarro_2022$time<"2022-10-12 17:20"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,600)
Oct12  #TA08  #Picarro on average 2.3 minutes ahead 

Oct12<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-12 17:50" & Picarro_2022$time<"2022-10-12 19:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,800)
Oct12 #TA14 #Picarro on average 2.2 minutes ahead 

Oct13<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-13 8:20" & Picarro_2022$time<"2022-10-13 9:30"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,1100)
Oct13 #TA21   #Picarro on average 2.25 minutes ahead

Oct13<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-13 9:40" & Picarro_2022$time<"2022-10-13 11:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") #+ ylim(420,1100)
Oct13  #TA20   #Picarro on average 2.3 minutes ahead

Oct13<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-13 11:15" & Picarro_2022$time<"2022-10-13 12:10"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,800)
Oct13 #TA13  #Picarro on average 2.4 minutes ahead

Oct13<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-13 14:35" & Picarro_2022$time<"2022-10-13 15:40"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,500)
Oct13  #TA06  #Picarro on average 2 minutes ahead

Oct13<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-13 17:00" & Picarro_2022$time<"2022-10-13 18:10"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,900)
Oct13 #TA17  #Picarro on average 2.2 minutes ahead

Oct14<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-14 10:00" & Picarro_2022$time<"2022-10-14 11:00"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,550)
Oct14 #TA22   #Picarro on average 2.5 minutes ahead

Oct14<- ggplot(data=Picarro_2022[which(Picarro_2022$time>"2022-10-14 11:20" & Picarro_2022$time<"2022-10-14 12:10"),],aes(time, CO2_dry))+ geom_point(size=0.5) + geom_line(size=0.1, alpha=0.5, colour="red") + scale_x_datetime(breaks=date_breaks("2 min"), date_labels= "%I:%M") + ylim(420,550)
Oct14 #TA24  #Picarro on average 2.2 minutes ahead

##########################################################
#### load the ancillary data ####
##########################################################

#Import ancillary data 

#Update as necessary to most recent file in the Ancillary data folder in case you make any corrections to the Google Drive sheet

ancil_dat <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Ancillary data/Data_entry_Tarentaine_2022_2023-07-26.csv", header=T)

str(ancil_dat) # 331 obs. of  37 variables

#Make site a factor
ancil_dat$site <- as.factor(ancil_dat$site)
ancil_dat$flow_state <- as.factor(ancil_dat$flow_state)

ancil_dat$water_temp_C <- as.numeric(ancil_dat$water_temp_C, na.rm=T)
ancil_dat$canopy_cover_. <-  as.numeric(ancil_dat$canopy_cover_., na.rm=T) 

#add an ID column  
ancil_dat <- ancil_dat %>%
  mutate(ID_unique = paste(date, site, transect, sep = "_")) %>%
  select(ID_unique, everything()) #make this ID_unique column come first

# Add date to the start and end time columns and format as as.POSITX

ancil_dat <- ancil_dat %>%
  mutate(
    start_time = ifelse(nchar(start_time) == 5, paste0(start_time, ":00"), start_time),
    end_time = ifelse(nchar(end_time) == 5, paste0(end_time, ":00"), end_time),
    datetime_start = as.POSIXct(paste(date, start_time), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris"),
    datetime_end = as.POSIXct(paste(date, end_time), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Paris")
  ) %>%
  select(ID_unique, campaign, date, site, flow_state, datetime_start, datetime_end, everything()) %>% #reorder
  mutate(date = as.POSIXct(date, format = "%Y-%m-%d", tz = "Europe/Paris")) #make date as.POSIXct too

# Add the nearest (in time) iButton air temperature from the iButton_air_temp.R script

air_temp <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/iButtons/ibutton_air_temp.csv")

#Subset the relevant rows and merge with ancil dat
air_temp  <- air_temp %>% select(-X)

#Merge with ancil_dat
ancil_dat <- merge (ancil_dat, air_temp , by="ID_unique")

# Add the air pressure data

pressure <- read.csv("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Pressure/Air_pressure_Tarentaine.csv")

#Subset just the relevatn columns

pressure <- pressure %>% 
  select(campaign, site, calculated.atm)

#merge with ancil data by site and campaign
ancil_dat <- merge(ancil_dat, pressure, by = c("site", "campaign"))

########################################################################

########## Adjust start and end times based on Picarro time drift #########

#On average LGR is 2 minutes ahead, I double checked each sampling date manually, so we have a specific time difference between the datasheet and Picarro time (minutes) in the columns "Picarro_start_offset_mins" and "Picarro_end_offset_mins"

ancil_dat <- ancil_dat %>%
  mutate(
    CO2datetime_start = datetime_start + (CO2_start_offset_mins*60),
    CO2datetime_end = datetime_end + (CO2_end_offset_mins*60) ,
    CH4datetime_start = datetime_start + (CH4_start_offset_mins*60),
    CH4datetime_end = datetime_end + (CH4_end_offset_mins*60) ) %>%
  select(ID_unique, campaign, date, site, flow_state, start_time, end_time, CO2datetime_start, CO2datetime_end, CO2_start_offset_mins, CO2_end_offset_mins,  CH4datetime_start, CH4datetime_end,  CH4_start_offset_mins, CH4_end_offset_mins,  everything()) 


#######################################################

#### Clip the Picarro data to the start and end times for CO2 #####

#If you clean CO2 and CH4 differently 

#To start, subset campaign 1/2: 
#ancil_dat_C3 <- ancil_dat %>%
 # filter(campaign %in% c("3"))

ID <- ancil_dat$ID_unique
startT<-ancil_dat$CO2datetime_start #start times
endT<-ancil_dat$CO2datetime_end  # end times 

for(i in 1:length(startT)){
  st<-startT[i]
  se<-endT[i]
  id<-ID[i]
  data<-Picarro_2022[Picarro_2022$time >= st & Picarro_2022$time <= se,]
  data$ID<-id
  
  if(i==1){
    assign(paste("data",i, sep="_"),data)
  } else {
    assign(paste("data",i, sep="_"),data)
    assign(paste("data",i, sep="_"),rbind(get(paste("data",i, sep="_")),get(paste("data",i-1, sep="_"))))
  }
}

Picarro_dat_CO2<-get(paste("data",length(startT),sep="_"))


str(Picarro_dat_CO2) #83828 obs. of  24 variables 

rm(list = ls()[grep("^data_", ls())]) #clear all of the clipped datasets

#######################################################

#### Clip the Picarro data to the start and end times for CH4 #####

#For now, subset campaign 1/2: 
#ancil_dat_C3 <- ancil_dat %>%
 # filter(campaign %in% c("3"))

ID <- ancil_dat$ID_unique
startT<-ancil_dat$CH4datetime_start #start times
endT<-ancil_dat$CH4datetime_end  # end times 

for(i in 1:length(startT)){
  st<-startT[i]
  se<-endT[i]
  id<-ID[i]
  data<-Picarro_2022[Picarro_2022$time >= st & Picarro_2022$time <= se,]
  data$ID<-id
  
  if(i==1){
    assign(paste("data",i, sep="_"),data)
  } else {
    assign(paste("data",i, sep="_"),data)
    assign(paste("data",i, sep="_"),rbind(get(paste("data",i, sep="_")),get(paste("data",i-1, sep="_"))))
  }
}

Picarro_dat_CH4<-get(paste("data",length(startT),sep="_"))


str(Picarro_dat_CH4) #84674 obs. of  24 variables 

rm(list = ls()[grep("^data_", ls())]) #clear all of the clipped datasets

#######################################################

#### Other data cleaning ####

#For CH4 there was a weird peak, so you will need to remove 2022-06-01_TA10_6

Picarro_dat_CH4 <- Picarro_dat_CH4 %>%
  filter(ID != "2022-06-01_TA10_6")

# CO2 bad flux 2022-10-12_TA14_2 and 2022-10-12_TA14_5 and 2022-10-14_TA24_5  and 2022-10-10_TA01_4 and 2022-10-12_TA14_5 and 2020-10-14_TA24_5 because the Picarro wasn't working/froze

Picarro_dat_CO2 <- Picarro_dat_CO2 %>%
  filter(!ID %in% c("2022-10-12_TA14_2", "2022-10-10_TA01_4", "2022-10-12_TA14_5", "2020-10-14_TA24_5" ))


#In C3 there are some negative CO2 values, 47 values. Need to remove them for gas fluxes. 

Picarro_dat_CO2 <- Picarro_dat_CO2 %>%
  filter(CO2_dry >= 0)

#######################################################

####  In order to set the initial time of each measurement to 0:  ####
####  start by adding a column for epoch time (expressed as seconds since Jan 1, 1970) ####
Picarro_dat_CO2$epoch_time <- as.integer(as.POSIXct(Picarro_dat_CO2$time), tz="Europe/Paris")
Picarro_dat_CH4$epoch_time <- as.integer(as.POSIXct(Picarro_dat_CH4$time), tz="Europe/Paris")

str(Picarro_dat_CO2)  #32264 obs. of  25 variables
str(Picarro_dat_CH4)  #32264 obs. of  25 variables 

#there are duplicated time rows in the data, delete them (can cause problems later)
dupsCO2 <- Picarro_dat_CO2[duplicated(epoch_time)]
dupsCH4 <- Picarro_dat_CH4[duplicated(epoch_time)]

str(dupsCO2)  # 248 obs. of  25 variables
str(dupsCH4)  # 457 obs. of 25

Picarro_dat_CO2 <- Picarro_dat_CO2 %>% 
  # Base the removal on the "epoch_time" column
  distinct(epoch_time, .keep_all = TRUE)

Picarro_dat_CH4 <- Picarro_dat_CH4 %>% 
   #Base the removal on the "epoch_time" column
  distinct(epoch_time, .keep_all = TRUE)

str(Picarro_dat_CO2)  # 32194 obs. of  25 variables
str(Picarro_dat_CH4)  # 32194 obs. of  25 variables

#then set  the initial time of each measure to 0h  (use Naiara's function to rest the min time to each time)
rescale <- function(x) (x-min(x))

#apply this function to all epoch_time (seconds) of each measure, 
#and divide by 3600 (hours)
Picarro_dat_CO2 <- setDT(Picarro_dat_CO2)[,c("flux_time"):=.(rescale(epoch_time/3600)),by=.(ID)]
Picarro_dat_CH4 <- setDT(Picarro_dat_CH4)[,c("flux_time"):=.(rescale(epoch_time/3600)),by=.(ID)]

## keep only the desired columns from Picarro data
Picarro_dat_CO2 <- subset(Picarro_dat_CO2, select = c( "CavityTemp","CO2_dry", "ID", "flux_time", "EPOCH_TIME"))
str(Picarro_dat_CO2) # 32194 obs. of  6 variables

Picarro_dat_CH4  <- subset(Picarro_dat_CH4, select = c( "CavityTemp", "CH4_dry",  "ID", "flux_time", "EPOCH_TIME"))
str(Picarro_dat_CH4)

#rename ID to ID unique
names(Picarro_dat_CO2)[names(Picarro_dat_CO2) == "ID"] <- "ID_unique"
names(Picarro_dat_CH4)[names(Picarro_dat_CH4) == "ID"] <- "ID_unique"

CO2_dat <- merge (Picarro_dat_CO2, ancil_dat , by="ID_unique", allow.cartesian=TRUE)
CH4_dat <- merge (Picarro_dat_CH4, ancil_dat , by="ID_unique", allow.cartesian=TRUE)

#flux time needs to be ordered within each ID (use ID unique for all of the data)
CO2_dat <- data.table(CO2_dat, key = c("ID_unique", "flux_time")) 
CH4_dat <- data.table(CH4_dat, key = c("ID_unique", "flux_time"))

str(CO2_dat) #82801 obs. of  52 variables
str(CH4_dat) #84075 obs. of  52 variables


#################################################################################

######  Convert Picarro_dat concentration from ppm to mg-CO2-C/L and ug CH4-C/L using the ideal gas law (PV=nRT) for input to gasfluxes package.  #############  Note that ppm = uL/L

#Need to update with this with the ibutton temperature and maybe get better air pressure values

#put CO2 and CH4 concentration (now in ppm) in mg/L
# mg/L = ((ppm  * molecular mass *1 atm )/1000) / (0.082 * 293K )
# ug/L = (ppm  * molecular mass *1 atm ) / (0.082 * 293K ) #devide whether in ug or mg

#Note here that we use the calculated atmospheric pressure and iButton datalogger temperature closest to the measurement time

CO2_dat$CO2_mg_L <- ((CO2_dat$CO2_dry  * 12.011 * CO2_dat$calculated.atm)/1000) / (0.08206 *(CO2_dat$temp_C + 273.15))
CH4_dat$CH4_mg_L <- ((CH4_dat$CH4_dry  * 12.011 * CH4_dat$calculated.atm)/1000) / (0.08206*(CH4_dat$temp_C + 273.15))


#Check the units
#V = L
#A = m2
# flux time = h
# concentration of CO2 / CH4 = mg/L

#[f0] = mg/m^2/h

mean(CO2_dat$chamber_volume) # 3.585717.....  L
mean(CO2_dat$chamber_area) # 0.04523893..... m2
median(CO2_dat$CO2_mg_L) # 0.224736........ mg/L   
median(CH4_dat$CH4_mg_L)  #0.00102........... mg/L
mean(CO2_dat$flux_time) # 0.04 h = 2.4 minutes
mean(CO2_dat$calculated.atm)  #  0.8956438 atm
mean(CO2_dat$temp_C)  #20.17568 C

###########################################################################

#### run gasfluxes package ####

#### for CO2 ####

setwd("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Tarentaine_GHGs/Flux_figures/gasfluxes/CO2")

#Run the package to calculate the gas flux rate
CO2.results <- gasfluxes(CO2_dat, .id = "ID_unique", .V = "chamber_volume", .A = "chamber_area",.times = "flux_time", .C = "CO2_mg_L",method = c("linear"), plot = F) #can turn plot to FALSE if the number of plots was getting out of hand

CO2.results #linear.f0 units are mg-CO2-C/m2/h
str(CO2.results) #  324 obs. of  10 variables


#Find out which ones are missing 
ancil_dat$ID_unique[!(ancil_dat$ID_unique %in% CO2.results$ID_unique)] # "2022-10-10_TA01_4" "2022-10-12_TA14_2" "2022-10-12_TA14_5" "2022-10-14_TA24_5"

#What we removed previously: "2022-10-10_TA01_4", "2022-10-12_TA14_2", "2022-10-12_TA14_5", "2022-10-14_TA24_5"

# Merge the flux data with the ancillary data

CO2.results<- subset(CO2.results, select = c( "ID_unique", "linear.f0"))

names(CO2.results)[names(CO2.results) == "linear.f0"] <- "CO2_C_mg_m2_h" #rename

CO2_fluxes <- full_join(CO2.results, ancil_dat, by= "ID_unique")

str(CO2_fluxes) #327 obs. of  48 variables (4 NA flux values)

##############################################################

#### for CH4 ####

 #  326 obs. of  10 variables

#Find out if any / which ones are missing 
ancil_dat$ID_unique[!(ancil_dat$ID_unique %in% CH4.results$ID_unique)]  #"2022-06-01_TA10_6"


# Merge the flux data with the ancillary data

CH4.results<- subset(CH4.results, select = c( "ID_unique", "linear.f0"))

names(CH4.results)[names(CH4.results) == "linear.f0"] <- "CH4_C_mg_m2_h" #rename

CH4_fluxes <- full_join(CH4.results, ancil_dat, by= "ID_unique")

str(CH4_fluxes) #327 obs. of  48 variables (1 NA value) 

#################################################################

#Merge the CH4 and CO2 outputs

#Merge the CO2 and CH4 results by ID
CO2.CH4.results <- merge (CO2.results, CH4.results, by= "ID_unique")
str(CO2.CH4.results) #322 obs. of  49 variables

#Merge with the ancillary data
CO2.CH4.fluxes <- merge (CO2.CH4.results, ancil_dat, by= "ID_unique")
str(CO2.CH4.fluxes) #252 obs. of  15 variables


# Save as csv

write.csv (CO2.CH4.fluxes, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Tarentaine_GHGs/CO2.CH4.fluxes.csv")

