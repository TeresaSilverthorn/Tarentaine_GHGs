#### Water level Tarentaine ####

#load necessary packages
library(tidyverse)
#########################


## Data from Vincent regarding the water level at each dam. 

#We can have an indirect idea of the water spilled from the dams using the water level monitoring. You’ll find attached an xls files with the water level in the 2 reservoirs. When it exceeds the elevation of the spillway (891.00 for Tarentaine and 887.50 for Eau Verte), water is flowing downstream (in addition to the environmental flow). It could help you to see the environmental conditions during and the days/weeks before your sampling campaigns to explain/discuss your results.


setwd("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Tarentaine_GHGs/Figures")


#### 1. Load water level data data ####

water_level <- read.csv ("C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Fieldwork_2022/Data/Water_level/Water_level.csv", header=T) 

#format date

water_level$X <-as.POSIXct(water_level$X, format="%Y-%m-%d %H:%M", tz = "Europe/Paris") #Including the timezone is important as it prevents ggplot from displaying GMT time

#Pivot for easier plotting

# Assuming your original data frame is called 'water_level_data'
water_level2 <- water_level %>%
  select(X, WL.Tarentaine, WL.Eau.Verte) %>%
  pivot_longer(cols = starts_with("WL"), names_to = "Reservoir", values_to = "WaterLevel") %>%
  mutate(Reservoir = gsub("^WL\\.", "", Reservoir),
         Reservoir = gsub("\\.", " ", Reservoir))

# View the first few rows of the new data frame
str(water_level2)




tiff("water_level_ts", units="in", width=12, height=6, res=300)

ts <- ggplot(water_level2, aes(x=X, y=WaterLevel)) +  theme_bw()  +  
  geom_rect(aes(xmin = as.POSIXct("2022-05-30"), xmax = as.POSIXct("2022-06-03"), ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.2) +
  geom_rect(aes(xmin = as.POSIXct("2022-07-18"), xmax = as.POSIXct("2022-07-22"), ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.2) +
  geom_rect(aes(xmin = as.POSIXct("2022-10-10"), xmax = as.POSIXct("2022-10-14"), ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.2) +
  geom_line(aes(group=Reservoir, colour=Reservoir) ) + 
  geom_hline(yintercept = 891, linetype = "dashed", colour='#5bd1c4') +   
  geom_hline(yintercept = 887.50 , linetype = "dashed", colour="#d15ba3") +
  theme(legend.position="top", axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),  panel.border = element_blank(),  axis.line = element_line(colour = "black"), text = element_text(size = 12), legend.title=element_blank(), axis.text = element_text(size = 12, colour="black"))  + 
  ylab(expression("Water level (masl)")) + scale_x_datetime(labels = scales::date_format("%b-%d"), breaks = scales::breaks_pretty(n = 12)) + scale_color_manual(values = c('#d15ba3',  "#5bd1c4"))
 ts

dev.off()

