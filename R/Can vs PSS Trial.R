#Aaron Lambert
#10/13/2021
#
#
############################################################################################
# This is code to plot 6 panels in a loop
#
# This was developed to look at Yukon PSS counts up to certain days in the season 
#  vs Candian EOS counts

library(tidyverse)
library(lubridate)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(viridis)
library(ggpmisc)

# save working directory path as an object
wd <- getwd()

# Define the location of our data desired subfolder
dir.data <- file.path(wd, "data")

# Import data sets
# Cumulative counts by day from 1995 -- 2020
# clean_counts <- readRDS(file = file.path(dir.data,"clean_counts.RDS"))
# names(clean_counts)[names(clean_counts) == 'year'] <- 'Year'

# Get PSS passage data
# PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))
PSS_hist <- readRDS(file = file.path(dir.data,"pss adfg 27April22.RDS"))
# EOS Canadian reonstructed abundance
# CAN_hist <- readRDS(file.path(dir.data,"/EOS Reconstructed Canadian.RDS"))
# CAN_hist <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
CAN_hist <- readRDS(file.path(dir.data,"Can Origin Reconstructed 2Feb22.RDS"))
# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))

# Plot Historic Counts by Canadian EOS for different days of interest #####

# First by unadjusted PSS counts...
# PSS_cum_Hist<-PSS_hist %>%  group_by(Year) %>% mutate(CummCount= cumsum(count))

# This is for the restricted years for comparison against the precise GSI adj
PSS_cum_Hist<-PSS_hist %>% 
  # subset(Year >= 2005) %>% 
  group_by(Year) %>%
  mutate(CummCount= cumsum(count))

# Create an empty list to store plot outputs in
p <-list()

# For loop to create plots
for(i in seq(from = 152, to = 212, by = 10)){
# control variable to change day of interest  
  MYDAY <- i
  
  test<-PSS_cum_Hist[PSS_cum_Hist$Day==MYDAY,]
  
  test2 <- merge(x = test, y = CAN_hist, by = "Year")
  
p[[paste("Plot",i)]]<- ggplot(test2, aes(x = CummCount, y = log(can.mean))) +
    geom_point(aes(color =factor(Year)), show.legend = TRUE) +
    geom_smooth(method = "lm") +
    # xlab("PSS Sonar Count (Unadjusted)") +
    xlab("")+
    stat_poly_eq(small.r = TRUE)+
    # ylab("Log Total Reconstructed Canadian Run") +
    ylab("")+
    theme_classic(base_size = 12, base_family = "serif")+
    scale_y_continuous(trans = "log", labels = label_number(accuracy = 0.01))+
    scale_x_continuous(limits = c(min(test2$CummCount),max(test2$CummCount)))+
  labs(color = "Year")
    # geom_text(aes(label = Year))
  
}


# Arrange plots in grid

 unadj<-ggarrange(p$`Plot 162`,p$`Plot 172`,p$`Plot 182`,p$`Plot 192`,p$`Plot 202`,p$`Plot 212`,
          labels = c("Jun 11", " Jun 21 ", "Jul 1", "Jul 11", "Jul 21", "Jun 31", "31"),
          font.label = list(size = 12, color = "black", family = "serif"),
          vjust = 1,
          hjust = -3,
          common.legend = TRUE,
          legend = "right")

annotate_figure(unadj,
                left =  text_grob("Log Total Reconstructed Canadian Run",
                                  family = "serif",
                                  rot = 90, 
                                  vjust = 2,
                                  size = 11,
                                  face = "plain"),
                bottom = text_grob("PSS Sonar Cummulative Passage (Unadjusted)",
                                   family = "serif",
                                   vjust = -1,
                                   size = 11,
                                   face = "plain")
                )

# GSI AVG adjustment ##############################################

# Create new column to pop with adj value
PSS_hist$MeanGSIadj <- as.numeric("")

# Use loop to make the adj by mean GSI
for(i in 1:length(PSS_hist$Day)){
  # i = 29
  PSS_hist$MeanGSIadj[i] <- ifelse(PSS_hist$Day[i] < 173, PSS_hist$count[i] * GSI_mean[1,3],
                                   ifelse(PSS_hist$Day[i] >= 173 & PSS_hist$Day[i] < 181, PSS_hist$count[i] *GSI_mean[2,3],
                                          ifelse( PSS_hist$Day[i] >= 181, PSS_hist$count[i] * GSI_mean[3,3], PSS_hist$count[i] == "NA")))
}


#  get cumulative count for 2005-2020
PSS_cum_Hist_Mean_adj<-PSS_hist  %>% 
  # subset(Year >= 2005) %>%   
  group_by(Year) %>% 
  mutate(CummMeanAdjCount= cumsum(MeanGSIadj))

# Create an empty list to store plot outputs in
p1 <-list()

head(PSS_cum_Hist_Mean_adj)
# For loop to create plots
for(i in seq(from = 152, to = 212, by = 10)){
  # control variable to change day of interest  
  MYDAY <- i
  
  test<-PSS_cum_Hist_Mean_adj[PSS_cum_Hist_Mean_adj$Day==MYDAY,]
  
  test2 <- merge(x = test, y = CAN_hist, by = "Year")
  
  p1[[paste("Plot",i)]]<- ggplot(test2, aes(x = CummMeanAdjCount, y = log(can.mean))) +
    geom_point(aes(color = factor(Year)), show.legend = TRUE) +
    geom_smooth(method = "lm") +
    # xlab("PSS Sonar Count (Mean Adj)") +
    # ylab("Log Total Reconstructed Canadian Run") +
    xlab("")+
    ylab("")+
    theme_classic(base_size = 12, base_family = "serif")+
    scale_y_continuous(trans = "log", labels = label_number(accuracy = 0.01))+
    scale_x_continuous(limits = 
                         c(min(test2$CummMeanAdjCount),max(test2$CummMeanAdjCount)))+
    stat_poly_eq(small.r = TRUE)+
    theme(plot.margin = unit(c(0.5,.5,.5,.5),"cm"))
  
}



# Arrange plots in grid

meanAdj<-ggarrange(p1$`Plot 162`,p1$`Plot 172`,p1$`Plot 182`,p1$`Plot 192`,p1$`Plot 202`,p1$`Plot 212`,
          labels = c("Jun 11", " Jun 21 ", "Jul 1", "Jul 11", "Jul 21", "Jun 31", "31"),
          font.label = list(size = 12, color = "black", family = "serif"),
          vjust = 1,
          hjust = -3,
          common.legend = TRUE,
          legend = "right")



annotate_figure(meanAdj,
                left =  text_grob("Log Total Reconstructed Canadian Run",
                                  family = "serif",
                                  rot = 90, 
                                  vjust = 2,
                                  size = 11,
                                  face = "plain"),
                bottom = text_grob("PSS Sonar Cummulative Passage (Mean Adjusted)",
                                   family = "serif",
                                   vjust = -1,
                                   size = 11,
                                   face = "plain")
)




# Plots to explore difference in GSI with precise GSI * pss counts #######

# Create new column to pop with adj value
PSS_hist$PreciseGSIadj <- as.numeric("")

# Subset the data for years included in GSI data
PSS_sub_hist <- PSS_hist[PSS_hist$Year >= 2005,]

# Rename for easier use in loop
GSI<-GSI_by_year

head(GSI)

# Use loop to make the adj by precise GSI
for(i in 1:length(PSS_sub_hist$Day)){
  # i = 1811
  PSS_sub_hist$PreciseGSIadj[i] <- PSS_sub_hist$count[i] * GSI$propCan[GSI$startday<=PSS_sub_hist$Day[i]  &
                                                                 GSI$endday>=PSS_sub_hist$Day[i] &
                                                                 GSI$year == PSS_sub_hist$Year[i]]
}

head(PSS_sub_hist)

#  get cumulative count
PSS_cum_Hist_Precise_adj<-PSS_sub_hist %>%  group_by(Year) %>% mutate(CummPrecAdjCount= cumsum(PreciseGSIadj))

# Create an empty list to store plot outputs in
p2 <-list()


# For loop to create plots
for(i in seq(from = 152, to = 212, by = 10)){
  # control variable to change day of interest  
  MYDAY <- i
  
  test<-PSS_cum_Hist_Precise_adj[PSS_cum_Hist_Precise_adj$Day == MYDAY,]
  
  test2 <- merge(x = test, y = CAN_hist, by = "Year")
  
  p2[[paste("Plot",i)]]<- ggplot(test2, aes(x = CummPrecAdjCount, y = log(can.mean))) +
    geom_point(aes(color = factor(Year)), show.legend = TRUE) +
    geom_smooth(method = "lm") +
    # xlab("PSS Sonar Count (Precise Adj)") +
    xlab("")+
    stat_poly_eq(small.r = TRUE)+
    # ylab("Log Total Reconstructed Canadian Run") +
    ylab("")  +
    theme_classic(base_size = 12, base_family = "serif")+
    scale_y_continuous(trans = "log", labels = label_number(accuracy = 0.01))+
    scale_x_continuous(limits = c(min(test2$PreciseGSIadj[i]),max(test2$PreciseGSIadj[i])))
  
}


# Arrange plots in grid

precAdj<-ggarrange(p2$`Plot 162`,p2$`Plot 172`,p2$`Plot 182`,p2$`Plot 192`,p2$`Plot 202`,p2$`Plot 212`,
          labels = c("Jun 11", " Jun 21 ", "Jul 1", "Jul 11", "Jul 21", "Jun 31", "31"),
          font.label = list(size = 12, color = "black", family = "serif"),
          vjust = 1,
          hjust = -3,
          common.legend = T,
          legend = "right")

annotate_figure(precAdj,
                left =  text_grob("Log Total Reconstructed Canadian Run",
                                             family = "serif",
                                             rot = 90, 
                                             vjust = 2,
                                             size = 11,
                                             face = "plain"),
                bottom = text_grob("PSS Sonar Cummulative Passage (Precise Adjusted)",
                                              family = "serif",
                                              vjust = -1,
                                              size = 11,
                                              face = "plain")
)


# GSI AVG over days across years


meanGSI_vect <- vector(length = length(152:252))
names(meanGSI_vect) <- c(152:252)
sdGSI_vect <- vector(length = length(152:252))
names(sdGSI_vect) <- c(152:252)

counter <- 1
for (d in 152:252) {
  meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  counter <- counter+1
}

gsiDF <- data.frame("prop"= meanGSI_vect, "day" = 152:252)

i=1
for (i in 1:length(PSS_hist$Day)) {
  
  PSS_hist$gsiadjcount[i] <- PSS_hist$count[i]* gsiDF$prop[gsiDF$day == PSS_hist$Day[i]]
}


#  get cumulative count for 2005-2020
PSS_cum_Hist_Mean_adj<-PSS_hist  %>% 
  # subset(Year >= 2005) %>%   
  group_by(Year) %>% 
  mutate(CummMeanAdjCount= cumsum(gsiadjcount))

# Create an empty list to store plot outputs in
p1 <-list()

head(PSS_cum_Hist_Mean_adj)
# For loop to create plots
for(i in seq(from = 152, to = 212, by = 10)){
  # control variable to change day of interest  
  MYDAY <- i
  
  test<-PSS_cum_Hist_Mean_adj[PSS_cum_Hist_Mean_adj$Day==MYDAY,]
  
  test2 <- merge(x = test, y = CAN_hist, by = "Year")
  
  p1[[paste("Plot",i)]]<- ggplot(test2, aes(x = CummMeanAdjCount, y = log(can.mean))) +
    geom_point(aes(color = factor(Year)), show.legend = TRUE) +
    geom_smooth(method = "lm") +
    # xlab("PSS Sonar Count (Mean Adj)") +
    # ylab("Log Total Reconstructed Canadian Run") +
    xlab("")+
    ylab("")+
    theme_classic(base_size = 12, base_family = "serif")+
    scale_y_continuous(trans = "log", labels = label_number(accuracy = 0.01))+
    scale_x_continuous(limits = 
                         c(min(test2$CummMeanAdjCount),max(test2$CummMeanAdjCount)))+
    stat_poly_eq(small.r = TRUE)+
    theme(plot.margin = unit(c(0.5,.5,.5,.5),"cm"))
  
}



# Arrange plots in grid

meanAdj<-ggarrange(p1$`Plot 162`,p1$`Plot 172`,p1$`Plot 182`,p1$`Plot 192`,p1$`Plot 202`,p1$`Plot 212`,
                   labels = c("Jun 11", " Jun 21 ", "Jul 1", "Jul 11", "Jul 21", "Jun 31", "31"),
                   font.label = list(size = 12, color = "black", family = "serif"),
                   vjust = 1,
                   hjust = -3,
                   common.legend = TRUE,
                   legend = "right")



annotate_figure(meanAdj,
                left =  text_grob("Log Total Reconstructed Canadian Run",
                                  family = "serif",
                                  rot = 90, 
                                  vjust = 2,
                                  size = 11,
                                  face = "plain"),
                bottom = text_grob("PSS Sonar Cummulative Passage (Mean Adjusted)",
                                   family = "serif",
                                   vjust = -1,
                                   size = 11,
                                   face = "plain")
)






# GSI Plots by year #########################################33
un_GSI_by_year <- readRDS(file = file.path(dir.data,"un_adj_GSI_by_year.RDS"))

head(un_GSI_by_year)
hist(un_GSI_by_year$propCan)
shapiro.test(un_GSI_by_year$propCan)

GSI_by_year$stratumLength <- GSI_by_year$endday - GSI_by_year$startday

un_GSI_by_year %>%  ggplot(aes(x = stratum, y = propCan,  col = factor(year)))+
  geom_line()+
  geom_point()+
  xlab("Stratum")+
  ylab("Propotion of PSS Canadian Stock")+
  labs(color = "Year")

un_GSI_by_year %>% ggplot(aes(xmin = ((year)-.2),xmax = ((year)+.2),
                           ymin =startday, ymax = endday, 
                           col = factor(stratum), 
                           fill = propCan))+
  geom_rect(alpha=0.9,size = 1) +
  scale_x_continuous(breaks=2005:2020,limits=c(2004,2021)) +
  scale_y_continuous(breaks = waiver(),limits=c(150,250)) +
  xlab("Years") +
  ylab("Days") +
  # theme_solarized()+
  # theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c("Red","Green","Blue", "yellow"))+
  scale_fill_viridis(option = "A")+
  coord_flip()+
  labs(fill = "Proportion Canadian", color = "Stratum")
  


