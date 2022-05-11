#Aaron Lambert
#10/13/2021
#
#
############################################################################################
# This is code to get summary plots
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

######### Import Data ###############
# Historical Canadian EOS reconstructed run
# CAN_hist <- readRDS(file.path(dir.data,"/EOS Reconstructed Canadian.RDS"))
CAN_hist <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
# PSS historical for Jun 1st = 1, to Aug 31 = 92
# Years 1995 - 2020
# PSScum_hist <- readRDS(file.path(dir.data,"/clean_counts.rds"))

# Read in PSS non-cummulative with date in date format
PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
pf_hist <- as.data.frame(readRDS(file.path(dir.data,"preseason forcast.RDS")))
pf_hist$Year <- as.double(pf_hist$Year)
# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))

full.hist <- inner_join(CAN_hist, pf_hist)

# Plot Hstorical EOS Canadian reconstructed abundance
ggplot(data = full.hist,
       aes(x = Year, y = can.mean,
           ymin = can.mean - can.sd,
           ymax = can.mean+can.sd) )+
  geom_col(fill = "red")+
  geom_errorbar(col = "blue")+
  scale_x_discrete(limits = c(1995:2019))+
  theme(axis.text.x = element_text( angle = 90))+
  geom_point(aes(y = Mean))
  
  
