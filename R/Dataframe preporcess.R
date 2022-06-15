#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Purpose: To preprocess Pilot Station Sonar for the model
#
#=================================================================================
#NOTES:
# This script will be where data is processed for entry into Yukon Inseason Forecast.R script
# 
#
# PSS historical counts
# Data imported from  
# https://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareayukon.salmon_escapement
# Import as xlsx file (excel)

library(tidyverse)
library(lubridate)
library(readxl)

# set to working directory
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

# Objects used to save/load data, outputs, or stan/R scripts
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

# Read in the file as xlsx
PSS <- read_xlsx(file.path(dir.data,"Yukon Escapement Daily 9June22.xlsx"),skip = 3)

# Add a day column with the appropriate day ranges
PSS$Day <- c(148:252)

# Convert to long format 
test2 <- PSS[,2:29] %>% pivot_longer(cols = -c(Day), names_to = "Year", values_to = "count") %>% as.data.frame()

# replace NA's with 0
test2[is.na(test2)]<-0

# Change Year to numeric
test2$Year <- as.numeric(test2$Year)

# Arrange data set by year
PSS_hist<- test2 %>% arrange((Year))

# make sure it worked
head(PSS_hist)
dim(PSS_hist)

