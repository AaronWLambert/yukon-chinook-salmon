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

curr_year <- 2022
# Read in the file as xlsx
# Pilot Station Sonar Daily Passage
PSS <- as.data.frame(read_xlsx(file.path(dir.data,"ADFG PSS Daily/Yukon Escapement Daily Final 2022.xlsx"),skip = 3))

# Eagle Sonar Daily Passage
Eagle <- as.data.frame(read_xlsx(file.path(dir.data,"ADFG Eagle Daily/Yukon Escapement Daily Eagle 23Nov22.xlsx"),skip = 3))
# prior.df <- readRDS(file = file.path(dir.data, "normalprior.RDS"))

# Add a day column with the appropriate day ranges
PSS.year <- length(1996:curr_year)+2
PSS$Day <- c(148:252)
PSS.mat <- PSS[2:PSS.year]

eagle.year <- length(2005:curr_year)+2
Eagle$Day <- c(178:243)
eagle.mat <- Eagle[2:eagle.year]


# Use this section if filling in missing days at the beginning of the Run #########################
# replacing NAs with estimates for days up to July 19
# 
# saveRDS(prior.df, file = file.path(dir.data,"normalprior.RDS"))
# years <- prior.df$Year
# days <- 148:200
# for(y in 1:length(prior.df$Year)){
# # y <- 1
# 
#   sigma <- prior.df$sd[prior.df$Year== years[y]]
#   mu <- prior.df$mu[prior.df$Year == years[y]]
#   alpha <- prior.df$alpha[prior.df$Year == years[y]]
# 
#   for(d in 1:length(days)){
#     # d <- 1
#     day <- days[d]
#     if(is.na(PSS.mat[d,y])){
# 
#         est <- alpha*(1/(sigma*sqrt(2*pi)) * exp(-0.5 *((day - mu)/sigma)^2))
# 
#           PSS.mat[d,y] <- round(est,digits = 0)
#           }# end of if
# 
#       } # end of day loop
# 
# }# end of year loop
# 
# # Pivot Longer
# test2 <- PSS.mat %>% 
#   pivot_longer(cols = -c(Day), names_to = "Year", values_to = "count") %>% 
#   as.data.frame()
# 
# # replace NA's with 0
# test2[is.na(test2)]<-0
# 
# # Change Year to numeric
# test2$Year <- as.numeric(test2$Year)
# 
# # Arrange data set by year
# PSS_hist_filled<- test2 %>% arrange((Year))
# 
# # Save 
# saveRDS(object = PSS_hist_filled, file = file.path(dir.data, "PSS filled values 21oct22.RDS"))
# 


# Use this section if keeping NA's as zeros ##############################################
# PSS
# Convert to long format 
test2 <- PSS[,2:PSS.year] %>% 
  pivot_longer(cols = -c(Day), names_to = "Year", values_to = "count") %>% 
  as.data.frame()

# replace NA's with 0
test2[is.na(test2)]<-0

# Change Year to numeric
test2$Year <- as.numeric(test2$Year)

# Arrange data set by year
PSS_hist<- test2 %>% arrange((Year))

# make sure it worked
head(PSS_hist)
dim(PSS_hist)

# Clean up the environment
rm(test2,PSS)

## Uncoment to save PSS_hist if needed 
# saveRDS(object = PSS_hist, file = file.path(dir.data, "PSS passage Final 2022.RDS"))


# Eagle #############################################################################
# Convert to long format 
test3 <- Eagle[,2:eagle.year] %>% 
  pivot_longer(cols = -c(Day), names_to = "Year", values_to = "count") %>% 
  as.data.frame()

# replace NA's with 0
test3[is.na(test3)]<-0

# Change Year to numeric
test3$Year <- as.numeric(test3$Year)

# Arrange data set by year
eagle_hist<- test3 %>% arrange((Year))

# make sure it worked
head(eagle_hist)
dim(eagle_hist)

saveRDS(object = eagle_hist, file = file.path(dir.data, "Eagle passage 23Nov22.RDS"))
