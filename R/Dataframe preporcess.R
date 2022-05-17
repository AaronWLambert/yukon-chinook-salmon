#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Purpose: To preprocess data frames for model
#
#   
#
#
#=================================================================================
#NOTES:
# This script will be where data is processed for entry into Ykon Inseason Forecast.R script
#
# 
# Next steps: 
# 
#
# PSS historical counts
# Data imported from  
# https://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareayukon.salmon_escapement
# Import as xlsx file (excel)

library(lubridate)
library(readxl)
test<- read_xlsx(file.path(dir.data,"Yukon Escapement Daily test2.xlsx"),skip = 3)
test$Day <- c(148:252)

test2<-test[,2:29] %>% pivot_longer(cols = -c(Day), names_to = "Year", values_to = "count") %>% as.data.frame()
test2[is.na(test2)]<-0
test2$Year <- as.numeric(test2$Year)
PSS_hist<- test2 %>% arrange((Year))
PSS_hist