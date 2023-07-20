# Aaron Lambert
# 2/21/2023
# Yukon Inseason Chinook Forecast
# Sea ICe
# https://cran.r-project.org/web/packages/heatwaveR/vignettes/OISST_preparation.htmlnegative

# NOAA OISST on ERDDAP server
# Define Workflow Paths ============================================

# set to working directory
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

# Objects used to save/load data, outputs, or stan/R scripts
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

all.timing.logistic <- read.csv(file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))
# Following scrape code from Robert W Schlegel and AJ Smit

# The packages we will use
library(dplyr) # A staple for modern data management in R
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing
library(ggpubr)


# The information for the NOAA OISST data
# rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# Note that there is also a version with lon values from 0 yo 360
rerddap::info(datasetid = "ncdcOisst21Agg", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# This function downloads and prepares data based on user provided start and end dates
OISST_sub_dl <- function(time_df){
  OISST_dat <- griddap(datasetx = "ncdcOisst21Agg", 
                       url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                       time = c(time_df$start, time_df$end), 
                       zlev = c(0, 0),
                       latitude = c(62, 64),
                       longitude = c(191,195),
                       fields = "ice")$data %>% 
    mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>%
    dplyr::rename(t = time) %>%
    select(longitude, latitude, t, ice) %>%
    na.omit()
}

# Date download range by start and end dates per year
dl_years <- data.frame(date_index = 1:4,
                       # start = as.Date(c("1982-01-01", "1990-01-01", 
                       #                   "1998-01-01", "2006-01-01", "2014-01-01")),
                       # end = as.Date(c("1989-12-31", "1997-12-31", 
                       #                 "2005-12-31", "2013-12-31", "2019-12-31")))
                       start = as.Date(c("1995-01-01", "2003-01-01", 
                                         "2011-01-01", "2020-01-01")),
                       end = as.Date(c("2002-12-31", "2010-12-31", 
                                       "2019-12-31", "2022-12-31")))

# Download all of the data with one nested request
# The time this takes will vary greatly based on connection speed
system.time(
  OISST_data <- dl_years %>% 
    group_by(date_index) %>% 
    group_modify(~OISST_sub_dl(.x)) %>% 
    ungroup() %>% 
    select(longitude, latitude, t, ice)
)



# Data analysis
yukon.ice  <- OISST_data %>% mutate(year = lubridate::year(t),
                                    month = lubridate::month(t),
                                    day = lubridate::day(t),
                                    jday = lubridate::yday(t))

norton.ice <- yukon.ice %>% 
  group_by(year, month, day) %>% 
  summarise(day.mean = mean(ice)) %>% 
  group_by(year,month) %>% 
  summarise(month.mean = mean(day.mean))


ggplot(norton.ice, aes(x = year, y = month.mean))+
  geom_point()+
  facet_wrap(~month)


mid.point <- all.timing.logistic[,c(2,5)]


norton.ice.mp <- left_join(norton.ice, mid.point)

norton.ice.mp <- norton.ice.mp[norton.ice.mp$year != 1996,]

ggplot(norton.ice.mp, aes(x = month.mean, y = mid))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~month, scales = "free_x")+
  stat_cor(method = "pearson", label.y = 182)+
  labs(x = "Mean Monthly PIC",
       y = "Observed Midpoint")+
  theme(text = element_text(size = 18))

# 




