#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.28.22
# 
# Purpose: This is a working script to generate figures for papers and presentations
library(tidyverse)
library(ggthemes)
library(wesanderson)
# Define Workflow Paths ============================================
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")

dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

######### Import Data ###############
# Historical Canadian EOS reconstructed run
# CAN_hist <- readRDS(file.path(dir.data,"/EOS Reconstructed Canadian.RDS"))
CAN_hist <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
str(CAN_hist)
# This is the reconstructed data from Curry for old reconstructed modeling procedure
# CAN_hist <- readRDS(file.path(dir.data,"Can Origin Reconstructed 2Feb22.RDS"))
# PSS historical for Jun 1st = 1, to Aug 31 = 92
# Years 1995 - 2020
# PSScum_hist <- readRDS(file.path(dir.data,"/clean_counts.rds"))

# Read in PSS non-cummulative with date in date format
PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))
# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
# pf_hist <- readRDS(file.path(dir.data,"inv_var_weighted_forcast_v3_Jan282022.RDS"))
# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))

# Read in harvest data for US portion of Yukon
usYukon <- readRDS(file = file.path(dir.data,"Alaska_Harvest_of_Yukon_River_Chinook_Salmon_1961_2020.RDS"))

# Read in Canadian Harvest of Yukon kings
canHarvest <- readRDS(file = (file.path(dir.data,"Canadian Harvest 1961_2020.RDS")))


# Get average harvest 
Can_avg <- canHarvest  %>% summarise(Mean2006 = mean(total)) 

pf_hist$Year <- as.double(pf_hist$Year)
# Join can_hist and pf_hist for plotting
full_abund <-full_join(CAN_hist,pf_hist, by = "Year")
# Plot to compare previous predictions to realized runsizes
ggplot(full_abund, aes(x = Year,
                     y = can.mean,
                     ymin = can.mean - can.sd,
                     ymax = can.mean + can.sd))+
  geom_point(fill = "red")+
  geom_errorbar()+
  geom_point( aes(x = Year, y = Weighted_Forecast), color = "red")

ggplot(pf_hist, aes(x = Year, ymin = Low/1000, ymax = High/1000, y = PostSeasonEstimate/1000))+
  geom_point(col = "red")+
  geom_errorbar()+
  labs(y = "1000's of Chinook Salmon")


# PLots for harvest
#
# Create Country variable for each df
canHarvest$Country <- "Canada"
usYukon$Country <- "USA"
canHarvest$Country <- as.factor(canHarvest$Country)
usYukon$Country <- as.factor(usYukon$Country)

longHarvest_US <- usYukon %>% pivot_longer(cols = -c(Year,Country)) %>% as.data.frame()
longHarvest_Can <- canHarvest%>% pivot_longer(cols = -c(Year,Country)) %>% as.data.frame()
longHarvest_Can
fullHarvestDF <- rbind(longHarvest_Can,longHarvest_US)

totalHarvDF <- fullHarvestDF %>% 
  group_by(Year,Country) %>% 
  summarise(TotalHarvest = sum(value)) %>% 
  as.data.frame()


ggplot(totalHarvDF, aes(x = Year, y = TotalHarvest/1000, fill = Country))+
  geom_col()+
  labs(y = "1,000's of Chinook Salmon",
       fill = "")+
  scale_fill_colorblind()+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank()
  )

# Average line calc
US_avg <- longHarvest_US %>% group_by(Year) %>% summarise(Total = sum(value))
US_avg <- mean(US_avg$Total)
ggplot(longHarvest_US, aes(x = Year, y = value/1000, fill = name))+
  geom_col()+
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  labs(y = "1,000's of Chinook Salmon",
       fill = "",
       title = "U.S. Harvest")+
  geom_hline(yintercept = US_avg/1000, linetype = 6, size = 2, color = "red")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle =  90)
          )+
  scale_x_continuous(breaks = seq(1961,2020,2))

# Get ave harvest for geom hline
Can_avg <- longHarvest_Can %>% group_by(Year) %>% summarise(Total = sum(value))
Can_avg <- mean(Can_avg$Total)
# Canadian harvest
ggplot(longHarvest_Can, aes(x = Year, y = value/1000, fill = name))+
  geom_col()+
  scale_fill_manual(values = wes_palette("IsleofDogs1"))+
  labs(y = "1,000's of Chinook Salmon",
       fill = "",
       title = "Can Harvest")+
  geom_hline(yintercept = Can_avg/1000, linetype = 6, size = 2, color = "red")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle =  90))+
  scale_x_continuous(breaks = seq(1961,2020,2))

ggplot(fullHarvestDF, aes(x = Year, y = value/1000, fill = factor(name)))+
  geom_col()+
  labs(x = "Year",
       y = "Number of Chinook Salmon",
       fill = "Harvest Catagory")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
  )+
  facet_wrap(~Country )










