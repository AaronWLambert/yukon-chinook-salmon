#**********************************************************************************
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Version: 1.0 to 3.4
# Purpose: To make yearly updates to data files:
#           GSI
#           EOS Can passage
#
#   1) Read in data
#   2) Preprocess data
#   3) Call to Stan model to generate inseason projection
#
#
#NOTES:
# This script is the general working script that is used to run single iterations 
#  of the model with a Stan file for all versions up to 3.4.
#   
#   This script is the rewritten code for leave-one-out retro testing
# Next steps: 
# 
#
# Packages #########################################################################
require(rstan)
require(bayesplot)
require(tidyverse)
require(mgcv)
require(ggthemes)
require(viridis)
require(shinystan)
require(lubridate)
require(reshape2)
require(dplyr)
require(tidybayes)


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


######### Import Data ###############
# Historical Canadian EOS reconstructed run
# CAN_hist <- readRDS(file.path(dir.data,"/EOS Reconstructed Canadian.RDS"))
# CAN_hist <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
# This is the reconstructed data from Curry for old reconstructed modeling procedure
CAN_hist <- readRDS(file.path(dir.data,"Can Origin Reconstructed 2Feb22.RDS"))

# PSS historical for Jun 1st = 152, to Aug 31 = 243
# Data obtained from ADFG website
# https://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareayukon.salmon_escapement
# NOTE!!!! See use "Dataframe preprocess.R" to get  most recent version.
#  If not using Dataframe preprocess.R, uncomment out the following to load relevant PSS historical data
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage 22Jun22.RDS"))

# Eagle passage
Eagle_hist <- readRDS(file = file.path(dir.data, "Eagle passage 23Nov22.RDS"))

# Back filled PSS matrix
# PSS_hist <- readRDS(file = file.path(dir.data,"PSS filled values 21oct22.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
# GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - 2021)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
# Read in historic preseason forecasts (2013 - current)
pf_hist <- readRDS(file.path(dir.data,"pf_ver3.1_14Aprl23_eagle+harvest.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs) and days not covered by GSI
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 21Mar23.RDS"))


# EOS Can_chinook ################################################################
# Add new row of observations for historic observed EOS Can-Chinook
head(CAN_hist)

# Year added

Year <- 2022

abund <- 13144

str(CAN_hist)
str(abund)
str(Year)
# Add to DF

CAN_hist[nrow(CAN_hist)+1,1] <- Year
CAN_hist[nrow(CAN_hist),2] <- abund

view(CAN_hist)

# saveRDS(CAN_hist, file = file.path(dir.data,"Can EOS Abund 3Mar23.RDS"))


# GSI Proportions ################################################################
tail(GSI_by_year,20)

# Fill in the relevant year, stratum, start and end days, and the observed prop and sd
GSI_add <- data.frame(year = rep(2022,3),
                      stratum = c(1:3),
                      startday = c(152, 174, 181),
                      endday = c(173, 180, 207),
                      propCan = c(0.67, 0.42, 0.35),
                      sd = c(NA,NA,NA)) # NOTE! The updated gsi is pulled from JTC report
                                        # No sd supplyed there, however, not used in model

# Add new df to previous years obs
GSI_by_year_updated <- rbind(GSI_by_year,GSI_add)

# saveRDS(GSI_by_year_updated, file = file.path(dir.data,"GSI by year unadj 4Apr23.RDS"))

# Preseason forecast
newYearPF <- data.frame("Year" = 2023, "mean" = 33,967)

Pf <- rbind(pf_hist,newYearPF)

saveRDS(object = Pf, file = file.path(dir.data, "pf_ver3.1_12June23_eagle+harvest.RDS"))

