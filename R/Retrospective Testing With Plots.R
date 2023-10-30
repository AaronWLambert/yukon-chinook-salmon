#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.14.22
# Version: This will run the selected model version in the function, and 
#           then generate plots and values.
# Purpose: This is a working script to generate all outputs and figures (at this time)
# 
#
#   1) Read in data
#   2) Call the function
#   3) Generate retrospective outputs such as figures and values.
#
#
#=================================================================================
# NOTES:
# 
#
# 
# Next steps: 

# 
#
#=================================================================================
library(rstan)
library(bayesplot)
library(tidyverse)
library(ggthemes)
library(viridis)
library(shinystan)
library(lubridate)
library(ggpubr)
library(gridExtra)
library(tidybayes)
library(wesanderson)
library(grid)
library(bbmle)

# Parralize for optimum model run time
rstan_options(auto_write = TRUE)
#
mc.cores = parallel::detectCores()
# mc.cores <-1

# Define Workflow Paths ============================================
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

# Functions ##################################################################################

# Functions to run model, create basic plots, and get day of year
source(file = file.path(dir.R,"model_run_function LOO3.r"))


# Function to get retrospective stats
source(file = file.path(dir.R,"Retro Function.R"))

######### Import Data ###############

# EOS Can-orig reconstructed counts (NEW METHOD)
CAN_hist_new <- readRDS(file.path(dir.data,"Canadian Passage RR 21Mar23.RDS"))

# This is the reconstructed data from Curry (Eagle + Harvest)
# CAN_hist_old <- readRDS(file.path(dir.data,"Can EOS Abund 3Mar23.RDS"))

# PSS Passage
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage Final 2022.RDS"))

# Back filled PSS matrix
# PSS_hist_filled <- readRDS(file = file.path(dir.data,"PSS filled values 21oct22.RDS"))

# Eagle passage
Eagle_hist <- readRDS(file = file.path(dir.data, "Eagle passage 23Nov22.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean <- readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
pf_hist <- readRDS(file.path(dir.data,"pf_ver3.1_14Aprl23_eagle+harvest.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 21Mar23.RDS"))

# Read in PSS observation error estimates
PSS_sd <- readRDS(file = file.path(dir.data,"PSS SD 1995_2021.RDS"))

# Shape parameters for fitting curves
logistic.all <- read.csv(file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))
# logistic.CAN <- read.csv(file = file.path(dir.output,"logistic curve parameters CAN Chinook 1995_2022.csv"))

normal.all <- read.csv(file = file.path(dir.output,"normal curve parameters All Chinook 1995_2022.csv"))
# normal.CAN <- read.csv(file = file.path(dir.output,"normal curve parameters CAN Chinook 1995_2022.csv"))

# SST May data
# norton.sst <- readRDS(file = file.path(dir.data,"norton.sst 6Apr23.RDS"))
norton.sst <- readRDS(file = file.path(dir.data,"SST Point -164875_62625 2000_2022.RDS"))

# Emmonak April air temp
Emmonak_Air_Temp <- readRDS(file = file.path(dir.data,"Emmonak Air Interpolated 2020 10Aug23.RDS"))

# Control Section ######################
model.version <- "6.prop.null"
# Range of years $$$$ to 2021
# myYear <- 2019

# Range of days 152 -243
# myDay <- 196

# MCMC Parameters
n.chains <- 4
n.iter <- 5000 #30000;#5e4
n.thin <- 2

# Days to use in retrospective testing runs ##############

# Test days used in full season run (every 5 days starting June 2)
testDays <- seq(from = 153, to = 243, by = 5)
# testDays <- seq(from = 198, to = 243, by = 5)

# Test days for runs starting June 12
# testDays_short <- seq(from = 153, to = 213, by = 5)

# Years included in retrospective testing

# Years included in full run (Used for retro with Eagle+Harvest Can-orig counts)
testYears <- c(2007:2022)
# testYears <- 2023

# Test Years used in retros for models only using recent PSS passage
# testYears_short <- c(2009:2010,2013:2021)

# Test years used for retro with New Can-orig recontructed counts
#  (currently only available to 2019 as of 11/5/22)
# testYears_short_new <- c(2007:2010,2013:2019)


# # Run the model for a single year and day of interest.
# # # # # #  Ajust the inputs for the model in the control section above
# model.output <-InSeasonProjection(model.version = model.version,
#                                   myYear = 2011,
#                                   myDay = 173,
#                                   n.chains = 4,
#                                   n.iter = 5000,
#                                   n.thin = n.thin,
#                                   CAN_hist = CAN_hist_new,
#                                   PSS_hist = PSS_hist,
#                                   Eagle_hist = Eagle_hist,
#                                   # GSI_mean = GSI_mean,
#                                   pf_hist = pf_hist,
#                                   GSI_by_year = GSI_by_year,
#                                   normal = FALSE,
#                                   logistic = TRUE,
#                                   startDayPSS = 148,
#                                   startYearPSS = 2005,
#                                   prior.df.log = logistic.all,
#                                   prior.df.norm = normal.all,
#                                   PSS_sd = PSS_sd
#                                     )


# Retrospective testing loop ######################

# # List to store outputs
outputList<-list()

# Increse memory limit for runs of model versions 4 and higher
# memory.limit(size = 60000)
options(warn = 1)
for(y in c(testYears)){
  for(d in c(testDays)){
    
    outputList[[paste("",y,"_",d, sep = "")]]<-InSeasonProjection(model.version = model.version,
                                                                  myYear = y,
                                                                  myDay = d,
                                                                  n.chains = n.chains,
                                                                  CAN_hist = CAN_hist_new,
                                                                  pf_hist = pf_hist,
                                                                  PSS_hist = PSS_hist,
                                                                  PSS_sd = PSS_sd,
                                                                  norton.sst = norton.sst, 
                                                                  Emmonak_Air_Temp = Emmonak_Air_Temp,
                                                                  n.thin = n.thin,
                                                                  n.iter = n.iter,
                                                                  GSI_by_year = GSI_by_year,
                                                                  Eagle_hist = Eagle_hist,
                                                                  normal = FALSE,
                                                                  logistic = FALSE,
                                                                  prior.df.log = logistic.all,
                                                                  prior.df.norm = normal.all,
                                                                  multiplier = 1,
                                                                  startDayPSS = 148,
                                                                  startYearPSS = 2005)
   print(paste("Day =",d,"Year =",y))
  } #dloop
  
  print(paste("Finally done with year",y))
} #yloop


# Save or read in outputlist #######################################

# Save output
saveRDS(object = outputList, file = file.path(dir.output, "Ver6propNull_final_27Sept23.RDS"))

# Read in outputs from retro testing 
# **Note** Increase memory limit if using versions 4 and higher
# memory.limit(size = 60000)
# outputlist_ver1 <- readRDS(file = file.path(dir.output,"Ver10_test_13Jan23.RDS")) This for days up to 213 and only to 2021

# outputlist_ver1 <- readRDS(file = file.path(dir.output, "Ver1_long_25Apr23.RDS")) # Day up to 243 and 2022

# outputlist_ver1RR <- readRDS(file = file.path(dir.output, "Ver1_RRCAN_long_27Apr23.RDS"))

# outputlist_ver1_T <- readRDS(file = file.path(dir.output,"Ver1_oldCAN_trunc_long_3May23.RDS"))

# outputlist_ver1RR_T <- readRDS(file = file.path(dir.output, "Ver1_rr_trunc_FINAL_5June23.RDS"))

outputlist_ver6propRR_T <- readRDS(file = file.path(dir.output, "Ver6prop2_rr_trunc_7Aug23.RDS"))

outputlist_ver1sstRR_T <- readRDS(file = file.path(dir.output, "Ver1sst_rr_trunc_2Aug23.RDS"))

outputlist_ver6402RR_T <- readRDS(file = file.path(dir.output, "Ver6402_rr_trunc_2Aug23.RDS"))

# outputlist_ver11RR_T <- readRDS(file = file.path(dir.output, "Ver11_rr_trunc_long_5May23.RDS"))
# 
# outputlist_ver101RR_T <- readRDS(file = file.path(dir.output, "Ver101_rr_trunc_long_5May23.RDS"))
# 
# outputlist_ver101ifRR_T <- readRDS(file = file.path(dir.output, "Ver101if_rr_trunc_long_5May23.RDS"))
# 
outputlist_ver102RR_T <- readRDS(file = file.path(dir.output, "Ver102_rr_trunc_long_5May23.RDS"))
# 
# outputlist_ver112RR_T <- readRDS(file = file.path(dir.output, "Ver112_rr_trunc_long_5May23.RDS"))
# 
# outputlist_ver6RR_T <- readRDS(file = file.path(dir.output, "Ver6_rr_trunc_long_5May23.RDS"))

# outputlist_ver2 <- readRDS(file = file.path(dir.output,"Ver2_long_25Apr23.RDS")) # With 2022 and till end Aug

outputlist_ver2RR <- readRDS(file = file.path(dir.output,"Ver2_RRCAN_long_1May23.RDS")) # With 2022 and till end Aug

# outputlist_ver2_T <- readRDS(file = file.path(dir.output,"Ver2_oldCAN_trunc_long_3May23.RDS")) # With 2022 and till end Aug

outputlist_ver2RR_T <- readRDS(file = file.path(dir.output,"Ver2_RRCAN_trunc_long_3May23.RDS")) # With 2022 and till end Aug

# outputlist_ver101 <- readRDS(file = file.path(dir.output,"Ver101_long_10Mar23.RDS"))

# outputlist_ver2c <- readRDS(file = file.path(dir.output,"Ver2c_long_28Feb23.RDS"))

# outputlist_ver2c <- readRDS(file = file.path(dir.output,"Ver2c_long_25Apr23.RDS")) # Day up to 243 and 2022

outputlist_ver2cRR <- readRDS(file = file.path(dir.output,"Ver2c_RRCAN_long_1May23.RDS")) # Day up to 243 and 2022

# outputlist_ver2c_TT <- readRDS(file = file.path(dir.output,"Ver2c_oldCAN_trunc_long_3May23.RDS")) # Day up to 243 and 2022

outputlist_ver2cRR_TT <- readRDS(file = file.path(dir.output,"Ver2c_RRCAN_trunc_long_3May23.RDS")) # Day up to 243 and 2022

# outputlist_ver2c1 <- readRDS(file = file.path(dir.output,"Ver2c1_bystrata_long_27Apr23.RDS")) # Day up to 243 and 2022 & mean across days

# outputlist_ver2c1strata <- readRDS(file = file.path(dir.output,"Ver2c1_long_25Apr23.RDS")) # Day up to 243 and 2022 & mean across days

# outputlist_ver2cSST <- readRDS(file = file.path(dir.output,"Ver2cSST_long_6Apr23.RDS"))

# outputlist_ver2c01 <- readRDS(file = file.path(dir.output,"Ver2c01_long_27Apr23.RDS")) # Day up to 243 and 2022

# outputlist_ver2c02 <- readRDS(file = file.path(dir.output,"Ver2c02_long_27Apr23.RDS")) # Day up to 243 and 2022

# outputlist_ver2c02sst <- readRDS(file = file.path(dir.output,"Ver2c02.sst_long_10Apr23.RDS")) 

# outputlist_ver2c01if <- readRDS(file = file.path(dir.output,"Ver2c01if_long_27Apr23.RDS"))# Day up to 243 and 2022

# outputlist_ver2c11 <- readRDS(file = file.path(dir.output,"Ver2c11_long_16Mar23.RDS"))

# outputlist_ver2csd<- readRDS(file = file.path(dir.output, "Ver2csd_long_10Mar23.RDS"))

# outputlist_ver2c <- readRDS(file = file.path(dir.output,"Ver2c_test_13Jan23.RDS"))

# outputlist_ver3 <- readRDS(file = file.path(dir.output,"Ver30_test_13Jan23.RDS"))

# outputlist_ver5 <- readRDS(file=file.path(dir.output, "Ver5_long_10Apr23.RDS"))


# outputlist_ver11 <- readRDS(file = file.path(dir.output,"Ver11_test_16Dec22.RDS"))

# outputlist_ver101 <- readRDS(file = file.path(dir.output,"Ver101_test_16Dec22.RDS"))



##### Calculations for retrospecitve testing #############################################
# Uses retrospective.function

# RetroList_ver1 <- retrospective.function(outputList = outputlist_ver1,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
#                                          pf = FALSE)

# RetroList_ver1_T <- retrospective.function(outputList = outputlist_ver1_T,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
#                                          pf = FALSE)



# RetroList_ver1RR <- retrospective.function(outputList = outputlist_ver1RR,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_new,
#                                          pf = FALSE)

RetroList_ver1RR_T <- retrospective.function(outputList = outputlist_ver1RR_T,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist_new,
                                           pf = FALSE)

RetroList_ver6RR_T <- retrospective.function(outputList = outputlist_ver6propRR_T,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             pf = FALSE)

RetroList_ver1sstRR_T <- retrospective.function(outputList = outputlist_ver1sstRR_T,
                                                testYears = testYears,
                                                testDays = testDays,
                                                CAN_hist = CAN_hist_new,
                                                pf = FALSE)

RetroList_ver6402RR_T <- retrospective.function(outputList = outputlist_ver6402RR_T,
                                                testYears = testYears,
                                                testDays = testDays,
                                                CAN_hist = CAN_hist_new,
                                                pf = FALSE)

RetroList_ver11RR_T <- retrospective.function(outputList = outputlist_ver11RR_T,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             pf = FALSE)

RetroList_ver101RR_T <- retrospective.function(outputList = outputlist_ver101RR_T,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             pf = FALSE)

RetroList_ver101ifRR_T <- retrospective.function(outputList = outputlist_ver101ifRR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

RetroList_ver102RR_T <- retrospective.function(outputList = outputlist_ver102RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

RetroList_ver112RR_T <- retrospective.function(outputList = outputlist_ver112RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

RetroList_ver6RR_T <- retrospective.function(outputList = outputlist_ver6RR_T,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             pf = FALSE)

# RetroList_ver11 <- retrospective.function(outputList = outputlist_ver11,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
#                                          pf = FALSE)
# 
# RetroList_ver101 <- retrospective.function(outputList = outputlist_ver101,
#                                           testYears = testYears,
#                                           testDays = testDays,
#                                           CAN_hist = CAN_hist_old,
#                                           pf = FALSE)

# RetroList_ver2 <- retrospective.function(outputList = outputlist_ver2,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
#                                          startYearRetro = 2007,
#                                          pf = FALSE)
# 
# RetroList_ver2_T <- retrospective.function(outputList = outputlist_ver2_T,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
#                                          startYearRetro = 2007,
#                                          pf = FALSE)
# 
# RetroList_ver2c <- retrospective.function(outputList = outputlist_ver2c,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
#                                          startYearRetro = 2007,
#                                          pf = FALSE)
# 
# RetroList_ver2c_T <- retrospective.function(outputList = outputlist_ver2c_TT,
#                                           testYears = testYears,
#                                           testDays = testDays,
#                                           CAN_hist = CAN_hist_old,
#                                           startYearRetro = 2007,
#                                           pf = FALSE)
# 
# 
# 
# RetroList_ver2RR <- retrospective.function(outputList = outputlist_ver2RR,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_new,
#                                          startYearRetro = 2007,
#                                          pf = FALSE)
# 
# RetroList_ver2RR_T <- retrospective.function(outputList = outputlist_ver2RR_T,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist_new,
#                                            startYearRetro = 2007,
#                                            pf = FALSE)
# 
# RetroList_ver2cRR <- retrospective.function(outputList = outputlist_ver2cRR,
#                                           testYears = testYears,
#                                           testDays = testDays,
#                                           CAN_hist = CAN_hist_new,
#                                           startYearRetro = 2007,
#                                           pf = FALSE)
# 
# RetroList_ver2cRR_T <- retrospective.function(outputList = outputlist_ver2cRR_TT,
#                                             testYears = testYears,
#                                             testDays = testDays,
#                                             CAN_hist = CAN_hist_new,
#                                             startYearRetro = 2007,
#                                             pf = FALSE)
# 
# RetroList_ver2c1 <- retrospective.function(outputList = outputlist_ver2c1,
#                                             testYears = testYears,
#                                             testDays = testDays,
#                                             CAN_hist = CAN_hist_old,
#                                             startYearRetro = 2007,
#                                             pf = FALSE)
# 
# RetroList_ver2c1_strata <- retrospective.function(outputList = outputlist_ver2c1strata,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist_old,
#                                            startYearRetro = 2007,
#                                            pf = FALSE)
# 
# RetroList_ver2cSST <- retrospective.function(outputList = outputlist_ver2cSST,
#                                           testYears = testYears,
#                                           testDays = testDays,
#                                           CAN_hist = CAN_hist_old,
#                                           startYearRetro = 2007,
#                                           pf = FALSE)
# 
# RetroList_ver2csd <- retrospective.function(outputList = outputlist_ver2csd,
#                                             testYears = testYears,
#                                             testDays = testDays,
#                                             CAN_hist = CAN_hist_old,
#                                             startYearRetro = 2007,
#                                             pf = FALSE)
# 
# RetroList_ver2c01 <- retrospective.function(outputList = outputlist_ver2c01,
#                                           testYears = testYears,
#                                           testDays = testDays,
#                                           CAN_hist = CAN_hist_old,
#                                           startYearRetro = 2007,
#                                           pf = FALSE)
# 
# RetroList_ver2c01if <- retrospective.function(outputList = outputlist_ver2c01if,
#                                             testYears = testYears,
#                                             testDays = testDays,
#                                             CAN_hist = CAN_hist_old,
#                                             startYearRetro = 2007,
#                                             pf = FALSE)
# 
# RetroList_ver2c02 <- retrospective.function(outputList = outputlist_ver2c02,
#                                               testYears = testYears,
#                                               testDays = testDays,
#                                               CAN_hist = CAN_hist_old,
#                                               startYearRetro = 2007,
#                                               pf = FALSE)
# 
# RetroList_ver2c02sst <- retrospective.function(outputList = outputlist_ver2c02sst,
#                                             testYears = testYears,
#                                             testDays = testDays,
#                                             CAN_hist = CAN_hist_old,
#                                             startYearRetro = 2007,
#                                             pf = FALSE)
# 
# RetroList_ver2c11 <- retrospective.function(outputList = outputlist_ver2c11,
#                                             testYears = testYears,
#                                             testDays = testDays,
#                                             CAN_hist = CAN_hist_old,
#                                             startYearRetro = 2007,
#                                             pf = FALSE)
# 
# RetroList_ver266<- retrospective.function(outputList = outputlist_ver266,
#                                          testYears = testYears,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
#                                          pf = FALSE)

# RetroList_ver2_short <- retrospective.function(outputList = outputlist_ver2_oldCanHist_short,
#                                          testYears = testYears_short,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_old,
# #                                          pf = FALSE)
# 
# RetroList_ver2_new <- retrospective.function(outputList = outputlist_ver2_newCanHist,
#                                          testYears = testYears_short_new,
#                                          testDays = testDays,
#                                          CAN_hist = CAN_hist_new,
#                                          pf = FALSE)
# 
# RetroList_ver2_new_filled <- retrospective.function(outputList = outputlist_ver2_newCanHist_filled,
#                                              testYears = testYears_short_new,
#                                              testDays = testDays,
#                                              CAN_hist = CAN_hist_new,
#                                              pf = FALSE)
# 
# RetroList_ver2_old_filled <- retrospective.function(outputList = outputlist_ver2_oldCanHist_filled,
#                                                     testYears = testYears_short_new,
#                                                     testDays = testDays,
#                                                     CAN_hist = CAN_hist_old,
#                                                     pf = FALSE)

# RetroList_ver2.1 <- retrospective.function(outputList = outputlist_ver2.1_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist,
#                                            pf = FALSE)

# RetroList_ver2.2 <- retrospective.function(outputList = outputlist_ver2.2_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist,
#                                            pf = FALSE)

# RetroList_ver201_old <- retrospective.function(outputList = outputlist_ver201_oldCanHist[1:143],
#                                                     testYears = testYears_short_new,
#                                                     testDays = testDays,
#                                                     CAN_hist = CAN_hist_old,
#                                                     pf = FALSE)
# 
# RetroList_ver201_new <- retrospective.function(outputList = outputlist_ver201_newCanHist,
#                                                testYears = testYears_short_new,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_new,
#                                                pf = FALSE)
# 
# RetroList_ver211_old <- retrospective.function(outputList = outputlist_ver211_oldCanHist,
#                                                testYears = testYears_short_new,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_old,
#                                                startYearRetro = 2007,
#                                                endYearRetro = 2019,
#                                                pf = FALSE)
# 
# RetroList_ver211_new <- retrospective.function(outputList = outputlist_ver211_newCanHist,
#                                                testYears = testYears_short_new,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_new,
#                                                startYearRetro = 2007,
#                                                endYearRetro = 2019,
#                                                pf = FALSE)
# 
# RetroList_ver2110_old <- retrospective.function(outputList = outputlist_ver2110_oldCanHist,
#                                                testYears = testYears[5:13],
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_old,
#                                                startYearRetro = 2013,
#                                                pf = FALSE)
# 
# RetroList_ver2110_new <- retrospective.function(outputList = outputlist_ver2110_newCanHist,
#                                                 testYears = testYears_short_new[3:11],
#                                                 testDays = testDays,
#                                                 CAN_hist = CAN_hist_old,
#                                                 startYearRetro = 2009,
#                                                 pf = FALSE)
# 
# RetroList_ver211_new <- retrospective.function(outputList = outputlist_ver211_newCanHist,
#                                                testYears = testYears_short_new,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_new,
#                                                pf = FALSE)

# RetroList_ver3.0 <- retrospective.function(outputList = outputlist_ver3,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist_old,
#                                            pf = FALSE)

# RetroList_ver3.2 <- retrospective.function(outputList = outputlist_ver3.2_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist_old,
#                                            startYearRetro = 2007,
#                                            pf = FALSE)

# RetroList_ver3.3 <- retrospective.function(outputList = outputlist_ver3.3_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist,
#                                            pf = FALSE)

# RetroList_ver3.4 <- retrospective.function(outputList = outputlist_ver3.4_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist,
#                                            pf = FALSE)

# RetroList_ver3.5 <- retrospective.function(outputList = outputlist_ver3.5_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist,
#                                            pf = FALSE)

# RetroList_ver4.0 <- retrospective.function(outputList = outputlist_ver4.0_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays_short,
#                                            CAN_hist = CAN_hist,
#                                            pf = FALSE)
# 
# RetroList_ver4.1 <- retrospective.function(outputList = outputlist_ver4.1_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist_old,
#                                            pf = FALSE)

# RetroList_ver4.2 <- retrospective.function(outputList = outputlist_ver4.2_oldCanHist,
#                                            testYears = testYears,
#                                            testDays = testDays_short,
#                                            CAN_hist = CAN_hist,
#                                            pf = FALSE)


RetroList_ver5.0 <- retrospective.function(outputList = outputlist_ver5,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist_old,
                                           pf = FALSE)


RetroList_pf_old <- retrospective.function(outputList = outputlist_ver1,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist_old,
                                           pf_hist = pf_hist,
                                           startYearRetro = 2007,
                                           endYearRetro = 2022,
                                           pf = TRUE)

RetroList_pf_new <- retrospective.function(outputList = outputlist_ver1RR_T,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist_new,
                                           pf_hist = pf_hist,
                                           startYearRetro = 2007,
                                           endYearRetro = 2022,
                                           pf = TRUE)


# RMSE

# Extract RMSE by day into dataframe for each version and posterior of interest
# rmseDF_ver1 <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver1$RMSE_by_day_vect,
#                           "Version" = "ver1")
# 
# rmseDF_ver1RR <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver1RR$RMSE_by_day_vect,
#                           "Version" = "ver1 RR")
# 
# rmseDF_ver1_T <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver1_T$RMSE_by_day_vect,
#                           "Version" = "ver1 Trunc")

rmseDF_ver1RR_T <- data.frame("Day" = testDays,
                            "RMSE" = RetroList_ver1RR_T$RMSE_by_day_vect,
                            "Version" = "ver1 RR Trunc")

rmseDF_ver6RR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver6RR_T$RMSE_by_day_vect,
                              "Version" = "ver6 RR Trunc")

rmseDF_ver1sstRR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver1sstRR_T$RMSE_by_day_vect,
                              "Version" = "ver1 SST RR Trunc")

rmseDF_ver6402RR_T <- data.frame("Day" = testDays,
                                 "RMSE" = RetroList_ver6402RR_T$RMSE_by_day_vect,
                                 "Version" = "ver6402 RR Trunc")

rmseDF_ver11RR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver11RR_T$RMSE_by_day_vect,
                              "Version" = "ver11 RR Trunc")

rmseDF_ver101RR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver101RR_T$RMSE_by_day_vect,
                              "Version" = "ver101 RR Trunc")

rmseDF_ver101ifRR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver101ifRR_T$RMSE_by_day_vect,
                                "Version" = "ver101 if RR Trunc")

rmseDF_ver102RR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver102RR_T$RMSE_by_day_vect,
                              "Version" = "ver102 RR Trunc")

rmseDF_ver112RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver112RR_T$RMSE_by_day_vect,
                                "Version" = "ver112 RR Trunc")

rmseDF_ver6RR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver6RR_T$RMSE_by_day_vect,
                              "Version" = "ver6 RR Trunc")

# rmseDF_ver11 <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver11$RMSE_by_day_vect,
#                           "Version" = "ver1.1")
# 
# rmseDF_ver101 <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver101$RMSE_by_day_vect,
#                            "Version" = "ver1.0.1")
# 
# rmseDF_ver2 <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver2$RMSE_by_day_vect,
#                           "Version" = "ver2")
# 
# rmseDF_ver2RR <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver2RR$RMSE_by_day_vect,
#                           "Version" = "ver2 RR")
# 
# rmseDF_ver2_T <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver2_T$RMSE_by_day_vect,
#                           "Version" = "ver2 trunc")
# 
# rmseDF_ver2RR_T <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_ver2RR_T$RMSE_by_day_vect,
#                             "Version" = "ver2 RR Trunc")
# 
# rmseDF_ver21 <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver21$RMSE_by_day_vect,
#                           "Version" = "ver2.1")
# 
# rmseDF_ver2c <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver2c$RMSE_by_day_vect,
#                           "Version" = "ver2.c")
# 
# rmseDF_ver2c_T <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver2c_T$RMSE_by_day_vect,
#                            "Version" = "ver2.c Trunc")
# 
# rmseDF_ver2cRR <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver2cRR$RMSE_by_day_vect,
#                            "Version" = "ver2.c RR")
# 
# rmseDF_ver2c_T <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver2c_T$RMSE_by_day_vect,
#                            "Version" = "ver2.c Trunc")
# 
# rmseDF_ver2cRR_T <- data.frame("Day" = testDays,
#                              "RMSE" = RetroList_ver2cRR_T$RMSE_by_day_vect,
#                              "Version" = "ver2.c RR Trunc")
# 
# rmseDF_ver2c1 <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver2c1$RMSE_by_day_vect,
#                            "Version" = "ver2.c.1")
# 
# rmseDF_ver2c1_strata <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_ver2c1_strata$RMSE_by_day_vect,
#                             "Version" = "ver2.c.1.strata")
# 
# rmseDF_ver2cSST <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver2cSST$RMSE_by_day_vect,
#                            "Version" = "ver2.c.SST")
# 
# rmseDF_ver2csd <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver2csd$RMSE_by_day_vect,
#                            "Version" = "ver2.c.sd")
# 
# rmseDF_ver2c01 <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver2c01$RMSE_by_day_vect,
#                            "Version" = "ver2.c.0.1")
# 
# rmseDF_ver2c02 <- data.frame("Day" = testDays,
#                              "RMSE" = RetroList_ver2c02$RMSE_by_day_vect,
#                              "Version" = "ver2.c.0.2")
# 
# rmseDF_ver2c02sst <- data.frame("Day" = testDays,
#                              "RMSE" = RetroList_ver2c02sst$RMSE_by_day_vect,
#                              "Version" = "ver2.c.0.2.sst")
# 
# rmseDF_ver2c01if <- data.frame("Day" = testDays,
#                              "RMSE" = RetroList_ver2c01if$RMSE_by_day_vect,
#                              "Version" = "ver2.c.0.1.if")
# 
# rmseDF_ver2c11 <- data.frame("Day" = testDays,
#                              "RMSE" = RetroList_ver2c11$RMSE_by_day_vect,
#                              "Version" = "ver2.c.1.1")

# rmseDF_ver2_new_filled <- data.frame("Day" = testDays,
#                               "RMSE" = RetroList_ver2_new_filled$RMSE_by_day_vect,
#                               "Version" = "ver2 new filled")
# 
# rmseDF_ver2_old_filled <- data.frame("Day" = testDays,
#                                      "RMSE" = RetroList_ver2_old_filled$RMSE_by_day_vect,
#                                      "Version" = "ver2 old filled")

# rmseDF_ver2_dyn <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver2_dyn$RMSE_by_day_vect,
#                           "Version" = "ver2 dyn")

# rmseDF_ver2_short <- data.frame("Day" = testDays, 
#                           "RMSE" = RetroList_ver2_short$RMSE_by_day_vect,
#                           "Version" = "ver2 short")

# rmseDF_ver201 <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver201_old$RMSE_by_day_vect,
#                           "Version" = "ver201")
# 
# rmseDF_ver201_new <- data.frame("Day" = testDays,
#                           "RMSE" = RetroList_ver201_new$RMSE_by_day_vect,
#                           "Version" = "ver201 new")
# 
# rmseDF_ver211 <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_ver211_old$RMSE_by_day_vect,
#                             "Version" = "ver211")
# 
# rmseDF_ver211_new <- data.frame("Day" = testDays,
#                                 "RMSE" = RetroList_ver211_new$RMSE_by_day_vect,
#                                 "Version" = "ver211 new")
# 
# rmseDF_ver2110 <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_ver2110_old$RMSE_by_day_vect,
#                             "Version" = "ver2110")
# 
# rmseDF_ver2110_new <- data.frame("Day" = testDays,
#                              "RMSE" = RetroList_ver2110_new$RMSE_by_day_vect,
#                              "Version" = "ver2110 new")
# 
# rmseDF_ver2.1 <- data.frame("Day" = testDays, 
#                             "RMSE" = RetroList_ver2.1$RMSE_by_day_vect,
#                             "Version" = "ver2.1")

# rmseDF_ver2.2 <- data.frame("Day" = testDays, 
#                             "RMSE" = RetroList_ver2.2$RMSE_by_day_vect,
#                             "Version" = "ver2.2")

# rmseDF_ver3.0 <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_ver3.0$RMSE_by_day_vect,
#                             "Version" = "ver3.0")
# 
# rmseDF_ver3.2 <- data.frame("Day" = testDays, 
#                             "RMSE" = RetroList_ver3.2$RMSE_by_day_vect,
#                             "Version" = "ver3.2")
# 
# rmseDF_ver3.3 <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_ver3.3$RMSE_by_day_vect,
#                             "Version" = "ver3.3")
# 
# # rmseDF_ver3.4 <- data.frame("Day" = testDays, 
# #                             "RMSE" = RetroList_ver3.4$RMSE_by_day_vect,
# #                             "Version" = "ver3.4")
# 
# # rmseDF_ver3.5 <- data.frame("Day" = testDays, 
# #                             "RMSE" = RetroList_ver3.5$RMSE_by_day_vect,
# #                             "Version" = "ver3.5")
# 
# # rmseDF_ver4.0 <- data.frame("Day" = testDays_short, 
# #                             "RMSE" = RetroList_ver4.0$RMSE_by_day_vect,
# #                             "Version" = "ver4.0")
# 
# rmseDF_ver4.1 <- data.frame("Day" = testDays, 
#                             "RMSE" = RetroList_ver4.1$RMSE_by_day_vect,
#                             "Version" = "ver4.1")
# 
# # rmseDF_ver4.2 <- data.frame("Day" = testDays_short, 
# #                             "RMSE" = RetroList_ver4.2$RMSE_by_day_vect,
# #                             "Version" = "ver4.2")

# rmseDF_ver5 <- data.frame("Day" = testDays,
#                            "RMSE" = RetroList_ver5.0$RMSE_by_day_vect,
#                            "Version" = "ver5.0")
# 
rmseDF_pf_new <- data.frame("Day" = testDays,
                        "RMSE" = RetroList_pf_new$RMSE_by_day_vect,
                        "Version" = "PF RR")

# rmseDF_pf_old <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_pf_old$RMSE_by_day_vect,
#                             "Version" = "PF Eagle+Harv")

# Combine data frames into one data.frame for plotting
full_rmseDF <- rbind(
                     # rmseDF_ver1,
                     # rmseDF_ver1RR,
                     # rmseDF_ver1_T,
                     rmseDF_ver1RR_T,
                     # rmseDF_ver11RR_T,
                     # rmseDF_ver101RR_T,
                     # rmseDF_ver101ifRR_T,
                     # rmseDF_ver102RR_T,
                     # rmseDF_ver112RR_T,
                     rmseDF_ver6RR_T,
                     rmseDF_ver1sstRR_T,
                     rmseDF_ver6402RR_T,
                     # rmseDF_ver11,
                     # rmseDF_ver101,
                     # rmseDF_ver2,
                     # rmseDF_ver2_T,
                     # rmseDF_ver2c,
                     # rmseDF_ver2c_T,
                     # rmseDF_ver2RR,
                     # rmseDF_ver2RR_T,
                     # rmseDF_ver21,
                     # rmseDF_ver2c,
                     # rmseDF_ver2cRR,
                     # rmseDF_ver2cRR_T,
                     # rmseDF_ver2c1,
                     # rmseDF_ver2c1_strata,
                     # rmseDF_ver2c01,
                     # rmseDF_ver2c01if,
                     # rmseDF_ver2c02,
                     # rmseDF_ver2c02sst,
                     # rmseDF_ver2cSST,
                     # rmseDF_ver2c11,
                     # rmseDF_ver2csd,
                     # rmseDF_ver2_new,
                     # rmseDF_ver2_new_filled,
                     # rmseDF_ver2_old_filled,
                     # rmseDF_ver2_dyn,
                     # rmseDF_ver2_short,
                     # rmseDF_ver201,
                     # rmseDF_ver201_new,
                     # rmseDF_ver211,
                     # rmseDF_ver211_new,
                     # rmseDF_ver2110,
                     # rmseDF_ver2110_new,
                     # rmseDF_ver2.1,
                     # rmseDF_ver2.2,
                     # rmseDF_ver3.0,
                     # rmseDF_ver3.2,
                     # rmseDF_ver3.3,
                     # rmseDF_ver3.4,
                     # rmseDF_ver3.5,
                     # rmseDF_ver4.0,
                     # rmseDF_ver4.1,
                     # rmseDF_ver5,
                     # rmseDF_pf_old
                     rmseDF_pf_new
                     )

# Make version a factor
full_rmseDF$Version <- as.factor(full_rmseDF$Version)

# Table for papers
# RMSE_table<-as.data.frame(pivot_wider(full_rmseDF,names_from = Version, values_from = RMSE))
# RMSE_table_round <- format(RMSE_table, digits = 2)

# write.table(x = RMSE_table_round, file = file.path(dir.output,"RMSE Table 15Jul22"), row.names = F)

# Save to PDF 
# **Uncomment out to save plots after this line to PDF
# pdf(file = file.path(dir.output,"retro_plots_29Sept22.pdf"))

# RMSE plot
ggplot(full_rmseDF, aes(x = Day, y = RMSE, col = Version, shape= Version, alpha = 0.7))+
  # geom_col(position = "dodge",width = 3) +
  geom_point(size = 5)+
  geom_line(aes(group = Version))+
  # labs(fill = "", col = "")+
  # geom_col(data = short_rmseDF, position = "dodge", width = 3)+
  coord_cartesian(ylim = c(
    min(full_rmseDF$RMSE)-500,
    max(full_rmseDF$RMSE)+500
    # 10000,15000
                           
                           ))+
  # scale_colour_manual(name = "",
  # #                     # labels = c("PF","Ver 2.c", "Ver 2.c.0.1", "Ver 2.c.sd"),
  #                     # labels = c("PF","No Eagle", "Eagle Prop Est", "GSI No Eagle","GSI & Eagle Prop Est"),
  #                     labels = c("PF","Ver 1", "SST Est Prop"),
  #                     values = c("red","gold","green","blue", "purple")) +
  # scale_shape_manual(name = "",
  #                    # labels = c("PF","No Eagle", "Eagle Prop Est", "GSI No Eagle","GSI & Eagle Prop Est"),
  # 
  #                    labels = c("PF","Ver 1", "SST Est Prop"),
  #                    values = c(16,17,18,19,20))+
  # # scale_fill_manual(values = wes_palette("IsleofDogs1"),
  #                   labels = c("PF (New; SS Recon.)",
  #                              "PF (Old; Eagle+ Harvest)",
  #                              "2.1 (New; SS Recon.)",
  #                              "2.1 (Old; Eagle + Harvest)"))+
  # labels = c("Ver 1", "Ver 1-PF", "Ver 2", "Ver2-PF"))+
  guides(alpha = "none")+
  scale_x_continuous(breaks = c(testDays), labels = c(
                                                 "June 2",
                                                 "June 7",
                                                 "June 12",
                                                 "June 17",
                                                 "June 22",
                                                 "June 27",
                                                 "July 2",
                                                 "July 7",
                                                 "July 12",
                                                 "July 17",
                                                 "July 22",
                                                 "July 27",
                                                 "Aug 1", #,
                                                 "Aug 6",
                                                 "Aug 11",
                                                 "Aug 16",
                                                 "Aug 21",
                                                 "Aug 26",
                                                 "Aug 31"
                                                 ))+
  # facet_wrap(~)+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18),
        axis.text = element_text(angle = 90))




# MAPE ####################################################################

# Create MAPE dataframe
mapeDF_ver1 <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_ver1$MAPE_vect,
                          "Version" = "Ver 1")

mapeDF_ver1RR <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_ver1RR$MAPE_vect,
                          "Version" = "Ver 1 RR")

mapeDF_ver1RR_T <- data.frame("Day" = testDays,
                            "MAPE" = RetroList_ver1RR_T$MAPE_vect,
                            "Version" = "Ver 1 RR Trunc")

mapeDF_ver6RR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver6RR_T$MAPE_vect,
                              "Version" = "Ver 6 RR Trunc")

mapeDF_ver1sstRR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver1sstRR_T$MAPE_vect,
                              "Version" = "Ver 1sst RR Trunc")

mapeDF_ver6402RR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver6402RR_T$MAPE_vect,
                              "Version" = "Ver 6402 RR Trunc")

mapeDF_ver11 <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_ver11$MAPE_vect,
                          "Version" = "Ver 1.1")

mapeDF_ver101 <- data.frame("Day" = testDays,
                           "MAPE" = RetroList_ver101$MAPE_vect,
                           "Version" = "Ver 1.0.1")

mapeDF_ver2 <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_ver2$MAPE_vect,
                          "Version" = "Ver 2")

mapeDF_ver2RR <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_ver2RR$MAPE_vect,
                          "Version" = "Ver 2 RR")

mapeDF_ver21<- data.frame("Day" = testDays,
                         "MAPE" = RetroList_ver21$MAPE_vect,
                         "Version" = "Ver 2.1")

mapeDF_ver2c <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_ver2c$MAPE_vect,
                          "Version" = "Ver 2c")

mapeDF_ver2cRR <- data.frame("Day" = testDays,
                           "MAPE" = RetroList_ver2cRR$MAPE_vect,
                           "Version" = "Ver 2c RR")

mapeDF_ver2c1 <- data.frame("Day" = testDays,
                           "MAPE" = RetroList_ver2c1$MAPE_vect,
                           "Version" = "Ver 2c1")

mapeDF_ver2c1_strata <- data.frame("Day" = testDays,
                            "MAPE" = RetroList_ver2c1_strata$MAPE_vect,
                            "Version" = "Ver 2c1 strata")

mapeDF_ver2cSST <- data.frame("Day" = testDays,
                           "MAPE" = RetroList_ver2cSST$MAPE_vect,
                           "Version" = "Ver 2cSST")

mapeDF_ver2csd <- data.frame("Day" = testDays,
                           "MAPE" = RetroList_ver2csd$MAPE_vect,
                           "Version" = "Ver 2csd")

mapeDF_ver2c01 <- data.frame("Day" = testDays,
                           "MAPE" = RetroList_ver2c01$MAPE_vect,
                           "Version" = "Ver 2c01")

mapeDF_ver2c01if <- data.frame("Day" = testDays,
                             "MAPE" = RetroList_ver2c01if$MAPE_vect,
                             "Version" = "Ver 2c01if")

mapeDF_ver2c02 <- data.frame("Day" = testDays,
                             "MAPE" = RetroList_ver2c02$MAPE_vect,
                             "Version" = "Ver 2c02")

mapeDF_ver2c02sst <- data.frame("Day" = testDays,
                             "MAPE" = RetroList_ver2c02sst$MAPE_vect,
                             "Version" = "Ver 2c02 sst")

mapeDF_ver2c11 <- data.frame("Day" = testDays,
                             "MAPE" = RetroList_ver2c11$MAPE_vect,
                             "Version" = "Ver 2c11")
# 
# mapeDF_ver2_new_filled <- data.frame("Day" = testDays, 
#                               "MAPE" = RetroList_ver2_new_filled$MAPE_vect,
#                               "Version" = "Ver 2 new filled")
# 
# mapeDF_ver2_old_filled <- data.frame("Day" = testDays, 
#                                      "MAPE" = RetroList_ver2_old_filled$MAPE_vect,
#                                      "Version" = "Ver 2 old filled")
# 
# mapeDF_ver2_dyn <- data.frame("Day" = testDays, 
#                               "MAPE" = RetroList_ver2_dyn$MAPE_vect,
#                               "Version" = "Ver 2 dyn")
# 
# mapeDF_ver201 <- data.frame("Day" = testDays,
#                             "MAPE" = RetroList_ver201_old$MAPE_vect,
#                             "Version" = "Ver 2.0.1")
# 
# mapeDF_ver201_new <- data.frame("Day" = testDays,
#                             "MAPE" = RetroList_ver201_new$MAPE_vect,
#                             "Version" = "Ver 2.0.1 new")
# 
# mapeDF_ver211 <- data.frame("Day" = testDays,
#                             "MAPE" = RetroList_ver211_old$MAPE_vect,
#                             "Version" = "Ver 2.1.1")
# 
# mapeDF_ver211_new <- data.frame("Day" = testDays,
#                                 "MAPE" = RetroList_ver211_new$MAPE_vect,
#                                 "Version" = "Ver 2.1.1 new")
# 
# mapeDF_ver2110 <- data.frame("Day" = testDays,
#                             "MAPE" = RetroList_ver2110_old$MAPE_vect,
#                             "Version" = "Ver 2.1.1.0")
# 
# # mapeDF_ver2.2 <- data.frame("Day" = testDays, 
# #                             "MAPE" = RetroList_ver2.2$MAPE_vect,
# #                             "Version" = "Ver 2.2")
# 
mapeDF_ver3.0 <- data.frame("Day" = testDays,
                            "MAPE" = RetroList_ver3.0$MAPE_vect,
                            "Version" = "Ver 3.0")
# 
# mapeDF_ver3.2 <- data.frame("Day" = testDays, 
#                             "MAPE" = RetroList_ver3.2$MAPE_vect,
#                             "Version" = "Ver 3.2")
# 
# mapeDF_ver3.3 <- data.frame("Day" = testDays,
#                             "MAPE" = RetroList_ver3.3$MAPE_vect,
#                             "Version" = "Ver 3.3")
# 
# # mapeDF_ver3.4 <- data.frame("Day" = testDays, 
# #                             "MAPE" = RetroList_ver3.4$MAPE_vect,
# #                             "Version" = "Ver 3.4")
# 
# # mapeDF_ver3.5 <- data.frame("Day" = testDays, 
# #                             "MAPE" = RetroList_ver3.5$MAPE_vect,
# #                             "Version" = "Ver 3.5")
# 
# # mapeDF_ver4.0 <- data.frame("Day" = testDays_short, 
# #                             "MAPE" = RetroList_ver4.0$MAPE_vect,
# #                             "Version" = "Ver 4.0")
# 
# mapeDF_ver4.1 <- data.frame("Day" = testDays, 
#                             "MAPE" = RetroList_ver4.1$MAPE_vect,
#                             "Version" = "Ver 4.1")
# 
# # mapeDF_ver4.2 <- data.frame("Day" = testDays, 
# #                             "MAPE" = RetroList_ver4.2$MAPE_vect,
# #                             "Version" = "Ver 4.2")
# 
# mapeDF_PF_new <- data.frame("Day" = testDays, 
#                         "MAPE" = RetroList_pf_new$MAPE_vect,
#                         "Version" = "PF new")

mapeDF_ver5 <- data.frame("Day" = testDays,
                               "MAPE" = RetroList_ver5.0$MAPE_vect,
                               "Version" = "Ver 5.0")

mapeDF_PF_old <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_pf_old$MAPE_vect,
                            "Version" = "PF")

mapeDF_PF_new <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_pf_new$MAPE_vect,
                            "Version" = "PF RR")

# Combine into one data-frame
full_MAPE.df <- rbind(
                      # mapeDF_ver1,
                      # mapeDF_ver1RR,
                      mapeDF_ver1RR_T,
                      mapeDF_ver1sstRR_T,
                      mapeDF_ver6402RR_T,
                      # mapeDF_ver11,
                      # mapeDF_ver101,
                      # mapeDF_ver2,
                      # mapeDF_ver2RR,
                      # mapeDF_ver21,
                      # mapeDF_ver2c,
                      # mapeDF_ver2cRR,
                      # mapeDF_ver2c1,
                      # mapeDF_ver2c1_strata,
                      # mapeDF_ver2cSST,
                      # mapeDF_ver2c01,
                      # mapeDF_ver2c02,
                      # mapeDF_ver2c02sst,
                      # mapeDF_ver2c01if,
                      # mapeDF_ver2c11,
                      # mapeDF_ver2csd,
                      # mapeDF_ver2_dyn,
                      # mapeDF_ver2_new_filled,
                      # mapeDF_ver2_old_filled,
                      # mapeDF_ver201,
                      # mapeDF_ver201_new,
                      # mapeDF_ver211,
                      # mapeDF_ver211_new,
                      # mapeDF_ver2110,
                      # mapeDF_ver2.1,
                      # mapeDF_ver2.2,
                      # mapeDF_ver3.0,
                      # mapeDF_ver3.2,
                      # mapeDF_ver3.3,
                      # mapeDF_ver3.4,
                      # mapeDF_ver3.5,
                      # mapeDF_ver4.0,
                      # mapeDF_ver4.1,
                      # mapeDF_ver4.2,
                      # mapeDF_PF_new,
                      # mapeDF_ver5,
                      mapeDF_ver6RR_T,
                      mapeDF_PF_new
                      )



# MAPE Table for papers
# MAPE_table<-as.data.frame(pivot_wider(full_MAPE.df,names_from = Version, values_from = MAPE))
# MAPE_table_round <- format(MAPE_table, digits = 3)

# write.table(x = RMSE_table_round, file = file.path(dir.output,"RMSE Table 15Jul22"),
#             row.names = F)


# MAPE plot
ggplot(full_MAPE.df, aes(x = Day,
                         y = MAPE*100, 
                         col = Version,
                         shape = Version))+
  # geom_col(position = "dodge",width = 3) + 
  
  geom_point(size = 3)+
  geom_line(aes(group = Version))+
  # geom_line(aes(col = Version), size = 2)+
  # labs(fill = "",
  #      y = "MAPE (%)")+
  coord_cartesian(ylim = c(min(full_MAPE.df$MAPE*100)-2,max(full_MAPE.df$MAPE*100)+2))+
  ylab('MAPE %')+
  # labs(fill = "", col = "")+
  # geom_col(data = short_rmseDF, position = "dodge", width = 3)+
  # scale_colour_manual(name = "",
  #                     # labels = c("PF","Ver 2.c", "Ver 2.c.0.1", "Ver 2.c.sd"),
  #                     # labels = c("PF","No GSI", "GSI across days", "GSI by strata"),
  #                     labels = c("PF","No Eagle", "Eagle Regress", "Eagle Regress Back-half Season","Eagle Prop Est"),
  #                     
  #                     values = c("red","blue", "green", "gold", "brown")) +
  # scale_shape_manual(name = "",
  #                    # labels = c("PF","Ver 2.c", "Ver 2.c.0.1", "Ver 2.c.sd"),
  #                    labels = c("PF","No Eagle", "Eagle Regress", "Eagle Regress Back-half Season","Eagle Prop Est"),
  #                    
  #                    # labels = c("PF","No GSI", "GSI across days", "GSI by strata"),
  #                    values = c(15, 16, 17, 18,19))+
  # scale_fill_manual(values = wes_palette("IsleofDogs1"),
  #                   labels = c("PF (New; SS Recon.)",
  #                              "PF (Old; Eagle+ Harvest)",
  #                              "2.1 (New; SS Recon.)",
  #                              "2.1 (Old; Eagle + Harvest)"))+
  # labels = c("Ver 1", "Ver 1-PF", "Ver 2", "Ver2-PF"))+
  scale_x_continuous(breaks = c(testDays), labels = c(
    "June 2",
    "June 7",
    "June 12",
    "June 17",
    "June 22",
    "June 27",
    "July 2",
    "July 7",
    "July 12",
    "July 17",
    "July 22",
    "July 27",
    "Aug 1",
    "Aug 6",
    "Aug 11",
    "Aug 16",
    "Aug 21",
    "Aug 26",
    "Aug 31"
    ))+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18),
        axis.text = element_text(angle = 90))


# Percent Error Box Plots

PE_ver1 <- as.data.frame(RetroList_ver1RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver1$version <- "Ver 1"

PE_ver6 <- as.data.frame(RetroList_ver6RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver6$version <- "Ver 6"

PE_ver11 <- as.data.frame(RetroList_ver11$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver11$version <- "Ver 1.1"

PE_ver101 <- as.data.frame(RetroList_ver101$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver101$version <- "Ver 1.0.1"

PE_ver2 <- as.data.frame(RetroList_ver2$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2$version <- "Ver 2"

PE_ver21 <- as.data.frame(RetroList_ver2$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver21$version <- "Ver 2.1"

PE_ver2c <- as.data.frame(RetroList_ver2c$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c$version <- "Ver 2c"

PE_ver2cRR <- as.data.frame(RetroList_ver2cRR$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2cRR$version <- "Ver 2c RR"

PE_ver2c1 <- as.data.frame(RetroList_ver2c1$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c1$version <- "Ver 2c1"

PE_ver2c1_strata <- as.data.frame(RetroList_ver2c1_strata$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c1_strata$version <- "Ver 2c1 strata"

PE_ver2cSST <- as.data.frame(RetroList_ver2cSST$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2cSST$version <- "Ver 2cSST"

PE_ver2csd <- as.data.frame(RetroList_ver2csd$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2csd$version <- "Ver 2csd"

PE_ver2c01 <- as.data.frame(RetroList_ver2c01$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c01$version <- "Ver 2c01"

PE_ver2c02 <- as.data.frame(RetroList_ver2c02$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c02$version <- "Ver 2c02"

PE_ver2c02sst <- as.data.frame(RetroList_ver2c02sst$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c02sst$version <- "Ver 2c02 sst"

PE_ver2c01if <- as.data.frame(RetroList_ver2c01if$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c01if$version <- "Ver 2c01if"

PE_ver2c11 <- as.data.frame(RetroList_ver2c11$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c11$version <- "Ver 2c11"

PE_ver2_new_filled <- as.data.frame(RetroList_ver2_new_filled$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2_new_filled$version <- "Ver 2 new filled"

PE_ver2_old_filled <- as.data.frame(RetroList_ver2_old_filled$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2_old_filled$version <- "Ver 2 old filled"

PE_ver2_dyn <- as.data.frame(RetroList_ver2_dyn$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2_dyn$version <- "Ver 2 dyn"

# PE_ver2.1 <- as.data.frame(RetroList_ver2.1$PE_mat) %>% 
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver2.1$version <- "Ver 2.1"
# 
# PE_ver2.2 <- as.data.frame(RetroList_ver2.2$PE_mat) %>% 
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver2.2$version <- "Ver 2.2"

PE_ver201 <- as.data.frame(RetroList_ver201_old$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver201$version <- "Ver 201 old"

PE_ver201_new <- as.data.frame(RetroList_ver201_new$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver201_new$version <- "Ver 201 new"

PE_ver211 <- as.data.frame(RetroList_ver211_old$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver211$version <- "Ver 211 old"

PE_ver211_new <- as.data.frame(RetroList_ver211_new$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver211_new$version <- "Ver 211 new"

PE_ver2110 <- as.data.frame(RetroList_ver2110_old$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2110$version <- "Ver 2110 old"

PE_ver2110_new <- as.data.frame(RetroList_ver2110_new$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2110_new$version <- "Ver 2110 new"

PE_ver3.0 <- as.data.frame(RetroList_ver3.0$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver3.0$version <- "Ver 3.0"

PE_ver3.2 <- as.data.frame(RetroList_ver3.2$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver3.2$version <- "Ver 3.2"
# 
# PE_ver3.3 <- as.data.frame(RetroList_ver3.3$PE_mat) %>%
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver3.3$version <- "Ver 3.3"
# 
# PE_ver3.4 <- as.data.frame(RetroList_ver3.4$PE_mat) %>% 
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver3.4$version <- "Ver 3.4"
# 
# PE_ver3.5 <- as.data.frame(RetroList_ver3.5$PE_mat) %>% 
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver3.5$version <- "Ver 3.5"

# PE_ver4.0 <- as.data.frame(RetroList_ver4.0$PE_mat) %>%
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver4.0$version <- "Ver 4.0"

# PE_ver4.1 <- as.data.frame(RetroList_ver4.1$PE_mat) %>%
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver4.1$version <- "Ver 4.1"

# PE_ver4.2 <- as.data.frame(RetroList_ver4.2$PE_mat) %>%
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver4.2$version <- "Ver 4.2"

PE_ver5.0 <- as.data.frame(RetroList_ver5.0$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver5.0$version <- "Ver 5.0"

PE_pf_new <- as.data.frame(RetroList_pf_new$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_pf_new$version <- "PF New"

PE_pf_old <- as.data.frame(RetroList_pf_old$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_pf_old$version <- "PF"

peDF_total <- rbind(
                    PE_ver1,
                    PE_ver6
                    # PE_ver101,
                    # PE_ver2,
                    # PE_ver21,
                    # PE_ver2c,
                    # PE_ver2cRR
                    # PE_ver2c1,
                    # PE_ver2c1_strata
                    # PE_ver2cSST
                    # PE_ver2c01,
                    # PE_ver2c02,
                    # PE_ver2c02sst
                    # PE_ver2c01if
                    # PE_ver2c11
                    # PE_ver2csd
                    # PE_ver2_old_filled,
                    # PE_ver2_new,
                    # PE_ver2_dyn,
                    # PE_ver201,
                    # PE_ver201_new,
                    # PE_ver211,
                    # PE_ver2110,
                    # PE_ver2110_new,
                    # PE_ver211_new
                    # PE_ver2.1,
                    # PE_ver2.2,
                    # PE_ver3.0
                    # PE_ver3.2
                    # PE_ver3.3,
                    # PE_ver3.4,
                    # PE_ver3.5,
                    # PE_ver4.0,
                    # PE_ver4.1,
                    # PE_ver4.2,
                    # PE_ver5.0
                    # PE_pf_new,
                    # PE_pf_old
                    )

# For plotting PF point shape 

newDF <- left_join(x = pf_hist, y = CAN_hist_new)
# newDF <- left_join(x = newDF, y = CAN_hist_old, "Year")

# newDF$PF_new <- (newDF$Weighted_Forecast-newDF$can.mean.x)/newDF$can.mean.x

newDF$Year <- as.factor(newDF$Year)

newDF$PF_old <- (newDF$mean-newDF$can.mean)/newDF$can.mean

# long_df <- newDF %>% pivot_longer(cols = c(PF_new,PF_old)) %>% as.data.frame()

ggplot(peDF_total, aes(x = Year, y = PE*100, fill =version))+
  geom_boxplot(aes()) + 
  labs(fill = "")+
  geom_point(data = newDF,
             aes( x = Year, y = PF_old, fill = "PF"), 
             shape = 25,
             position = position_nudge(x = -.3),
             size = 3)+
  # scale_fill_colorblind()+
  # scale_fill_manual(values = wes_palette("IsleofDogs1"),
  #                   labels = c(
  #                              "2.1 (New SS Recon.)",
  #                              "2.1 (Old Eagle + Harvest)"))+

  # scale_colour_manual(name = "",
  #                     labels = c("No GSI", "GSI across days", "GSI by strata"),
  #                     values = c("blue", "green", "gold","purple"))+
  # scale_fill_colorblind(name = "",
  #                     labels = c("No GSI", "GSI across days", "GSI by strata"))+
  # 
  geom_hline(yintercept  = 0, linetype = 2, color = "orange", size = 1)+
  # theme(legend.position = "top",
  #       panel.grid.major.y  =element_line(linetype = 2, size = .5, color = "white"),
  #       panel.grid.major.x = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       panel.border = element_blank(),
  #       plot.background = element_blank(),
  #       panel.background = element_rect(fill = "grey"),
  #       text = element_text(size = 18))+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 12),
        axis.text = element_text(angle = 90))+
  facet_wrap(~Year, scales = "free")
  # ggtitle("New Can_hist ver 2")


# Turn off PDF save
# dev.off()

# Extract divergent transition warnings
day.vect  <- vector(length = length(outputlist_ver3))
year.vect <- vector(length = length(outputlist_ver3))
div_trans <- vector(length = length(outputlist_ver3))

for (i in 1:length(outputlist_ver3)) {
  
  # i = 1
  
  day.vect[i] <- outputlist_ver3[[i]]$myDay
  year.vect[i] <- outputlist_ver3[[i]]$myYear
  div_trans[i] <- outputlist_ver3[[i]]$divergent_trans
  
}

div.df <- data.frame("day" = day.vect,
                     "year" = year.vect,
                     "div_trans" = div_trans)

div.df
#


# Plot of ~ model weights of PF and PSS over a year #################################
sigma_func <- function(mod, testYears, testDays){
  
  # Vector for storing sigmas
  sigma_vect <- vector(length=length(mod))
  
  # Loop for geting vector sigmas 
  for (p in 1:length(sigma_vect)) {
    sigma_vect[p]<- median(mod[[p]]$pars$sigma)
  } # End loop
  
  # Matrix of PSS prediction
  sigma_mat <- matrix(sigma_vect,
                      nrow = length(testDays),
                      ncol = length(testYears),
                      byrow = FALSE)
  # Label with years
  colnames(sigma_mat)<- c(testYears)
  
  # Turn matrix into DF
  sigma_DF <- as.data.frame((sigma_mat))
  
  # Put days into DF
  sigma_DF$day <- testDays
  
  # Pivot df longer
  sigma_DF <- as.data.frame(pivot_longer(sigma_DF,cols = -day))
  
  # Name columns
  names(sigma_DF)<- c("Day", "Year", "PSSpred")
  
  # Change year to double
  sigma_DF$Year<- as.double(sigma_DF$Year)
  
  # Return df
  return(sigma_DF)
  
} # End Function


sig_mat<-sigma_func(mod = mod,
                    testYears = testYears,
                    testDays = testDays)

sig_2019<- sig_mat[sig_mat$Year == 2019,]

pf_sigma

weight_pf <- 1/pf_sigma^2

weight_PSS <- 1/sig_2019$PSSpred^2

sig_DF <- data.frame("PSS" = weight_PSS, "PF" = weight_pf)

sig_DF$PSS_stand <- sig_DF$PSS/(sig_DF$PSS+sig_DF$PF)
sig_DF$PF_stand <- sig_DF$PF/(sig_DF$PSS+sig_DF$PF)
sig_DF$Day <- testDays
sig_DF <- sig_DF[,c(3,4,5)]

sig_long <- pivot_longer(sig_DF, cols = -Day)

ggplot(sig_long, aes(x = Day, y = value, fill = name))+
  geom_area()





##### Density Plots with plot function #####################################################
gg <- outPlots(outputList = model.output, 
         CAN_hist = CAN_hist_old,
         GSI = FALSE, 
         Retrospective = TRUE,
         eagle = TRUE) 
# gg +  geom_point(data = data.frame("x" = sum(Eagle_hist$count[Eagle_hist$Day<=208 &
#                                                       Eagle_hist$Year == 2007]),
#                                "y" = CAN_hist_old$can.mean[CAN_hist_old$Year == 2007]),
#              aes( x = x/1000, y = y/1000), col = "red", size = 20)

# List to hold plots
plots_list <- list()
for (i in 1:length(ver102_RR_T)){
  
  
  plots_list[[i]] <- density.func(outputList = ver102_RR_T[[i]], CAN_hist = CAN_hist_new)
}

# Vector of names form plot list
test<-names(ver102_RR_T)
names(plots_list)<-test



figure<-ggarrange(plots_old_can_ver3.5$`2021_158`$DensPlot,
                  plots_old_can_ver3.5$`2021_178`$DensPlot,
                  plots_old_can_ver3.5$`2021_198`$DensPlot,
                  plots_old_can_ver3.5$`2019_153`$DensPlot,
                  plots_old_can_ver3.5$`2019_178`$DensPlot,
                  plots_old_can_ver3.5$`2019_198`$DensPlot,
                  plots_old_can_ver3.5$`2017_153`$DensPlot,
                  plots_old_can_ver3.5$`2017_178`$DensPlot,
                  plots_old_can_ver3.5$`2017_198`$DensPlot,
                  ncol = 3, 
                  nrow = 3, 
                  common.legend = TRUE,
                  label.x = .5,
                  label.y = 1,
                  labels = c("June 2 2021 ","June 27 2021","July 17 2021",
                             "June 2 2019","June 27 2019","July 17 2019",
                             "June 2 2017", "June 27 2017", "July 17 2017"),
                  font.label = list(size = 10, 
                                    color = "black",
                                    family = "serif"))

figure<-ggarrange(plots_list$`2022_158` + labs(x = "", y ="")+
                    theme_light(),
                  plots_list$`2022_168`$PredPlot + labs(x = "", y ="")+
                    theme_light(),
                  plots_list$`2022_183`$PredPlot + labs(x = "", y ="")+
                    theme_light(),
                  plots_list$`2022_213`$PredPlot + labs(x = "", y ="")+
                    theme_light(),
                  ncol = 4, 
                  nrow = 1, 
                  common.legend = TRUE,
                    
                  label.x = .4,
                  label.y = 1,
                  labels = c("June 7 2022",
                             "June 17 2022",
                             "July 2 2022",
                             "Aug 1 2022"),
                  font.label = list(size = 18, 
                                    color = "black",
                                    family = "serif"))

annotate_figure(figure,
                left = textGrob("Relative Probability", 
                                        rot = 90,
                                        vjust = 1, 
                                gp = gpar(cex = 2)),
                bottom = textGrob("Thousands of Chinook Salmon", 
                                  gp = gpar(cex = 2)),
                fig.lab.size = 50 )







figure<-ggarrange(plots_old_can_ver4.0$`2021_178`$EaglePlot,
                  plots_old_can_ver4.0$`2021_193`$EaglePlot,
                  plots_old_can_ver4.0$`2021_213`$EaglePlot,
                  ncol = 3, 
                  nrow = 1, 
                  common.legend = TRUE,
                  
                  label.x = .35,
                  label.y = 1,
                  labels = c("June 27 2021 ","July 12 2021","August 1 2021"),
                  font.label = list(size = 18, 
                                    color = "black",
                                    family = "serif"))

annotate_figure(figure,
                left = textGrob("Run Size (Eagle+Harvest)", 
                                rot = 90,
                                vjust = 1, 
                                gp = gpar(cex = 2)),
                bottom = textGrob("Eagle Chinook Passage (1000's)", 
                                  gp = gpar(cex = 2)),
                fig.lab.size = 50 )



plots_list$`2020_153`$PredPlot+
  xlab("")

