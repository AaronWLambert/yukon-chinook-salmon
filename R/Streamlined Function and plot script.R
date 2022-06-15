#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.14.22
# Version: This will run the selected model version in the function
# Purpose: This is a working script to generate all outputs and figures (at this time)
# 
#
#   1) Read in data
#   2) Call the function
#   3) Generate outputs as "pars" and generated figures from the posterior outputs
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
require(rstan)
require(rstanarm)
require(bayesplot)
require(tidyverse)
require(ggthemes)
require(viridis)
require(shinystan)
require(lubridate)
require(ggpubr)
require(gridExtra)
require(tidybayes)
library(wesanderson)
require(grid)


# Parralize for optimum model run time
rstan_options(auto_write = TRUE)
#
mc.cores = parallel::detectCores()
# mc.cores <-1

# Define Workflow Paths ###########################################################
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

# Call in functions from model_run_funciton.r
source(file = file.path(dir.R,"model_run_function.r"))

######### Import Data ##############################################################
# Historical Canadian EOS reconstructed run
# This is the reconstructed data from Curry for old reconstructed modeling procedure
CAN_hist <- readRDS(file.path(dir.data,"Can Origin Reconstructed 2Feb22.RDS"))
# PSS historical for Jun 1st = 1, to Aug 31 = 92
# Years 1995 - 2020
# PSScum_hist <- readRDS(file.path(dir.data,"/clean_counts.rds"))

# Read in PSS non-cummulative with date in date format
# PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean <- readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
pf_hist <- readRDS(file.path(dir.data,"inv_var_weighted_forcast_v3_June072022.RDS"))
# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 27April22.RDS"))


# Control Section ##################################################################
model.version <- "3.2"

# Range of years 1995 to 2022 (1996 is not available)
myYear <- 2022

# Range of days 152 (June 1) - 243 (August 31)
myDay <- 165

# MCMC Parameters
n.chains <- 4
n.iter <- 30000;#5e4
n.thin <- 2


# Run the model for a single year and day of interest.
#  Ajust the inputs for the model in the control section above
model.output <-InSeasonProjection(model.version = model.version,
                                  myYear = myYear,
                                  myDay = myDay,
                                  n.chains = n.chains,
                                  n.iter = n.iter,
                                  n.thin = n.thin,
                                  CAN_hist = CAN_hist,
                                  PSS_hist = PSS_hist,
                                  # GSI_mean = GSI_mean,
                                  pf_hist = pf_hist,
                                  GSI_by_year = GSI_by_year,
                                  savefit = TRUE)


##### Density Plots with plot function #############################################

# Use the outPlots function to generate two different kinds of plots
#     1. Density plot with preseason forecast draws, PSS prediction draws,
#        and the projected Canadian-origin abundance draws.
#     2. Plot of PSS passage against EOS reconstructed Canadian abundance with
#        High Density Intervals (HDI) of 50 and 95 %.

# Select GSI = TRUE if running model versions 3
# Select Retrospective = FALSE if running inseason projection for current year
outPlots(outputList = model.output ,
         CAN_hist = CAN_hist, GSI = TRUE, Retrospective = FALSE)

# Generate intervals *Note that pars may be changed to any parameter of interest
mcmc_intervals(x = model.output$fit, pars = "RunSize", prob_outer = .9, prob = .5 )


# Scratch work section #########################################################
quant <-quantile(x = model.output$pars$RunSize, probs = c(.1,.9))

gg <- outPlots(outputList = model.output ,
         CAN_hist = CAN_hist,GSI = FALSE, Retrospective = FALSE)

df1<-data.frame("value" = model.output$pars$RunSize, "par" = "full")

pp <- ggplot(data = df1, aes(x = value))+
  geom_density(fill = "red", alpha = .5)

d <- ggplot_build(pp)$data[[1]]

pp+
  geom_area(data = subset(d, x > quant[[1]] & x< quant[[2]]), 
                          aes(x = x, y = y), alpha = .5, fill = "green")

xx<-model.output$pars$RunSize

ggplot(as.data.frame(xx), aes(x = xx/1000))+
  stat_ecdf()
