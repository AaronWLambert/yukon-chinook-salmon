# Aaron Lambert
# 10/4/22
# Plots to conduct retrospective testing


library(wesanderson)
library(tidyverse)
library(ggthemes)
# library(plyr)

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
# CAN_hist <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
# This is the reconstructed data from Curry for old reconstructed modeling procedure
CAN_hist_old <- readRDS(file.path(dir.data,"Can Origin Reconstructed 2Feb22.RDS"))
CAN_hist_new <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
# PSS historical for Jun 1st = 1, to Aug 31 = 92
# Years 1995 - 2020
# PSScum_hist <- readRDS(file.path(dir.data,"/clean_counts.rds"))

# Read in PSS non-cummulative with date in date format
# PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage 22Jun22.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean <- readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
pf_hist <- readRDS(file.path(dir.data,"temp_PF_insamp_12Dec22.RDS"))
# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 27April22.RDS"))



# Year included for retro testing
testYears <- c(2007:2022)
# Years included for retro testing with new reconstructed canadian total
# Only avialable (currently) for years up to 2019
# testYears_short_new <- c(2007:2010,2013:2019)
# # Test Years for regression including only 2005 and onward
# testYears_short <- c(2009:2010,2013:2021)
# Test days for a full year of testing
testDays <- seq(from = 153, to = 243, by = 5)
# Test days for abbreviated cumulative counts
# Used to look at if zeros are true zeros
# testDays_short <- seq(from = 168, to = 213, by = 5)

# Read in data list 
# mod_old <- readRDS(file = file.path(dir.output,"OutPut_ver2_oldCan_2007_2010_2013_2021_5apr22.RDS"))
# mod_new <- readRDS(file = file.path(dir.output, "OutPut_ver2_NEWCan_2007_2010_2013_2019_12oct22.RDS"))
# mod_new_2.1 <- readRDS(file = file.path(dir.output, "OutPut_ver21_NEWCan_2007_2010_2013_2019_12oct22.RDS"))
# mod_new_filled <- readRDS(file = file.path(dir.output,"OutPut_ver2_NEWCan_frontfilled_2007_2010_2013_2019_12oct22.RDS"))
# mod_old_filled <- readRDS(file = file.path(dir.output,"OutPut_ver2_oldCan_frontfilled_2007_2010_2013_2019_12oct22.RDS"))
# mod_short_old <- readRDS(file = file.path(dir.output,"OutPut_ver2.short_oldCan_2007_2010_2013_2021_12oct22.RDS"))
# mod_dyn <- readRDS(file = file.path(dir.output,"OutPut_ver2_oldCan_dynamicbeta_2007_2010_2013_2019_12oct22.RDS"))
# mod_old_201 <- readRDS(file = file.path(dir.output,"OutPut_ver201_oldCan_2007_2010_2013_2019_27oct22.RDS"))
# mod_new_201 <- readRDS(file = file.path(dir.output,"OutPut_ver201_NEWCan_2007_2010_2013_2019_27oct22.RDS"))
# mod_new_211 <- readRDS(file = file.path(dir.output,"OutPut_ver211_NEWCan_2007_2010_2013_2019_4Nov22.RDS"))
# mod_old_211 <- readRDS(file = file.path(dir.output,"OutPut_ver211_oldCan_2007_2010_2013_2021_3Nov22.RDS"))
# mod_old_2110 <- readRDS(file = file.path(dir.output,"OutPut_ver2110_oldCan_2007_2010_2013_2021_28Nov22.RDS"))
# mod_old_32 <- readRDS(file = file.path(dir.output,"OutPut_ver32_oldCan_2007_2010_2013_2021_22Nov22.RDS"))

mod_2c_long <- readRDS(file = file.path(dir.output,"Ver2c_long_28Feb23.RDS"))

mod_2c01_long <- readRDS(file = file.path(dir.output,"Ver201c_long_28Feb23.RDS"))

mod_2csd_long <- readRDS(file = file.path(dir.output,"Ver2csd_long_10Mar23.RDS"))


# pdf(file = file.path(dir.output,"sigma_retro_plots2newvsold_13Oct22.pdf"))

# Functions for extracting parameter values ################################################

runsize_trend_func <- function(mod, testYears, testDays, CAN_hist, PF){
  
# Matrix of Runsize estimates
RunSize_vect <- vector(length=length(mod))

for (p in 1:length(RunSize_vect)) {
  RunSize_vect[p]<- median(mod[[p]]$pars$RunSize)
}

RunSize_mat <-matrix(RunSize_vect, 
                     nrow = length(testDays), 
                     ncol = length(testYears),
                     byrow = FALSE)
colnames(RunSize_mat)<- c(testYears)
RunSize_DF <- as.data.frame((RunSize_mat))
# RunSize_DF$day <- testDays
RunSize_DF$day <- testDays
RunSize_DF <- as.data.frame(pivot_longer(RunSize_DF,cols = -day))
names(RunSize_DF)<- c("Day", "Year", "RunSize")
RunSize_DF$Year<- as.double(RunSize_DF$Year)
RunSize_DF <- left_join(RunSize_DF,CAN_hist)
RunSize_DF <- left_join(RunSize_DF, pf_hist)

# Matrix of 2007 Runsize estimates
PSSpred_vect <- vector(length=length(mod))

for (p in 1:length(RunSize_vect)) {
  PSSpred_vect[p]<- median(mod[[p]]$pars$post_curr_predPSS)
}

# Matrix of PSS prediction
PSSpred_mat <-matrix(PSSpred_vect,
                     nrow = length(testDays),
                     ncol = length(testYears),
                     byrow = FALSE)
# PSSpred_mat <-matrix(PSSpred_vect, 
#                      nrow = length(testDays_short), 
#                      ncol = length(testYears),
#                      byrow = FALSE) 
colnames(PSSpred_mat)<- c(testYears)
PSSpred_DF <- as.data.frame((PSSpred_mat))
# PSSpred_DF$day <- testDays
PSSpred_DF$day <- testDays
PSSpred_DF <- as.data.frame(pivot_longer(PSSpred_DF,cols = -day))
names(PSSpred_DF)<- c("Day", "Year", "PSSpred")
PSSpred_DF$Year<- as.double(PSSpred_DF$Year)
RunSize_DF <- left_join(PSSpred_DF,RunSize_DF)
RunSize_DF$Version <- as.factor(mod[[1]]$version)
return(RunSize_DF)
} # End of function

# shortVer2 <- runsize_trend_func(mod = mod, testYears = testYears_short, testDays = testDays)
# ver2_old <- runsize_trend_func(mod = mod_old,
#                                testYears = testYears,
#                                testDays = testDays ,
#                                CAN_hist = CAN_hist_old)
# ver2_new <- runsize_trend_func(mod = mod_new,
#                                testYears = testYears_short_new,
#                                testDays = testDays,
#                                CAN_hist = CAN_hist_new)
# ver2.1_new <- runsize_trend_func(mod = mod_new_2.1,
#                                testYears = testYears_short_new,
#                                testDays = testDays,
#                                CAN_hist = CAN_hist_new)
# ver2_new_filled <- runsize_trend_func(mod = mod_new_filled, 
#                                testYears = testYears_short_new,
#                                testDays = testDays,
#                                CAN_hist = CAN_hist_new)
# ver2_old_filled <- runsize_trend_func(mod = mod_old_filled, 
#                                       testYears = testYears,
#                                       testDays = testDays,
#                                       CAN_hist = CAN_hist_old)
# ver2_dyn <- runsize_trend_func(mod = mod_dyn,
#                                  testYears = testYears,
#                                  testDays = testDays,
#                                CAN_hist = CAN_hist_old)
# ver2_old_short <- runsize_trend_func(mod = mod_short_old ,
#                                testYears = testYears_short,
#                                testDays = testDays ,
#                                CAN_hist = CAN_hist_old)
# ver201_old <- runsize_trend_func(mod = mod_old_201,
#                                      testYears = testYears,
#                                      testDays = testDays ,
#                                      CAN_hist = CAN_hist_old)
ver2c_old <- runsize_trend_func(outputlist_ver2c,
                                  testYears = testYears,
                                  testDays = testDays ,
                                  CAN_hist = CAN_hist_old)

ver2cSST_old <- runsize_trend_func(outputlist_ver2cSST,
                                testYears = testYears,
                                testDays = testDays ,
                                CAN_hist = CAN_hist_old)

ver2csd_old <- runsize_trend_func(mod = mod_2csd_long,
                                testYears = testYears,
                                testDays = testDays ,
                                CAN_hist = CAN_hist_old)

ver2c01_old <- runsize_trend_func(mod = outputlist_ver2c01,
                                 testYears = testYears,
                                 testDays = testDays ,
                                 CAN_hist = CAN_hist_old)

ver2c02_old <- runsize_trend_func(mod = outputlist_ver2c02,
                                  testYears = testYears,
                                  testDays = testDays ,
                                  CAN_hist = CAN_hist_old)

ver2c02sst_old <- runsize_trend_func(mod = outputlist_ver2c02sst,
                                  testYears = testYears,
                                  testDays = testDays ,
                                  CAN_hist = CAN_hist_old)

ver2c01if_old <- runsize_trend_func(mod = outputlist_ver2c01if,
                                  testYears = testYears,
                                  testDays = testDays ,
                                  CAN_hist = CAN_hist_old)

ver5 <- runsize_trend_func(mod = outputlist_ver5,
                           testYears = testYears,
                           testDays = testDays,
                           CAN_hist = CAN_hist_old)
# ver201_new <- runsize_trend_func(mod = mod_new_201,
#                                   testYears = testYears_short_new,
#                                   testDays = testDays ,
#                                   CAN_hist = CAN_hist_new)
# ver211_new <- runsize_trend_func(mod = mod_new_211,
#                                  testYears = testYears_short_new,
#                                  testDays = testDays ,
#                                  CAN_hist = CAN_hist_new)
# 
# ver211_old <- runsize_trend_func(mod = mod_old_211,
#                                  testYears = testYears,
#                                  testDays = testDays ,
#                                  CAN_hist = CAN_hist_old)
# 
# ver2110_old <- runsize_trend_func(mod = mod_old_2110,
#                                  testYears = testYears,
#                                  testDays = testDays ,
#                                  CAN_hist = CAN_hist_old)
# 
# ver32_old <- runsize_trend_func(mod = mod_old_32,
#                                   testYears = testYears,
#                                   testDays = testDays ,
#                                   CAN_hist = CAN_hist_old)
# 



# ver2_old$Model <- "Old Can"
# ver2_old$Model <- as.factor(ver2_old$Model)
# ver2_old$CanMethod <- "Old Can"
# ver2_old$CanMethod <- as.factor(ver2_old$CanMethod)
# 
# ver2_old_short$Model <- "Old Can Short"
# ver2_old_short$Model <- as.factor(ver2_old_short$Model)
# ver2_old_short$CanMethod <- "Old Can"
# ver2_old_short$CanMethod <- as.factor(ver2_old_short$CanMethod)
# 
# ver2_new$Model <- "New Can"
# ver2_new$Model <- as.factor(ver2_new$Model)
# ver2_new <- ver2_new[,c(1:5,7,8)]
# ver2_new$CanMethod <- 'New Can'
# ver2_new$CanMethod <- as.factor(ver2_new$CanMethod)
# 
# ver2.1_new$Model <- "New Can 2.1"
# ver2.1_new$Model <- as.factor(ver2.1_new$Model)
# ver2.1_new <- ver2.1_new[,c(1:5,7,8)]
# ver2.1_new$CanMethod <- 'New Can'
# ver2.1_new$CanMethod <- as.factor(ver2.1_new$CanMethod)
# 
# ver2_new_filled$Model <- "New Can Filled"
# ver2_new_filled$Model <- as.factor(ver2_new_filled$Model)
# ver2_new_filled <- ver2_new_filled[,c(1:5,7,8)]
# ver2_new_filled$CanMethod <- 'New Can'
# ver2_new_filled$CanMethod <- as.factor(ver2_new_filled$CanMethod)
# 
# ver2_old_filled$Model <- "Old Can Filled"
# ver2_old_filled$Model <- as.factor(ver2_old_filled$Model)
# # ver2_old_filled <- ver2_old_filled[,c(1:5,7,8)]
# ver2_old_filled$CanMethod <- 'Old Can'
# ver2_old_filled$CanMethod <- as.factor(ver2_old_filled$CanMethod)
# 
# ver2_dyn$Model <- "Dynamic"
# ver2_dyn$Model <- as.factor(ver2_dyn$Model)
# ver2_dyn$CanMethod <- "Old Can"
# ver2_dyn$CanMethod <- as.factor(ver2_dyn$CanMethod)
# 
# ver201_old$Model <- "Old Can 201"
# ver201_old$Model <- as.factor(ver201_old$Model)
# ver201_old$CanMethod <- "Old Can"
# ver201_old$CanMethod <- as.factor(ver201_old$CanMethod)

ver2c_old$Model <- "Old Can 2c"
ver2c_old$Model <- as.factor(ver2c_old$Model)
# ver2c_old$CanMethod <- "Old Can"
# ver2c_old$CanMethod <- as.factor(ver2c_old$CanMethod)

ver2cSST_old$Model <- "Old Can 2cSST"
ver2cSST_old$Model <- as.factor(ver2cSST_old$Model)
# ver2cSST_old$CanMethod <- "Old Can"
# ver2cSST_old$CanMethod <- as.factor(ver2cSST_old$CanMethod)

ver2csd_old$Model <- "Old Can 2csd"
ver2csd_old$Model <- as.factor(ver2csd_old$Model)
ver2csd_old$CanMethod <- "Old Can"
ver2csd_old$CanMethod <- as.factor(ver2csd_old$CanMethod)

ver2c01_old$Model <- "Old Can 2c01"
ver2c01_old$Model <- as.factor(ver2c01_old$Model)
# ver2c01_old$CanMethod <- "Old Can"
# ver2c01_old$CanMethod <- as.factor(ver2c01_old$CanMethod)

ver2c02_old$Model <- "Old Can 2c02"
ver2c02_old$Model <- as.factor(ver2c02_old$Model)
# ver2c02_old$CanMethod <- "Old Can"
# ver2c02_old$CanMethod <- as.factor(ver2c02_old$CanMethod)

ver2c02sst_old$Model <- "Old Can 2c02 sst"
ver2c02sst_old$Model <- as.factor(ver2c02sst_old$Model)
# ver2c02_old$CanMethod <- "Old Can"
# ver2c02_old$CanMethod <- as.factor(ver2c02_old$CanMethod)

ver2c01if_old$Model <- "Old Can 2c01if"
ver2c01if_old$Model <- as.factor(ver2c01if_old$Model)
# ver2c01if_old$CanMethod <- "Old Can"
# ver2c01if_old$CanMethod <- as.factor(ver2c01if_old$CanMethod)

# ver201_new$Model <- "New Can 201"
# ver201_new$Model <- as.factor(ver201_new$Model)
# ver201_new <- ver201_new[,c(1:5,7,8)]
# ver201_new$CanMethod <- 'New Can'
# ver201_new$CanMethod <- as.factor(ver201_new$CanMethod)
# 
# ver211_new$Model <- "New Can 211"
# ver211_new$Model <- as.factor(ver211_new$Model)
# ver211_new <- ver211_new[,c(1:5,7,8)]
# ver211_new$CanMethod <- 'New Can'
# ver211_new$CanMethod <- as.factor(ver211_new$CanMethod)
# 
# ver211_old$Model <- "Old Can 211"
# ver211_old$Model <- as.factor(ver211_old$Model)
# ver211_old$CanMethod <- "Old Can"
# ver211_old$CanMethod <- as.factor(ver211_old$CanMethod)
# 
# ver2110_old$Model <- "Old Can 2110"
# ver2110_old$Model <- as.factor(ver2110_old$Model)
# ver2110_old$CanMethod <- "Old Can"
# ver2110_old$CanMethod <- as.factor(ver2110_old$CanMethod)
# 
# ver32_old$Model <- "Old Can 32"
# ver32_old$Model <- as.factor(ver32_old$Model)
# ver32_old$CanMethod <- "Old Can"
# ver32_old$CanMethod <- as.factor(ver32_old$CanMethod)

ver5$Model <- "Old Can 5.0"
ver5$Model <- as.factor(ver5$Model)
ver5$CanMethod <- "Old Can"
ver5$CanMethod <- as.factor(ver5$CanMethod)

comb.df <- rbind(
                 # ver2_new,
                 # ver2_old,
                 # ver2_dyn,
                 # ver2_new_filled,
                 # ver2_old_filled
                 # ver2_old_short,
                 # ver2.1_new
                 # ver201_new,
                 ver2c_old,
                 # ver2cSST_old
                 ver2c01_old,
                 ver2c02_old,
                 # ver2c02sst_old
                 ver2c01if_old
                 # ver2csd_old
                 # ver211_new,
                 # ver211_old,
                 # ver2110_old
                 # ver32_old
                 # ver5
                 )

ggplot(comb.df, aes(x = Day, 
                    y = RunSize/1000,
                    color = Model))+
  geom_line(size = 2)+
  # geom_line( aes(x = Day, y = PSSpred/1000, color = "PSS Pred" ))+
  # geom_line(data = PSSpred_DF_mult, aes(x = Day, y = PSSpred/1000, color = "PSS Pred 2x"))+
  # geom_line(data = RunSize_DF_mult, aes(x = Day, y = RunSize/1000, color = "RunSize 2x"))+
  geom_hline(aes(yintercept = can.mean/1000, 
                 # linetype = CanMethod,
                 ), size = 1, color="purple")+
  geom_hline(aes(yintercept= mean/1000, color = "PF"),size = 1)+
  # geom_hline(data = CAN_hist_new[CAN_hist_new$Year >= 2007,],aes(yintercept = can.mean/1000))+
  facet_wrap(~Year, scales = "free_y", nrow = 3)+
  theme(legend.position = "top")+
  ggtitle(paste("Model Version 2.0"))+
  scale_color_colorblind()

ggplot(comb.df, aes(x = Day, y = PSSpred/1000, color = Model))+
  geom_line(size = 2)+
  # geom_line( aes(x = Day, y = PSSpred/1000, color = "PSS Pred" ))+
  # geom_line(data = PSSpred_DF_mult, aes(x = Day, y = PSSpred/1000, color = "PSS Pred 2x"))+
  # geom_line(data = RunSize_DF_mult, aes(x = Day, y = RunSize/1000, color = "RunSize 2x"))+
  geom_hline(aes(yintercept = can.mean/1000, 
                 color = "EOS Can-origin Chinook"
  ), size = 1.5,linetype = 2)+
  geom_hline(aes(yintercept= mean/1000, color = "PF"))+
  facet_wrap(~Year, scales = "free_y", nrow = 3)+
  theme(legend.position = "top")+
  ggtitle(paste("Model Version 2.0"))


ggplot(comb.df, aes(x = Day, y = RunSize/1000, color = Model))+
  geom_line(size = 2)+
  geom_line( aes(x = Day, y = PSSpred/1000, color = Model ), size = .75, linetype = 1)+
  # geom_line(data = PSSpred_DF_mult, aes(x = Day, y = PSSpred/1000, color = "PSS Pred 2x"))+
  # geom_line(data = RunSize_DF_mult, aes(x = Day, y = RunSize/1000, color = "RunSize 2x"))+
  geom_hline(aes(yintercept = can.mean/1000, 
                 color = "EOS Can-origin Chinook"
  ), size = 1.5,linetype = 2)+
  geom_hline(aes(yintercept= mean/1000, color = "PF"))+
  facet_wrap(~Year, scales = "free_y")+
  theme(legend.position = "top")+
  ggtitle(paste("Model Version 2.0"))

# Plot the error
comb.df$error <- comb.df$RunSize - comb.df$can.mean

ggplot(comb.df, aes(x = Day, 
                    y = error/1000,
                    color = Model
                    
                    ))+
  geom_line(size = 1.5, show.legend = T)+
  # geom_line( aes(x = Day, y = PSSpred/1000, color = "PSS Pred" ))+
  # geom_line(data = PSSpred_DF_mult, aes(x = Day, y = PSSpred/1000, color = "PSS Pred 2x"))+
  # geom_line(data = RunSize_DF_mult, aes(x = Day, y = RunSize/1000, color = "RunSize 2x"))+
  geom_hline(aes(yintercept = 0 , col = "Truth"), size = 1, linetype = 2)+
  facet_wrap(~Year ,
             nrow = 4,
             # scales = "free_y"
             )+
  theme(legend.position = "top")+
  labs(color = "", alpha = "")+
  ylab("Predicted - Observed (1000's Chinook)")+
  # scale_color_colorblind(name = "",
  #                        labels = c(
  #                                   "PSS Only ",
  #                                   "PSS and Eagle",
  #                                   "EOS Run Size"))+
  scale_x_continuous(breaks = c(153,173,193,213,232), labels = c(
    "June 2",
    
    # "June 12",
    
    "June 22",
    
    # "July 2",
    
    "July 12",
    
    # "July 22",
    
    "Aug 1",
    
    "Aug 20"))+
  # scale_color_colorblind(
  #                   labels = c(
  #                     "2.c",
  #                       "2.c.0.1 (w/ Eagle)",
  #                     "2.c.SD",
  #                     "Proverbial Truth"
  #                     ))+
  scale_colour_manual(name = "",
                      # labels = c("PF","Ver 2.c", "Ver 2.c.0.1", "Ver 2.c.sd"),
                      # labels = c("PF","No GSI", "GSI across days", "GSI by strata"),
                      labels = c("No Eagle", "Eagle Regress", "Eagle Regress Back-half Season","Eagle Prop Est","Truth"),
                      
                      values = c("blue", "green", "gold", "brown","black")) +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        
        legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90,
                                   size = 12))
  # guides(alpha = "none")

# Parameter Interval Analysis #######################################################
# Function to get intervals of parameters
interval_function <- function(mod, parameter = "", testDays, year = "", counter){
  # mod = list of model outputs from retrospective testing
  # parameter = parameter to get values for
  # testDays = days used for retrospective testing
  # year = year that parameter values are extracted for      
     xx.mat <- matrix(nrow = length(mod[[1]]$pars[[parameter]]),
                      ncol = length(testDays), byrow = F) 
  for (d in 1:length(testDays)) {
    for(v in 1:length(mod[[1]]$pars[[parameter]]))
      
      xx.mat[v,d]<- mod[[d+(counter)]]$pars[[parameter]][[v]]
    
  } #end of loop
  
  quant.xx <- apply(X = xx.mat, 
                    MARGIN = 2, 
                    FUN = quantile, 
                    probs=c(0.025, 0.25, 0.5, 0.75, 0.975),)
  
  pred.xx.df <- data.frame(testDays,
                           t(quant.xx))
  
  names(pred.xx.df) <- c("day","low95","low50","median",
                         "up50","up95")
  
  pred.xx.df$Year <- year
  
  # pred.xx.df$Year <- as.factor(pred.xx.df$Year)
  
  return(pred.xx.df)
}# End of function

###################################################

# Sigma from regression fitting
# Loop through years to calculate intervals for each year and day
years <- c(testYears)
n.years <- length(years)
count <- 0
for (y in 1:n.years) {
  # y = 2
  year <- years[y]
  df <-interval_function(mod = outputlist_ver1RR_T, 
                         parameter = "sigma", # Parameter of interest
                         testDays = testDays,
                         year = year,
                         counter = count)
  if(y == 1){
  final.df <- df}else{
    
    final.df <- rbind(final.df,df)
    
  }
  count <- count+ length(testDays)
}

sigmaDF <- final.df
sigmaDF$param <- "regression"



# Extract empirical sigma 
# Loop through years to calculate intervals for each year and day

count <- 0
for (y in 1:n.years) {
  # y = 2
  year <- years[y]
  df <-interval_function(mod = outputlist_ver1RR_T, 
                         parameter = "sigma_predPSS", # Parameter of interest
                         testDays = testDays,
                         year = year,
                         counter = count)
  if(y == 1){
    final.df <- df}else{
      
      final.df <- rbind(final.df,df)
      
    }
  count <- count+ length(testDays)
}

# Add param catagory
final.df$param <- "empirical"

# Call new name
sigma.predDF <- final.df

# Bind the df
sigma_plotDF <- rbind(sigma.predDF,sigmaDF)
















# Plot of RunSize with intervals
ggplot(new.df[new.df$Year==2017,], aes(x = day, by = param))+
  # geom_point() +
  # geom_ribbon(aes(ymin=low95, ymax=up95,fill = model), alpha=0.25) +
  # geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.2) +
  geom_ribbon(aes(ymin=low50/1000, 
                  ymax=up50/1000, 
                  fill = "HDI 50"), 
              alpha=0.5) +
  geom_line(aes(y=median/1000, 
                color = "Median"),
            size = 2, 
            show.legend = T) +
  labs(x = "Day of Year",
       y = "Thousands of Chinook Salmon")+
  geom_hline(aes(yintercept = can.mean/1000,
                 color = "Realized EOS Abundance"),
             size = 1.5,
             linetype = 2)+
  geom_hline(data = pf_hist, aes(yintercept = mean/1000, col = "PF"), 
             size = 1.5,
             linetype = 3)+
  # facet_wrap(~Year,nrow = 3)+
  scale_fill_calc()+
  scale_color_colorblind()+
  scale_x_continuous(breaks = c(153,173,193,213,232), labels = c(
    "June 2",
    
    # "June 12",
    
    "June 22",
    
    # "July 2",
    
    "July 12",
    
    # "July 22",
    
    "Aug 1",
    
    "Aug 20"
    ))+
  # ggtitle(paste("Model Version", mod$`2007_153`$version,"\n Parameter = sigma normal curve"))+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        text = element_text(size = 18))


# # Plot of RunSize with intervals
# ggplot(new.df, aes(x = day))+
#   # geom_point() +
#   # geom_ribbon(aes(ymin=low95, ymax=up95,fill = model), alpha=0.25) +
#   geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "pink"),
#               
#               alpha=0.7,
#               show.legend = T) +
#   geom_ribbon(aes(ymin=low50/1000, 
#                   ymax=up50/1000, 
#                   fill = "red"), 
#               alpha=0.3) +
#   geom_line(aes(y=median/1000, 
#                 col = "darkred"),
#             size = 2, 
#             show.legend = T) +
#   labs(x = "Day of Year",
#        y = "Thousands of Chinook Salmon")+
#   # geom_hline(aes(yintercept = can.mean/1000,
#   #                col = "black"),
#   #            size = 1.5,
#   #            linetype = 2)+
#   
#   geom_hline(data = pf_hist, aes(yintercept = mean/1000), 
#              size = 1.5,
#              linetype = 2)+
#   
#   facet_wrap(~Year, scales="free_y")+
#   # scale_fill_calc()+
#   # scale_color_colorblind()+
#   # scale_color_continuous()
#   scale_x_continuous(breaks = c(153,173,193,213,232), labels = c(
#     "June 2",
#     
#     # "June 12",
#     
#     "June 22",
#     
#     # "July 2",
#     
#     "July 12",
#     
#     # "July 22",
#     
#     "Aug 1",
#     
#     "Aug 20"
#   ))+
#   # ggtitle(paste("Model Version", mod$`2007_153`$version,"\n Parameter = sigma normal curve"))+
#   theme_light()+
#   theme(legend.position = "top",
#         legend.title = element_blank(),
#         axis.text.x = element_text(angle = 90),
#         text = element_text(size = 18))

# Sigma for regression in normal curve fitting################

ggplot(sigma_plotDF[sigma_plotDF$Year==2015,], aes(x = day, y = median))+
  geom_ribbon(aes(ymin = median-0.003, ymax = median+0.003, fill = "Median")) +
  # geom_point(aes(
  #   col = param),
  #   size = 2, 
  #   show.legend = T)+
  geom_ribbon(aes(ymin=low95, ymax=up95, fill = "HDI 95"), alpha = .2,
              show.legend = T) +
  geom_ribbon(aes(ymin=low50,
                  ymax=up50, fill = "HDI 50"), alpha = .5) +
  facet_wrap(~param)+
  labs(x = "Day",
       y = bquote(sigma^2),
       fill = "")+
  scale_x_continuous(breaks = c(152,162,172,182,192,202),
                     labels = c("June-1",
                                "June-11",
                                "June-21",
                                "July-1",
                                "July-11",
                                "July-22"))+
  scale_fill_colorblind()+
  coord_cartesian(xlim = c(152,210))+
  theme_light(base_size = 15)+
  theme(axis.title    = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90))

outputlist_ver1RR_T[[1]]$pars$ln_prior_pf


years <- c(2022)
n.years <- length(years)
count <- 0
for (y in 1:n.years) {
  # y = 2
  year <- years[y]
  df <-interval_function(mod = outputList, 
                         parameter = "sigma_predPSS", # Parameter of interest
                         testDays = testDays,
                         year = year,
                         counter = count)
  if(y == 1){
    final.df <- df}else{
      
      final.df <- rbind(final.df,df)
      
    }
  count <- count+ length(testDays)
}


ggplot(final.df, aes(x = day, y = median))+
  geom_line(aes()) +
  geom_point(aes(),
    size = 2, 
    show.legend = T)+
  geom_ribbon(aes(ymin=low95, ymax=up95), alpha = .2,
              show.legend = T) +
  geom_ribbon(aes(ymin=low50,
                  ymax=up50), alpha = .5)
