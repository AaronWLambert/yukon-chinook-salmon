#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 05.25.23
#          
# Purpose: Script to generate retrospective plots from outputs for manuscript
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

# Functions #############################################################################

# Functions to run model, create basic plots, and get day of year
# source(file = file.path(dir.R,"model_run_function LOO3.r"))


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

# Days used in retrospective testing runs #################################################

# Test days used in full season run (every 5 days starting June 2)
testDays <- seq(from = 153, to = 243, by = 5)

# Years included in retrospective testing
testYears <- c(2007:2022)

# Read in outputlist ######################################################################

# Read in outputs from retro testing 
# **Note** Increase memory limit if using versions 4 and higher
memory.limit(size = 60000)


# Chapter 1 outputs #######################################################################
# Days tested are day 153 to 243 by 5 days
# Years Tested include 2007-2022
# name = version_reconstructionMethod_truncated years in regression
ver1_RR_T <- readRDS(file = file.path(dir.output, "Ver1_rr_trunc_FINAL_5June23.RDS"))

# ver1_RR_Full <- readRDS(file = file.path(dir.output, "Ver1_RRCAN_long_27Apr23.RDS"))

ver11_RR_T <- readRDS(file = file.path(dir.output, "Ver11_rr_trunc_final_13May23.RDS"))

ver101_RR_T <- readRDS(file = file.path(dir.output, "Ver101_rr_trunc_final_13May23.RDS"))

ver101if_RR_T <- readRDS(file = file.path(dir.output, "Ver101if_rr_trunc_final_13May23.RDS"))

ver102_RR_T <- readRDS(file = file.path(dir.output, "Ver102_rr_trunc_17July23.RDS"))

# ver112_RR_T <- readRDS(file = file.path(dir.output, "Ver112_rr_trunc_final_13May23.RDS"))

ver2_RR_T <- readRDS(file = file.path(dir.output,"Ver2_rr_trunc_final_13May23.RDS")) 

ver2c_RR_T <- readRDS(file = file.path(dir.output,"Ver2c_rr_trunc_final_13May23.RDS")) 

# Version fitting one logistic curve to get proportions
ver602_RR_T <- readRDS(file = file.path(dir.output,"Ver6propNull_final_27Sept23.RDS")) 

# Version fitting normal curve to all years
ver402_RR_T <- readRDS(file=file.path(dir.output,"Ver402_rr_trunc_FINAL_13July23.RDS"))

# Version fitting logistic curve to all years
ver502_RR_T <- readRDS(file = file.path(dir.output,"Ver502_rr_trunc_final_31Julu23.RDS"))

# ver2_RR <- readRDS(file = file.path(dir.output,"Ver2_RRCAN_long_1May23.RDS")) 

# ver2c_RR <- readRDS(file = file.path(dir.output,"Ver2c_RRCAN_long_1May23.RDS")) 


##### Calculations for retrospecitve testing ############################################
# Uses retrospective.function

#
RetroList_ver1RR_T <- retrospective.function(outputList = ver1_RR_T,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             pf = FALSE)

RetroList_ver1RR_F <- retrospective.function(outputList = ver1_RR_Full,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             pf = FALSE)

RetroList_ver11RR_T <- retrospective.function(outputList = ver11_RR_T,
                                              testYears = testYears,
                                              testDays = testDays,
                                              CAN_hist = CAN_hist_new,
                                              pf = FALSE)

RetroList_ver101RR_T <- retrospective.function(outputList = ver101_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

RetroList_ver101ifRR_T <- retrospective.function(outputList = ver101if_RR_T,
                                                 testYears = testYears,
                                                 testDays = testDays,
                                                 CAN_hist = CAN_hist_new,
                                                 pf = FALSE)

RetroList_ver102RR_T <- retrospective.function(outputList = ver102_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

RetroList_ver2RR_T <- retrospective.function(outputList = ver2_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

RetroList_ver2cRR_T <- retrospective.function(outputList = ver2c_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

RetroList_ver602RR_T <- retrospective.function(outputList = ver602_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

# rm(ver602_RR_T)

RetroList_ver402RR_T <- retrospective.function(outputList = ver402_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)
rm(ver402_RR_T)

RetroList_ver502RR_T <- retrospective.function(outputList = ver502_RR_T,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist_new,
                                           pf = FALSE)
rm(ver502_RR_T)

# RetroList_pf_old <- retrospective.function(outputList = outputlist_ver1,
#                                            testYears = testYears,
#                                            testDays = testDays,
#                                            CAN_hist = CAN_hist_old,
#                                            pf_hist = pf_hist,
#                                            startYearRetro = 2007,
#                                            endYearRetro = 2022,
#                                            pf = TRUE)

RetroList_pf_new <- retrospective.function(outputList = ver102_RR_T,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist_new,
                                           pf_hist = pf_hist,
                                           startYearRetro = 2007,
                                           endYearRetro = 2022,
                                           pf = TRUE)

# RMSE
# Extract RMSE by day into dataframe for each version and posterior median RunSize

rmseDF_ver1RR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver1RR_T$RMSE_by_day_vect,
                              "Version" = "1.0")

rmseDF_ver1RR_F <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver1RR_F$RMSE_by_day_vect,
                              "Version" = "1.0 PSS 1995-2022")

rmseDF_ver11RR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver11RR_T$RMSE_by_day_vect,
                              "Version" = "1.1")

rmseDF_ver101RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver101RR_T$RMSE_by_day_vect,
                                "Version" = "1.0.1")

rmseDF_ver101ifRR_T <- data.frame("Day" = testDays,
                                  "RMSE" = RetroList_ver101ifRR_T$RMSE_by_day_vect,
                                  "Version" = "1.0.1.IF")

rmseDF_ver102RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver102RR_T$RMSE_by_day_vect,
                                "Version" = "1.0.2")

rmseDF_ver112RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver112RR_T$RMSE_by_day_vect,
                                "Version" = "1.1.2")

rmseDF_ver2RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver2RR_T$RMSE_by_day_vect,
                                "Version" = "2.0")

rmseDF_ver2cRR_T <- data.frame("Day" = testDays,
                              "RMSE" = RetroList_ver2cRR_T$RMSE_by_day_vect,
                              "Version" = "2c.0")

rmseDF_ver402RR_T <- data.frame("Day" = testDays,
                            "RMSE" = RetroList_ver402RR_T$RMSE_by_day_vect,
                             "Version" = "3.0.2")

rmseDF_ver502RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver502RR_T$RMSE_by_day_vect,
                                "Version" = "4.0.2")

rmseDF_ver602RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver602RR_T$RMSE_by_day_vect,
                                "Version" = "5.0.2")
# 
rmseDF_pf_new <- data.frame("Day" = testDays,
                            "RMSE" = RetroList_pf_new$RMSE_by_day_vect,
                            "Version" = "PF")

# rmseDF_pf_old <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_pf_old$RMSE_by_day_vect,
#                             "Version" = "PF Eagle+Harv")

# Combine data frames into one data.frame for plotting
full_rmseDF <- rbind(
                     rmseDF_ver1RR_T,
                     rmseDF_ver2RR_T, 
                     rmseDF_ver2cRR_T,
                     rmseDF_ver11RR_T,
                     rmseDF_ver101RR_T,
                     rmseDF_ver101ifRR_T,
                     rmseDF_ver102RR_T,
                     rmseDF_ver402RR_T,
                     rmseDF_ver502RR_T,
                     rmseDF_ver602RR_T
                     # rmseDF_pf_new
                     )

# Make version a factor
full_rmseDF$Version <- factor(full_rmseDF$Version, levels = c("1.0",
                                                              "2.0",
                                                              "2c.0",
                                                              "1.1",
                                                              "1.0.1",
                                                              "1.0.1.IF",
                                                              "1.0.2",
                                                              "3.0.2",
                                                              "4.0.2",
                                                              "5.0.2"
                                                              # "PF"
                                                              ))

# Table for papers
RMSE_table <-as.data.frame(pivot_wider(full_rmseDF,names_from = Day, values_from = RMSE))

RMSE_table_round <- format(RMSE_table, digits = 2)

# write.table(x = RMSE_table_round, file = file.path(dir.output,"RMSE Table 13Jul23"), row.names = F)

# write.csv(x = RMSE_table_round, file = file.path(dir.output,"RMSE Table 17Jul23 2"), row.names = FALSE)
ggplot(full_rmseDF, aes(x = factor(Day), y = factor(Version), fill = RMSE))+
  geom_tile()+
  geom_text(aes(label = round(RMSE)))+
  scale_fill_distiller(palette = "Spectral")+
  labs(x = "Day of Year",
       y = "Model Versions")+
  scale_x_discrete(breaks = c(testDays), labels = c(
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
  theme(text = element_text(size = 14)
        )

# Add PF for plotting
full_rmseDF$PF <- rmseDF_pf_new$RMSE


# png(filename = file.path(dir.figs,"RMSE Plot 5Oct23.png"), width = 1000, height = 800)
# RMSE plot
ggplot(full_rmseDF, aes(x = Day, y = RMSE
                        # col = Version,
                        # shape= Version, 
                        # alpha = 0.7
                        ))+
  # geom_col(position = "dodge",width = 3) +
  geom_point()+
  geom_line()+
  # labs(fill = "", col = "")+
  # geom_col(data = short_rmseDF, position = "dodge", width = 3)+
  coord_cartesian(ylim = c(
    min(full_rmseDF$RMSE)-500,
    max(full_rmseDF$RMSE)+500
    # 10000,15000
    
  ))+
  # scale_color_colorblind(name = "",
  #                        labels = c("No Eagle", "Eagle Regression", "Eagle Regression Back-half","Eagle Prop Est","PF"))+
  # # scale_colour_manual(name = "",
  # #                     #                     # labels = c("PF","Ver 2.c", "Ver 2.c.0.1", "Ver 2.c.sd"),
  # #                     # labels = c("PF","No Eagle", "Eagle Prop Est", "GSI No Eagle","GSI & Eagle Prop Est"),
  # #                     labels = c("PF","Ver 1", "SST Est Prop"),
  # #                     values = c("red","gold","green","blue", "purple")) +
  # # scale_shape_manual(name = "",
  # #                    labels = c("No Eagle", "Eagle Regression", "Eagle Regression Back-half","Eagle Prop Est","PF"),
  # #                    # labels = c("PF","Ver 1", "SST Est Prop"),
  # #                    values = c(16,17,18,19,20))+
  # # # scale_fill_manual(values = wes_palette("IsleofDogs1"),
  # #                   labels = c("PF (New; SS Recon.)",
  # #                              "PF (Old; Eagle+ Harvest)",
  # #                              "2.1 (New; SS Recon.)",
  # #                              "2.1 (Old; Eagle + Harvest)"))+
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
  geom_hline(aes(yintercept = PF, color = factor(PF)),size = 1,linetype = 2, show.legend = T)+
  facet_wrap(~Version, ncol = 3)+
  scale_color_discrete(labels = "Preseason Forecast",
                       name = "")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 0.5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90))


# dev.off()


# RMSE Plot comparing full vs truncated PSS data set

# Combine data frames into one data.frame for plotting
full_rmseDF_supp1 <- rbind(
  rmseDF_ver1RR_T,
  rmseDF_ver1RR_F
  # rmseDF_pf_new
)
full_rmseDF_supp1$PF <- rmseDF_pf_new$RMSE
# png(filename = file.path(dir.figs,"RMSE Plot 31July23.png"), width = 1000, height = 800)

# RMSE plot
ggplot(full_rmseDF_supp1, aes(x = Day, y = RMSE
                        # col = Version,
                        # shape= Version, 
                        # alpha = 0.7
))+
  # geom_col(position = "dodge",width = 3) +
  geom_point()+
  geom_line()+
  # labs(fill = "", col = "")+
  # geom_col(data = short_rmseDF, position = "dodge", width = 3)+
  coord_cartesian(ylim = c(
    min(full_rmseDF$RMSE)-500,
    max(full_rmseDF$RMSE)+500
    # 10000,15000
    
  ))+
  # scale_color_colorblind(name = "",
  #                        labels = c("No Eagle", "Eagle Regression", "Eagle Regression Back-half","Eagle Prop Est","PF"))+
  # # scale_colour_manual(name = "",
  # #                     #                     # labels = c("PF","Ver 2.c", "Ver 2.c.0.1", "Ver 2.c.sd"),
  # #                     # labels = c("PF","No Eagle", "Eagle Prop Est", "GSI No Eagle","GSI & Eagle Prop Est"),
  # #                     labels = c("PF","Ver 1", "SST Est Prop"),
  # #                     values = c("red","gold","green","blue", "purple")) +
  # # scale_shape_manual(name = "",
  # #                    labels = c("No Eagle", "Eagle Regression", "Eagle Regression Back-half","Eagle Prop Est","PF"),
  # #                    # labels = c("PF","Ver 1", "SST Est Prop"),
  # #                    values = c(16,17,18,19,20))+
# # # scale_fill_manual(values = wes_palette("IsleofDogs1"),
# #                   labels = c("PF (New; SS Recon.)",
# #                              "PF (Old; Eagle+ Harvest)",
# #                              "2.1 (New; SS Recon.)",
# #                              "2.1 (Old; Eagle + Harvest)"))+
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
  geom_hline(aes(yintercept = PF, color = factor(PF)),size = 1,linetype = 2, show.legend = T)+
  facet_wrap(~Version)+
  scale_color_discrete(labels = "Preseason Forecast",
                       name = "")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 0.5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90))


# dev.off()
# MAPE ####################################################################

# Create MAPE dataframe
mapeDF_ver1RR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver1RR_T$MAPE_vect,
                              "Version" = "1.0")

mapeDF_ver11RR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver11RR_T$MAPE_vect,
                              "Version" = "1.1")

mapeDF_ver101RR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver101RR_T$MAPE_vect,
                              "Version" = "1.0.1")

mapeDF_ver101ifRR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver101ifRR_T$MAPE_vect,
                              "Version" = "1.0.1.IF")

mapeDF_ver102RR_T <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver102RR_T$MAPE_vect,
                              "Version" = "1.0.2")

mapeDF_ver2RR <- data.frame("Day" = testDays,
                            "MAPE" = RetroList_ver2RR_T$MAPE_vect,
                            "Version" = "2.0")

mapeDF_ver2cRR <- data.frame("Day" = testDays,
                             "MAPE" = RetroList_ver2cRR_T$MAPE_vect,
                             "Version" = "2c.0")

mapeDF_ver402RR <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver402RR_T$MAPE_vect,
                              "Version" = "3.0.2")

mapeDF_ver502RR <- data.frame("Day" = testDays,
                          "MAPE" = RetroList_ver502RR_T$MAPE_vect,
                          "Version" = "4.0.2")

mapeDF_ver602RR <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver602RR_T$MAPE_vect,
                              "Version" = "5.0.2")

mapeDF_PF_new <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_pf_new$MAPE_vect,
                            "Version" = "PF")

# Combine into one data-frame
full_MAPE.df <- rbind(
                      mapeDF_ver1RR_T,
                      mapeDF_ver11RR_T,
                      mapeDF_ver101RR_T,
                      mapeDF_ver101ifRR_T,
                      mapeDF_ver102RR_T,
                      mapeDF_ver2RR,
                      mapeDF_ver2cRR,
                      mapeDF_ver402RR,
                      mapeDF_ver502RR,
                      mapeDF_ver602RR
                      # mapeDF_PF_new
                      )


# Make version a factor
full_MAPE.df$Version <- factor(full_MAPE.df$Version, levels = c("1.0",
                                                              "2.0",
                                                              "2c.0",
                                                              "1.1",
                                                              "1.0.1",
                                                              "1.0.1.IF",
                                                              "1.0.2",
                                                              "3.0.2",
                                                              "4.0.2",
                                                              "5.0.2"
                                                              # "PF"
))

# MAPE Table for papers
# MAPE_table<-as.data.frame(pivot_wider(full_MAPE.df,names_from = Day, values_from = MAPE))
# MAPE_table_round <- format(MAPE_table, digits = 3)

# write.table(x = MAPE_table_round, file = file.path(dir.output,"MAPE Table 1Jun23"),
            # row.names = F,sep = ",")

# Add PF for plotting
full_MAPE.df$PF <- mapeDF_PF_new$MAPE

ggplot(full_MAPE.df, aes(x = factor(Day), y = factor(Version), fill = MAPE*100))+
  geom_tile()+
  geom_text(aes(label = round(MAPE*100,1)))+
  scale_fill_distiller(palette = "Spectral", name = "MAPE")+
  scale_x_discrete(breaks = c(testDays), labels = c(
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
  labs(x = "Day of Year",
       y = "Model Versions")+
  theme(text = element_text(size = 14))

# png(file = file.path(dir.figs,"MAPE Plot 5Oct23.png"),width = 1000, height = 800)
# MAPE plot
ggplot(full_MAPE.df, aes(x = Day,
                         y = MAPE*100, 
                         # col = Version,
                         # shape = Version
                         ))+
  # geom_col(position = "dodge",width = 3) + 
  
  geom_point(size = 1.5)+
  geom_line(aes(group = Version))+
  geom_hline(aes(yintercept = PF*100, color = factor(PF)),size = 1,linetype = 2, show.legend = T)+
  coord_cartesian(ylim = c(min(full_MAPE.df$MAPE*100)-2,max(full_MAPE.df$MAPE*100)+2))+
  ylab('MAPE (%)')+
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
  facet_wrap(~Version, ncol = 3)+
  scale_color_discrete(labels = "Preseason Forecast",
                       name = "")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = .5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90))

# dev.off()
# Percent Error Box Plots

PE_ver1 <- as.data.frame(RetroList_ver1RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver1$version <- "1.0"

PE_ver11 <- as.data.frame(RetroList_ver11RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver11$version <- "1.1"

PE_ver101 <- as.data.frame(RetroList_ver101RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver101$version <- "1.0.1"

PE_ver101if <- as.data.frame(RetroList_ver101ifRR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver101if$version <- "1.0.1.IF"

PE_ver102 <- as.data.frame(RetroList_ver102RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver102$version <- "1.0.2"

PE_ver2 <- as.data.frame(RetroList_ver2RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2$version <- "2.0"

PE_ver2c <- as.data.frame(RetroList_ver2cRR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2c$version <- "2c.0"

PE_ver402 <- as.data.frame(RetroList_ver402RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver402$version <- "3.0.2"

PE_ver502 <- as.data.frame(RetroList_ver502RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver502$version <- "4.0.2"

PE_ver602 <- as.data.frame(RetroList_ver602RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver602$version <- "5.0.2"

PE_pf_new <- as.data.frame(RetroList_pf_new$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_pf_new$version <- "PF New"

peDF_total <- rbind(
                    PE_ver1,
                    PE_ver11,
                    PE_ver101,
                    PE_ver101if,
                    PE_ver102,
                    PE_ver2,
                    PE_ver2c,
                    PE_ver402,
                    PE_ver502,
                    PE_ver602
  
)

# For plotting PF point shape 

newDF <- left_join(x = pf_hist, y = CAN_hist_new)
# newDF <- left_join(x = newDF, y = CAN_hist_old, "Year")

# newDF$PF_new <- (newDF$Weighted_Forecast-newDF$can.mean.x)/newDF$can.mean.x

newDF$Year <- as.factor(newDF$Year)

newDF$PF <- (newDF$mean-newDF$can.mean)/newDF$can.mean



peDF_total<-left_join(peDF_total, newDF[,c(2,4)])

peDF_total$version <- factor(peDF_total$version, levels = c("1.0",
                                                              "2.0",
                                                              "2c.0",
                                                              "1.1",
                                                              "1.0.1",
                                                              "1.0.1.IF",
                                                              "1.0.2",
                                                              "3.0.2",
                                                              "4.0.2",
                                                            "5.0.2"))
# Colors for plot
cust.col.pe <- c("black",
              "#E69F00",
              "#56B4E9",
              "#009E73",
              "#F0E442",
              "#0072B2",
              "#D55E00",
              "#CC79A7",
              "darkgray",
              "white",
              "red")
# peDF_total <- rbind(peDF_total, newDF[,c(2,4,5)])
# png(file = file.path(dir.figs,"PE Plot 5Oct23.png"),width = 1000, height = 800)
# long_df <- newDF %>% pivot_longer(cols = c(PF_new,PF_old)) %>% as.data.frame()
ggplot(peDF_total, aes(x = Year, y = PE*100, 
                       fill =version))+
  guides(alpha = "none")+
  geom_boxplot(alpha = .65, position = position_dodge(width = 1.01)) + 
  labs(fill = "Model Version",
       y = "Percent Error",
       x = "",
       color = "")+
  geom_point(
             aes( x = Year, y = PF*100, col = "PF"),fill = "red",
             shape = 21,
             position = position_nudge(x = -.55),
             size = 3)+
  # scale_fill_colorblind()+
  scale_fill_manual(values = cust.col.pe)+
  # scale_color_colorblind(name = "")+
geom_hline(yintercept  = 0, linetype = 2, color = "orange", size = 1)+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = .1, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text.x =   element_blank(),
        axis.ticks.x = element_blank())+
  facet_wrap(~Year, scales = "free")
# ggtitle("New Can_hist ver 2")
# dev.off()
# CV precision of posterior ##################################

# Function to calculate the CV of the 
# posterior run size projection

#' Title 
#'CV of the posterior run size projection
#' @param mod a model output from retrospective testing 
#'
#' @return a matrix with the cv of run size estimate, the year, 
#' the day of the projection, and the model version.
#' 
cv.func <- function(mod){
  
  cv.mat <- matrix(nrow = length(mod), ncol = 4)
  colnames(cv.mat) <- c("cv","year", "day", "model")
  
  for(i in 1:length(mod)){
    
    cv.mat[i,1] <- sd(mod[[i]]$pars$RunSize)/median(mod[[i]]$pars$RunSize)
    cv.mat[i,2] <- mod[[i]]$myYear
    cv.mat[i,3] <- mod[[i]]$myDay
    cv.mat[i,4] <- (mod[[i]]$version)
  }
  
  cv.df <- as.data.frame(cv.mat)
  
  return(cv.df)
}

cv.df.10 <- cv.func(mod = ver1_RR_T)

cv.df.2c <- cv.func(mod = ver2c_RR_T)

cv.df.2 <- cv.func(mod = ver2_RR_T)

cv.df.11 <- cv.func(mod = ver11_RR_T)

cv.df.101 <- cv.func(mod = ver101_RR_T)

cv.df.101if <- cv.func(mod = ver101if_RR_T)

cv.df.102 <- cv.func(mod = ver102_RR_T)

ver402_RR_T <- readRDS(file=file.path(dir.output,"Ver402_rr_trunc_FINAL_29May23.RDS"))

cv.df.402 <- cv.func(mod = ver402_RR_T)

rm(ver402_RR_T)

ver502_RR_T <- readRDS(file = file.path(dir.output,"Ver502_rr_trunc_final_31Julu23.RDS"))

cv.df.502 <- cv.func(mod = ver502_RR_T)

rm(ver502_RR_T)

# Version fitting one logistic curve to get proportions
ver602_RR_T <- readRDS(file = file.path(dir.output,"Ver6propNull_final_27Sept23.RDS")) 

cv.df.602 <- cv.func(mod = ver602_RR_T)

cv.df.101if[,4] <- "1.0.1.IF"

cv.df.402[,4] <- "3.0.2"

cv.df.502[,4] <- "4.0.2"

cv.df.602[,4] <- "5.0.2"

cv.total <- rbind(
                  cv.df.10,
                  cv.df.2,
                  cv.df.2c,
                  cv.df.11,
                  cv.df.101,
                  cv.df.101if,
                  cv.df.102,
                  cv.df.402,
                  cv.df.502,
                  cv.df.602)

# str(cv.df.102)

cv.total$cv <- as.numeric(cv.total$cv)
cv.total$year <- as.numeric(cv.total$year)
cv.total$day <- as.numeric(cv.total$day)
cv.total$model <- factor(cv.total$model, levels = c("1.0",
                                                    "2.0",
                                                    "2.c",
                                                    "1.1",
                                                    "1.0.1",
                                                    "1.0.1.IF",
                                                    "1.0.2",
                                                    "3.0.2",
                                                    "4.0.2",
                                                    "5.0.2"))


# png(file = file.path(dir.figs,"CV Plot 5Oct23.png"),width = 1300, height = 800)

cust.col <- c("black",
              "#E69F00",
              "#56B4E9",
              "#009E73",
              "#F0E442",
              "#0072B2",
              "#D55E00",
              "#CC79A7",
              "darkgray",
              "purple")

ggplot(cv.total, aes(x = day, y = cv*100, 
                     col = model))+
  geom_line(size = .75)+
  geom_point(size = 2)+
  labs(x = "Day",
       y = "CV (%)")+
  # scale_color_colorblind()+
  scale_color_manual(values = cust.col)+
  # facet_grid(model~year)+
  facet_wrap(~year)+
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
        panel.grid.major.y  =element_line(linetype = 2, size = .1, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text.x =   element_text(angle = 90, size = 16)
        # axis.ticks.x = element_blank()
  )

# dev.off()

ggplot(cv.total, aes(x = day, y = model, fill = cv))+
  geom_tile()+
  geom_text(aes(label = round(cv*100,1)),size = 2)+
  scale_fill_distiller(palette = "Spectral", name = "CV")+
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
  labs(x = "Day of Year",
       y = "Model Versions")+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90))+
  facet_wrap(~year)
# Turn off PDF save
# dev.off()

# Plots comparing PSS prediction, full run size projection and pf projection ###############################
#
RetroList_ver1RR_T_PSS <- retrospective.function(outputList = ver1_RR_T,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             param = "post_curr_predPSS",
                                             pf = FALSE)

RetroList_ver1RR_F_PSS <- retrospective.function(outputList = ver1_RR_Full,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             param = "post_curr_predPSS",
                                             pf = FALSE)

RetroList_ver11RR_T_PSS <- retrospective.function(outputList = ver11_RR_T,
                                              testYears = testYears,
                                              testDays = testDays,
                                              CAN_hist = CAN_hist_new,
                                              param = "post_curr_predPSS",
                                              pf = FALSE)

RetroList_ver101RR_T_PSS <- retrospective.function(outputList = ver101_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               param = "post_curr_predPSS",
                                               pf = FALSE)

RetroList_ver101ifRR_T_PSS <- retrospective.function(outputList = ver101if_RR_T,
                                                 testYears = testYears,
                                                 testDays = testDays,
                                                 CAN_hist = CAN_hist_new,
                                                 param = "post_curr_predPSS",
                                                 pf = FALSE)

RetroList_ver102RR_T_PSS <- retrospective.function(outputList = ver102_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               param = "post_curr_predPSS",
                                               pf = FALSE)

RetroList_ver2RR_T_PSS <- retrospective.function(outputList = ver2_RR_T,
                                             testYears = testYears,
                                             testDays = testDays,
                                             CAN_hist = CAN_hist_new,
                                             param = "post_curr_predPSS",
                                             pf = FALSE)

RetroList_ver2cRR_T_PSS <- retrospective.function(outputList = ver2c_RR_T,
                                              testYears = testYears,
                                              testDays = testDays,
                                              CAN_hist = CAN_hist_new,
                                              param = "post_curr_predPSS",
                                              pf = FALSE)

RetroList_ver602RR_T_PSS <- retrospective.function(outputList = ver602_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               param = "post_curr_predPSS",
                                               pf = FALSE)

# rm(ver602_RR_T)

RetroList_ver402RR_T_PSS <- retrospective.function(outputList = ver402_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               param = "post_curr_predPSS",
                                               pf = FALSE)
# rm(ver402_RR_T)

RetroList_ver502RR_T_PSS <- retrospective.function(outputList = ver502_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               param = "post_curr_predPSS",
                                               pf = FALSE)


# Create MAPE dataframe for PSS posterior estimates
mapeDF_ver1RR_T_PSS <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver1RR_T_PSS$MAPE_vect,
                              "Version" = "1.0")

mapeDF_ver11RR_T_PSS <- data.frame("Day" = testDays,
                               "MAPE" = RetroList_ver11RR_T_PSS$MAPE_vect,
                               "Version" = "1.1")

mapeDF_ver101RR_T_PSS <- data.frame("Day" = testDays,
                                "MAPE" = RetroList_ver101RR_T_PSS$MAPE_vect,
                                "Version" = "1.0.1")

mapeDF_ver101ifRR_T_PSS <- data.frame("Day" = testDays,
                                  "MAPE" = RetroList_ver101ifRR_T_PSS$MAPE_vect,
                                  "Version" = "1.0.1.IF")

mapeDF_ver102RR_T_PSS <- data.frame("Day" = testDays,
                                "MAPE" = RetroList_ver102RR_T_PSS$MAPE_vect,
                                "Version" = "1.0.2")

mapeDF_ver2RR_PSS <- data.frame("Day" = testDays,
                            "MAPE" = RetroList_ver2RR_T_PSS$MAPE_vect,
                            "Version" = "2.0")

mapeDF_ver2cRR_PSS <- data.frame("Day" = testDays,
                             "MAPE" = RetroList_ver2cRR_T_PSS$MAPE_vect,
                             "Version" = "2c.0")

mapeDF_ver402RR_PSS <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver402RR_T_PSS$MAPE_vect,
                              "Version" = "3.0.2")

mapeDF_ver502RR_PSS <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver502RR_T_PSS$MAPE_vect,
                              "Version" = "4.0.2")

mapeDF_ver602RR_PSS <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver602RR_T_PSS$MAPE_vect,
                              "Version" = "5.0.2")

mapeDF_PF_new <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_pf_new$MAPE_vect,
                            "Version" = "PF")

# Combine into one data-frame
full_MAPE.PSS.df <- rbind(
  mapeDF_ver1RR_T_PSS,
  mapeDF_ver11RR_T_PSS,
  mapeDF_ver101RR_T_PSS,
  mapeDF_ver101ifRR_T_PSS,
  mapeDF_ver102RR_T_PSS,
  mapeDF_ver2RR_PSS,
  mapeDF_ver2cRR_PSS,
  mapeDF_ver402RR_PSS,
  mapeDF_ver502RR_PSS,
  mapeDF_ver602RR_PSS
  # mapeDF_PF_new
)


# Make version a factor
full_MAPE.PSS.df$Version <- factor(full_MAPE.PSS.df$Version, levels = c("1.0",
                                                                "2.0",
                                                                "2c.0",
                                                                "1.1",
                                                                "1.0.1",
                                                                "1.0.1.IF",
                                                                "1.0.2",
                                                                "3.0.2",
                                                                "4.0.2",
                                                                "5.0.2"
                                                                # "PF"
))

# Add column with PSS or RunSize
full_MAPE.df$Estimate <- "Integrated Run Size"

full_MAPE.PSS.df$Estimate <- "PSS Run Size"

# Join the df
full_MAPE_comp.df <- rbind(full_MAPE.df[,-4],full_MAPE.PSS.df)


# Add PF for plotting
full_MAPE_comp.df$PF <- mapeDF_PF_new$MAPE

# png(file = file.path(dir.figs,"MAPE Comp Plot 5Oct23.png"),width = 1000, height = 800)
# MAPE plot
ggplot(full_MAPE_comp.df, aes(x = Day,
                         y = MAPE*100, 
                         col = Estimate,
                         # shape = Version
))+
  # geom_col(position = "dodge",width = 3) + 
  
  geom_point(size = 1.5)+
  geom_line()+
  geom_hline(aes(yintercept = PF*100, color = factor(PF)),size = 1,linetype = 2, show.legend = T)+
  coord_cartesian(ylim = c(min(full_MAPE_comp.df$MAPE*100)-2,max(full_MAPE_comp.df$MAPE*100)+2))+
  ylab('MAPE (%)')+
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
  facet_wrap(~Version, ncol = 3)+
  scale_color_colorblind(labels = c("Preseason Forecast",
                                    "Posterior Integrated Run Size", 
                                    "Posterior PSS Prediction"),
                       name = "")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = .5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90))

# dev.off()

# Plots of median run size prediction over each year ##################################################

# Use the runsize_trend_func to extract median projections from models

runTrend_ver1 <- runsize_trend_func(ver1_RR_T,
                                    testYears = testYears,
                                    testDays = testDays ,
                                    CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver11 <- runsize_trend_func(ver11_RR_T,
                                     testYears = testYears,
                                     testDays = testDays ,
                                     CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver101 <- runsize_trend_func(ver101_RR_T,
                                      testYears = testYears,
                                      testDays = testDays ,
                                      CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver101if <- runsize_trend_func(ver101if_RR_T,
                                        testYears = testYears,
                                        testDays = testDays ,
                                        CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver102 <- runsize_trend_func(ver102_RR_T,
                                      testYears = testYears,
                                      testDays = testDays ,
                                      CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver2   <- runsize_trend_func(ver2_RR_T,
                                      testYears = testYears,
                                      testDays = testDays ,
                                      CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver2c  <- runsize_trend_func(ver2c_RR_T,
                                      testYears = testYears,
                                      testDays = testDays ,
                                      CAN_hist = CAN_hist_new,PF = pf_hist)

ver402_RR_T <- readRDS(file=file.path(dir.output,"Ver402_rr_trunc_FINAL_29May23.RDS"))

runTrend_ver402 <- runsize_trend_func(ver402_RR_T,
                                      testYears = testYears,
                                      testDays = testDays ,
                                      CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver402$Version <- "3.0.2"

rm(ver402_RR_T)
gc()

ver502_RR_T <- readRDS(file = file.path(dir.output,"Ver502_rr_trunc_final_31Julu23.RDS"))

runTrend_ver502 <- runsize_trend_func(ver502_RR_T,
                                      testYears = testYears,
                                      testDays = testDays ,
                                      CAN_hist = CAN_hist_new,PF = pf_hist)

runTrend_ver502$Version <- "4.0.2"

rm(ver502_RR_T)


comb.df <- rbind(
                  runTrend_ver1,
                  runTrend_ver2,
                  runTrend_ver2c,
                  runTrend_ver11,
                  runTrend_ver101,
                  runTrend_ver101if,
                  runTrend_ver102,
                  runTrend_ver402,
                  runTrend_ver502
)

comb.df$Version <- factor(comb.df$Version, levels = c("1.0",
                                                     "2.0",
                                                     "2.c",
                                                     "1.1",
                                                     "1.0.1",
                                                     "1.0.1.if",
                                                     "1.0.2",
                                                     "3.0.2",
                                                     "4.0.2"))

# Plot the error
comb.df$error <- comb.df$RunSize - comb.df$can.mean

# Custom colors for run trend plots
cust.col2 <- c("black",
              "#E69F00",
              "#56B4E9",
              "#009E73",
              "#F0E442",
              "#0072B2",
              "#D55E00",
              "#CC79A7",
              "darkgray",
              "red")

# Plot of median run size predictions across days for each year
ggplot(comb.df, aes(x = Day, 
                    y = RunSize/1000,
                    color = Version
                    
))+
  geom_line(size = 1.5, show.legend = T, alpha = 0.75)+
  # geom_line( aes(x = Day, y = PSSpred/1000, color = "PSS Pred" ))+
  # geom_line(data = PSSpred_DF_mult, aes(x = Day, y = PSSpred/1000, color = "PSS Pred 2x"))+
  # geom_line(data = RunSize_DF_mult, aes(x = Day, y = RunSize/1000, color = "RunSize 2x"))+
  geom_hline(aes(yintercept = can.mean/1000 , col = "Truth"), size = 1, linetype = 2)+
  facet_wrap(~Year ,
             nrow = 4,
             scales = "free_y"
  )+
  theme(legend.position = "top")+
  labs(color = "", alpha = "")+
  ylab("Candadian-origin Chinook Salmon (1000s)")+
  # scale_color_colorblind(name = "",
  #                        labels = c(
  #                                   "PSS Only ",
  #                                   "PSS and Eagle",
  #                                   "EOS Run Size"))+
  # scale_color_colorblind()+
  scale_color_manual(values = cust.col2)+
  scale_x_continuous(breaks = c(153,173,193,213,232), labels = c(
    "June 2",
    
    # "June 12",
    
    "June 22",
    
    # "July 2",
    
    "July 12",
    
    # "July 22",
    
    "Aug 1",
    
    "Aug 20"))+
 
theme(legend.position = "right",
      panel.grid.major.y  =element_line(linetype = 2, size = .1, 
                                        color = "grey"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      plot.background = element_blank(),
      panel.background = element_blank(),
      text = element_text(size = 12),
      axis.text.x =   element_text(angle = 90),
      # axis.ticks.x = element_blank()
      )

# Plot the erverror
ggplot(comb.df, aes(x = Day, 
                    y = (error/1000),
                    color = Version
                    
))+
  geom_line(size = 1.5, show.legend = T, alpha = 0.5)+
  geom_hline(aes(yintercept = 0 , col = "Truth"), size = 1, linetype = 2)+
  facet_wrap(~Year ,
             nrow = 4,
             scales = "free_y"
  )+
  theme(legend.position = "top")+
  labs(color = "", alpha = "")+
  ylab("Predicted-Observed (1000s)")+
  # scale_color_colorblind(name = "",
  #                        labels = c(
  #                                   "PSS Only ",
  #                                   "PSS and Eagle",
  #                                   "EOS Run Size"))+
  # scale_color_colorblind()+
  
  scale_color_manual(values = cust.col2)+
  scale_x_continuous(breaks = c(153,173,193,213,232), labels = c(
    "June 2",
    
    # "June 12",
    
    "June 22",
    
    # "July 2",
    
    "July 12",
    
    # "July 22",
    
    "Aug 1",
    
    "Aug 20"))+
  
  theme(legend.position = "right",
        panel.grid.major.y  =element_line(linetype = 2, size = .1, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 12),
        axis.text.x =   element_text(angle = 90),
        # axis.ticks.x = element_blank()
  )


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
                    probs=c(0.025, 0.1,0.25,0.5, 0.75,0.9, 0.975),)
  
  pred.xx.df <- data.frame(testDays,
                           t(quant.xx))
  
  names(pred.xx.df) <- c("day","low95","low80","low50","median",
                         "up50","up80","up95")
  
  pred.xx.df$Year <- year
  
  # pred.xx.df$Year <- as.factor(pred.xx.df$Year)
  
  return(pred.xx.df)
}# End of function

###################################################


# Loop through years to calculate intervals for each year and day
# for run size

# ver 1.0.2
years <- c(testYears)
n.years <- length(years)
count <- 0
for (y in 1:n.years) {
  # y = 2
  year <- years[y]
  df <-interval_function(mod = ver102_RR_T, 
                         parameter = "RunSize", # Parameter of interest
                         testDays = testDays,
                         year = year,
                         counter = count)
  if(y == 1){
    final.df <- df}else{
      
      final.df <- rbind(final.df,df)
      
    }
  count <- count+ length(testDays)
}

final.df.102 <- left_join(final.df,CAN_hist_new)
final.df.102$version <- "1.0.2"

# ver 4.0.2
years <- c(testYears)
n.years <- length(years)
count <- 0
for (y in 1:n.years) {
  # y = 2
  year <- years[y]
  df <-interval_function(mod = ver402_RR_T, 
                         parameter = "RunSize", # Parameter of interest
                         testDays = testDays,
                         year = year,
                         counter = count)
  if(y == 1){
    final.df <- df}else{
      
      final.df <- rbind(final.df,df)
      
    }
  count <- count+ length(testDays)
}

final.df.402 <- left_join(final.df,CAN_hist_new)
final.df.402$version <- "4.0.2"


master.df <- rbind(final.df.102, final.df.402)

# Sigma from PSS regression################

ggplot(master.df, aes(x = day, y = median/1000))+
  # geom_ribbon(aes(ymin = median-0.003, ymax = median+0.003, fill = "Median")) +
  geom_line(aes(x = day, y = median/1000))+
  geom_ribbon(aes(ymin=low80/1000, ymax=up80/1000, fill = "HDI 80"), alpha = .5,
              show.legend = T) +
  geom_ribbon(aes(ymin=low50/1000,
                  ymax=up50/1000, fill = "HDI 50"), alpha = .5) +
  geom_hline(aes(yintercept = can.mean/1000, col = "Realized Run Size"), linetype = 2)+
  # facet_wrap(~Year, scales = "free_y")+
  facet_grid(version~Year)+
  labs(x = "Day",
       y = "Run Size (1000s Chinook Salmon)",
       fill = "",
       col = ""
       # title = "4.0.2"
       )+
  scale_x_continuous(breaks = c(152,162,172,182,192,202),
                     labels = c("June-1",
                                "June-11",
                                "June-21",
                                "July-1",
                                "July-11",
                                "July-22"))+
  scale_fill_colorblind()+
  coord_cartesian(xlim = c(152,210))+
  # theme_minimal(base_size = 15)+
  # theme(axis.title    = element_text(size = 15),
  #       legend.text = element_text(size = 12),
  #       axis.text.x = element_text(angle = 90),
  #       legend.position = "top")
  theme(legend.position = "top",
        # panel.grid.major.y  =element_line(linetype = 2, size = .1, 
        #                                   color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 18),
        axis.text.x =   element_text(angle = 90, size = 16))


# Rhat and Neff

# Matrix for values
summary.mat <- matrix(nrow = 9, ncol = 3)

# Name the rows and columns
colnames(summary.mat) <- c("Model", "Min n_eff", "Max Rhat")
rownames(summary.mat) <- c("Ver 1.0", 
                           "Ver 1.1",
                           "Ver 2.0",
                           "Ver 2c.0",
                           "Ver 1.0.1",
                           "Ver 1.0.1.if",
                           "Ver 1.0.2",
                           "Ver 4.0.2",
                           "Ver 5.0.2")

# Version 1.0
mod.name <- c("ver1_RR_T", "ver101_RR_T", "ver101if_RR_T", "ver102_RR_T")

summary.mat["Ver 1.0","Min n_eff"] <- min(ver1_RR_T[[]]$summary$summary[,"n_eff"],na.rm = TRUE)

which(ver1_RR_T[[1]]$summary$summary[,"n_eff"]==min.value)

i <-(mod.name[2])
call(i)
ver1_RR_T

for(i in 100:130){
  
  min.value <-min(ver502_RR_T[[i]]$summary$summary[,"n_eff"],na.rm = TRUE)
  
  
  pp<- which(ver502_RR_T[[i]]$summary$summary[,"n_eff"]==min.value)
  
  print(paste0(ver502_RR_T[[i]]$myYear,"    ",ver502_RR_T[[i]]$myDay,"      ", round(min.value),"    ",names(pp)))
  
}
# SCRATCH ##################################################################
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
sigma_func <- function(mod, testYears, testDays, sigma = sigma){
  
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
for (i in 1:length(outputlist_ver1_T)){
  
  
  plots_list[[i]] <- outPlots(outputList = outputlist_ver1RR_T[[i]],
                              CAN_hist = CAN_hist_new,
                              GSI = F,
                              Retrospective = TRUE,
                              eagle = FALSE)
}

# Vector of names form plot list
test<-names(outputlist_ver1_T)
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

figure<-ggarrange(plots_list$`2022_158`$PredPlot + labs(x = "", y ="")+
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




sd(ver102_RR_T$"2022_153"$pars$RunSize)/median(ver102_RR_T$"2022_153"$pars$RunSize)
