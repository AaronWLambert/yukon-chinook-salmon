#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 08.9.23 
# Chapter 2
#          
# Purpose: Script to generate retrospective plots from outputs for manuscript Chapter 2
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
source(file = file.path(dir.R,"Plot Function 9Aug23.R"))
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


# Chapter 2 outputs #######################################################################
# Days tested are day 153 to 243 by 5 days
# Years Tested include 2007-2022
ver102_RR_T <- readRDS(file = file.path(dir.output, "Ver102_rr_trunc_17July23.RDS"))

RetroList_ver102RR_T <- retrospective.function(outputList = ver102_RR_T,
                                                  testYears = testYears,
                                                  testDays = testDays,
                                                  CAN_hist = CAN_hist_new,
                                                  pf = FALSE)


# Model Version 102 with sst
ver102sst_RR_T <- readRDS(file = file.path(dir.output, "Ver102sstpoint_13Sept23.RDS"))

RetroList_ver102sstRR_T <- retrospective.function(outputList = ver102sst_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)

ver602_RR_T <- readRDS(file = file.path(dir.output, "Ver6prop2sstpoint_14Sept23.RDS"))

RetroList_ver602sstRR_T <- retrospective.function(outputList = ver602_RR_T,
                                                  testYears = testYears,
                                                  testDays = testDays,
                                                  CAN_hist = CAN_hist_new,
                                                  pf = FALSE)

ver602null_RR_T <- readRDS(file = file.path(dir.output, "6prop_null_5Sept23_point_sst.RDS"))

RetroList_ver602_null_RR_T <- retrospective.function(outputList = ver602null_RR_T,
                                                     testYears = testYears,
                                                     testDays = testDays,
                                                     CAN_hist = CAN_hist_new,
                                                     pf = FALSE)

ver402_RR_T <- readRDS(file=file.path(dir.output,"Ver402_rr_trunc_final_13July23.RDS"))

RetroList_ver402RR_T <- retrospective.function(outputList = ver402_RR_T,
                                               testYears = testYears,
                                               testDays = testDays,
                                               CAN_hist = CAN_hist_new,
                                               pf = FALSE)
rm(ver402_RR_T)

ver6402_RR_T <- readRDS(file = file.path(dir.output, "Ver6402sstpoint_14Sept23.RDS"))

RetroList_ver6402RR_T <- retrospective.function(outputList = ver6402_RR_T,
                                                testYears = testYears,
                                                testDays = testDays,
                                                CAN_hist = CAN_hist_new,
                                                pf = FALSE)
rm(ver6402_RR_T)

ver6403_RR_T <- readRDS(file = file.path(dir.output, "Ver6403sstpoint_14Sept23.RDS"))

RetroList_ver6403RR_T <- retrospective.function(outputList = ver6403_RR_T,
                                                testYears = testYears,
                                                testDays = testDays,
                                                CAN_hist = CAN_hist_new,
                                                pf = FALSE)
rm(ver6403_RR_T)

##### Calculations for retrospecitve testing ############################################
# Uses retrospective.function

# RetroList_ver102RR_T <- retrospective.function(outputList = ver102_RR_T,
#                                                testYears = testYears,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_new,
#                                                pf = FALSE)

# RetroList_ver102sstRR_T <- retrospective.function(outputList = ver1sst_RR_T,
#                                                testYears = testYears,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_new,
#                                                pf = FALSE)

# RetroList_ver102sstintRR_T <- retrospective.function(outputList = ver1sstint_RR_T,
#                                                   testYears = testYears,
#                                                   testDays = testDays,
#                                                   CAN_hist = CAN_hist_new,
#                                                   pf = FALSE)

# RetroList_ver602sstRR_T <- retrospective.function(outputList = ver602_RR_T,
#                                                   testYears = testYears,
#                                                   testDays = testDays,
#                                                   CAN_hist = CAN_hist_new,
#                                                   pf = FALSE)

RetroList_ver602sstpointRR_T <- retrospective.function(outputList = ver602point_RR_T,
                                                  testYears = testYears,
                                                  testDays = testDays,
                                                  CAN_hist = CAN_hist_new,
                                                  pf = FALSE)

# RetroList_ver602_null_RR_T <- retrospective.function(outputList = ver602null_RR_T,
#                                                        testYears = testYears,
#                                                        testDays = testDays,
#                                                        CAN_hist = CAN_hist_new,
#                                                        pf = FALSE)


# RetroList_ver402RR_T <- retrospective.function(outputList = ver402_RR_T,
#                                                testYears = testYears,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_new,
#                                                pf = FALSE)
# rm(ver402_RR_T)

# RetroList_ver6402RR_T <- retrospective.function(outputList = ver6402_RR_T,
#                                                testYears = testYears,
#                                                testDays = testDays,
#                                                CAN_hist = CAN_hist_new,
#                                                pf = FALSE)

# RetroList_ver6403RR_T <- retrospective.function(outputList = ver6403_RR_T,
#                                                 testYears = testYears,
#                                                 testDays = testDays,
#                                                 CAN_hist = CAN_hist_new,
#                                                 pf = FALSE)

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

rmseDF_ver102RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver102RR_T$RMSE_by_day_vect,
                                "Version" = "1.0.2",
                                "Method" = "Regression",
                                "Hyp" = "Null")

rmseDF_ver102sstRR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver102sstRR_T$RMSE_by_day_vect,
                                "Version" = "1.0.2.sst",
                                "Method" = "Regression",
                                "Hyp" = "SST")

rmseDF_ver602sstRR_T <- data.frame("Day" = testDays,
                                   "RMSE" = RetroList_ver602sstRR_T$RMSE_by_day_vect,
                                   "Version" = "6.0.2.sst",
                                   "Method" = "Proportion",
                                   "Hyp" = "SST")

rmseDF_ver602nullRR_T <- data.frame("Day" = testDays,
                                    "RMSE" = RetroList_ver602_null_RR_T$RMSE_by_day_vect,
                                    "Version" = "6.0.2 null",
                                    "Method" = "Proportion",
                                    "Hyp" = "Null")

rmseDF_ver402RR_T <- data.frame("Day" = testDays,
                                "RMSE" = RetroList_ver402RR_T$RMSE_by_day_vect,
                                "Version" = "4.0.2",
                                "Method" = "Normal Dist",
                                "Hyp" = "Null")

# rmseDF_ver6402RR_T <- data.frame("Day" = testDays,
#                                 "RMSE" = RetroList_ver6402RR_T$RMSE_by_day_vect,
#                                 "Version" = "6.4.0.2",
#                                 "Method" = "Normal Dist",
#                                 "Hyp" = "SST")

rmseDF_ver6403RR_T <- data.frame("Day" = testDays,
                                 "RMSE" = RetroList_ver6403RR_T$RMSE_by_day_vect,
                                 "Version" = "4.0.2 SST",
                                 "Method" = "Normal Dist",
                                 "Hyp" = "SST")
# 
rmseDF_pf_new <- data.frame("Day" = testDays,
                            "RMSE" = RetroList_pf_new$RMSE_by_day_vect,
                            "Version" = "PF",
                            "Method" = "PF")

# rmseDF_pf_old <- data.frame("Day" = testDays,
#                             "RMSE" = RetroList_pf_old$RMSE_by_day_vect,
#                             "Version" = "PF Eagle+Harv")

# Combine data frames into one data.frame for plotting
full_rmseDF <- rbind(
  rmseDF_ver102RR_T,
  rmseDF_ver102sstRR_T,
  rmseDF_ver602sstRR_T,
  rmseDF_ver602nullRR_T,
  rmseDF_ver402RR_T,
  # rmseDF_ver6402RR_T,
  rmseDF_ver6403RR_T
)

# Make version a factor
# full_rmseDF$Version <- factor(full_rmseDF$Version, levels = c("1.0",
#                                                               "2.0",
#                                                               "2c.0",
#                                                               "1.1",
#                                                               "1.0.1",
#                                                               "1.0.1.if",
#                                                               "1.0.2",
#                                                               "4.0.2",
#                                                               "5.0.2"
#                                                               # "PF"
# ))

# Vector of dates for labeling x axis
 x.date <- c(
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
 )

 # Table for papers
# RMSE_table <-as.data.frame(pivot_wider(full_rmseDF,names_from = Day, values_from = RMSE))

# RMSE_table_round <- format(RMSE_table, digits = 2)

# write.table(x = RMSE_table_round, file = file.path(dir.output,"RMSE Table 13Jul23"), row.names = F)

# write.csv(x = RMSE_table_round, file = file.path(dir.output,"RMSE Table 17Jul23 2"), row.names = FALSE)
ggplot(full_rmseDF, aes(x = factor(Day), 
                        y = factor(Version), 
                        fill = RMSE))+
  geom_tile()+
  geom_text(aes(label = round(RMSE)))+
  scale_fill_distiller(palette = "Spectral")+
  labs(x = "Day of Year",
       y = "Model Versions")+
  scale_x_discrete(breaks = c(testDays), labels = x.date
  )+
  theme(text = element_text(size = 14)
  )
# Save to PDF 
# **Uncomment out to save plots after this line to PDF
# pdf(file = file.path(dir.output,"retro_plots_29Sept22.pdf"))
# Add PF for plotting
full_rmseDF$PF <- rmseDF_pf_new$RMSE


# png(filename = file.path(dir.figs,"RMSE Plot 31July23.png"), width = 1000, height = 800)
# RMSE plot
(rmse.plot <- ggplot(full_rmseDF, aes(x = Day, y = RMSE,
                        col = Hyp
))+
  geom_point()+
  geom_line()+
  coord_cartesian(ylim = c(
    min(full_rmseDF$RMSE)-500,
    max(full_rmseDF$RMSE)+500
    # 10000,15000
  
  ))+
    xlab("")+
  scale_color_colorblind(name = "",
                         labels = c("Preseason Forecast", "Null", "SST","SST 6403"))+

# labels = c("Ver 1", "Ver 1-PF", "Ver 2", "Ver2-PF"))+
guides(alpha = "none")+
  scale_x_continuous(breaks = c(testDays), labels = x.date)+
  geom_hline(aes(yintercept = PF,
                 color = factor(PF)
                 ),size = 1,linetype = 2, show.legend = T)+
  facet_wrap(~Method)+
  # scale_color_discrete(labels = c("Preseason Forecast", "Null", "SST","SST 6403"),
  #                      name = "")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 0.5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90),
        
        axis.text.x = element_blank(),)+
    theme(plot.margin = unit(c(1,1,0,1), 'lines')))


# dev.off()

# MAPE ####################################################################

# Create MAPE dataframe

mapeDF_ver102RR_T <- data.frame("Day" = testDays,
                                "MAPE" = RetroList_ver102RR_T$MAPE_vect,
                                "Version" = "1.0.2",
                                "Method" = "Regression",
                                "Hyp" = "Null")

mapeDF_ver102sstRR_T <- data.frame("Day" = testDays,
                                "MAPE" = RetroList_ver102sstRR_T$MAPE_vect,
                                "Version" = "1.0.2.sst",
                                "Method" = "Regression",
                                "Hyp" = "SST")

mapeDF_ver602sstRR_T <- data.frame("Day" = testDays,
                                   "MAPE" = RetroList_ver602sstRR_T$MAPE_vect,
                                   "Version" = "6.0.2 SST",
                                   "Method" = "Proportion",
                                   "Hyp" = "SST")

mapeDF_ver602nullRR_T <- data.frame("Day" = testDays,
                                    "MAPE" = RetroList_ver602_null_RR_T$MAPE_vect,
                                    "Version" = "6.0.2 null",
                                    "Method" = "Proportion",
                                    "Hyp" = "Null")

mapeDF_ver402RR <- data.frame("Day" = testDays,
                              "MAPE" = RetroList_ver402RR_T$MAPE_vect,
                              "Version" = "4.0.2",
                              "Method" = "Normal Dist",
                              "Hyp" = "Null")

# mapeDF_ver6402RR <- data.frame("Day" = testDays,
#                               "MAPE" = RetroList_ver6402RR_T$MAPE_vect,
#                               "Version" = "6.4.0.2",
#                               "Method" = "Normal Dist",
#                               "Hyp" = "SST")

mapeDF_ver6403RR <- data.frame("Day" = testDays,
                               "MAPE" = RetroList_ver6403RR_T$MAPE_vect,
                               "Version" = "4.0.2 SST",
                               "Method" = "Normal Dist",
                               "Hyp" = "SST")

mapeDF_PF_new <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_pf_new$MAPE_vect,
                            "Version" = "PF",
                            "Method" = "PF")

# Combine into one data-frame
full_MAPE.df <- rbind(
  mapeDF_ver102RR_T,
  mapeDF_ver102sstRR_T,
  # mapeDF_ver102sstintRR_T,
  mapeDF_ver602sstRR_T,
  # mapeDF_ver602pointsstRR_T,
  mapeDF_ver602nullRR_T,
  mapeDF_ver402RR,
  # mapeDF_ver6402RR,
  mapeDF_ver6403RR
  # mapeDF_PF_new
)

ver.names <- unique(full_MAPE.df$Version)

# # Make version a factor
full_MAPE.df$Version <- factor(full_MAPE.df$Version, levels = c(
                                                                ver.names
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
  scale_x_discrete(breaks = c(testDays), labels = x.date)+
  labs(x = "Day of Year",
       y = "Model Versions")+
  theme(text = element_text(size = 14))


# png(file = file.path(dir.figs,"MAPE Plot 31July23.png"),width = 1000, height = 800)
# MAPE plot
(mape.plot <-ggplot(full_MAPE.df, aes(x = Day,
                         y = MAPE*100,
                         col = Hyp
                         # shape = Version
))+
  geom_point(size = 1.5)+
  geom_line(aes())+
  geom_hline(aes(yintercept = PF*100, color = factor(PF)),size = 1,linetype = 2, show.legend = T)+
  coord_cartesian(ylim = c(min(full_MAPE.df$MAPE*100)-2,max(full_MAPE.df$MAPE*100)+2))+
  ylab('MAPE (%)')+
  # xlab("")+
  scale_x_continuous(breaks = c(testDays), labels = x.date)+
  facet_wrap(~Method)+
  scale_color_colorblind(name = "",
                         labels = c("Preseason Forecast", "Null", "SST","SST 6403"))+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = .5, 
                                          color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(angle = 90),
        # axis.text.x = element_blank(),
        strip.text = element_blank()
        )+
    theme(plot.margin = unit(c(-1.25,1,1,1), 'lines'))
  )

# dev.off()

png(file = file.path(dir.figs,"MAPE_RMSE Plot 20Sept23.png"),width = 1000, height = 800)

# Combine mape and rmse plots into one plot
ggarrange(rmse.plot,mape.plot, 
          ncol = 1, common.legend =  T)

dev.off()

# Percent Error Box Plots ##############################################################

PE_ver102 <- as.data.frame(RetroList_ver102RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver102$version <- "1.0.2 Null"


PE_ver102sst <- as.data.frame(RetroList_ver102sstRR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver102sst$version <- "1.0.2 SST"

PE_ver602 <- as.data.frame(RetroList_ver602sstRR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver602$version <- "6.0.2 SST"

PE_ver602null <- as.data.frame(RetroList_ver602_null_RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver602null$version <- "6.0.2 Null"

PE_ver402 <- as.data.frame(RetroList_ver402RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver402$version <- "4.0.2 Null"

# PE_ver6402 <- as.data.frame(RetroList_ver6402RR_T$PE_mat) %>%
#   pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
#   as.data.frame()
# PE_ver6402$version <- "6.4.0.2 SST"

PE_ver6403 <- as.data.frame(RetroList_ver6403RR_T$PE_mat) %>%
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver6403$version <- "4.0.2 SST"

PE_pf_new <- as.data.frame(RetroList_pf_new$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_pf_new$version <- "PF New"

peDF_total <- rbind(
  PE_ver102,
  PE_ver102sst,
  PE_ver602null,
  PE_ver602,
  PE_ver402,
  # PE_ver6402,
  PE_ver6403
  
)

# For plotting PF point shape 

newDF <- left_join(x = pf_hist, y = CAN_hist_new)
# newDF <- left_join(x = newDF, y = CAN_hist_old, "Year")

# newDF$PF_new <- (newDF$Weighted_Forecast-newDF$can.mean.x)/newDF$can.mean.x

newDF$Year <- as.factor(newDF$Year)

newDF$PF <- (newDF$mean-newDF$can.mean)/newDF$can.mean



peDF_total<-left_join(peDF_total, newDF[,c(2,4)])

peDF_total$version <- factor(peDF_total$version, levels = c("1.0.2 Null",
                                                            "1.0.2 SST",
                                                            "6.0.2 Null",
                                                            "6.0.2 SST",
                                                            "4.0.2 Null",
                                                            "4.0.2 SST"))

# peDF_total <- rbind(peDF_total, newDF[,c(2,4,5)])
# png(file = file.path(dir.figs,"PE Plot Chapt2 20Sept23.png"),width = 1000, height = 800)
# long_df <- newDF %>% pivot_longer(cols = c(PF_new,PF_old)) %>% as.data.frame()
ggplot(peDF_total, aes(x = Year, y = PE*100, 
                       fill =version))+
  guides(alpha = "none")+
  geom_boxplot(alpha = .65, position = position_dodge(width = 1.01)) + 
  labs(fill = "Model Version",
       y = "Percent Error",
       x = "Model",
       color = "")+
  geom_point(
    aes( x = Year, y = PF*100, fill = "PF"),
    shape = 21,
    position = position_nudge(x = -.55),
    size = 3)+
  scale_fill_colorblind(name = "")+
  scale_color_colorblind(name = "")+
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

cv.df.102 <- cv.func(mod = ver102_RR_T)

cv.df.102$model <- "1.0.2 Null"

cv.df.102sst <- cv.func(mod = ver102sst_RR_T)

cv.df.102sst$model <- "1.0.2 SST"

cv.df.602sst <- cv.func(mod = ver602_RR_T)

cv.df.602sst$model <- "6.0.2 SST"

cv.df.602null <- cv.func(mod = ver602null_RR_T)

cv.df.602null$model <- "6.0.2 Null"

ver402_RR_T <- readRDS(file=file.path(dir.output,"Ver402_rr_trunc_FINAL_29May23.RDS"))

cv.df.402 <- cv.func(mod = ver402_RR_T)

cv.df.402$model <- "4.0.2 Null"

rm(ver402_RR_T)

# ver6402_RR_T <- readRDS(file = file.path(dir.output, "Ver6402sstpoint_14Sept23.RDS"))
# 
# cv.df.6402 <- cv.func(mod = ver6402_RR_T)
# 
# cv.df.6402$model <- "4.0.2 SST"

# rm(ver6402_RR_T)

ver6403_RR_T <- readRDS(file = file.path(dir.output, "Ver6403sstpoint_14Sept23.RDS"))

cv.df.6403 <- cv.func(mod = ver6403_RR_T)

cv.df.6403$model <- "4.0.2 SST"

rm(ver6403_RR_T)

cv.total <- rbind(
  cv.df.102,
  cv.df.102sst,
  cv.df.602sst,
  cv.df.602null,
  cv.df.402,
  # cv.df.6402,
  cv.df.6403
  )

# str(cv.df.102)

cv.total$cv <- as.numeric(cv.total$cv)
cv.total$year <- as.numeric(cv.total$year)
cv.total$day <- as.numeric(cv.total$day)
cv.total$model <- factor(cv.total$model, c("1.0.2 Null",
                                           "1.0.2 SST",
                                           "6.0.2 Null",
                                           "6.0.2 SST",
                                           "4.0.2 Null",
                                           "4.0.2 SST"))


# png(file = file.path(dir.figs,"CV Plot Chapt 2 20Sept23.png"),width = 1300, height = 800)

cust.col <- c("black",
              "#E69F00",
              "#56B4E9",
              "#009E73",
              "#F0E442",
              "#0072B2",
              "#D55E00",
              "#CC79A7",
              "darkgray")

ggplot(cv.total, aes(x = day, y = cv*100, 
                     col = model))+
  geom_line(size = .75)+
  geom_point(size = 2)+
  labs(x = "Day",
       y = "CV (%)")+
  # scale_color_colorblind()+
  scale_color_manual(values = cust.col,
                     name = "")+
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

# Plots of median run size prediction over each year

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
rm(ver402_RR_T)
gc()

ver502_RR_T <- readRDS(file = file.path(dir.output,"Ver502_rr_trunc_final_31Julu23.RDS"))

runTrend_ver502 <- runsize_trend_func(ver502_RR_T,
                                      testYears = testYears,
                                      testDays = testDays ,
                                      CAN_hist = CAN_hist_new,PF = pf_hist)
rm(ver502_RR_T)


comb.df <- rbind(
  # runTrend_ver1,
  # runTrend_ver11,
  # runTrend_ver101,
  # runTrend_ver101if,
  runTrend_ver102,
  runTrend_ver402,
  runTrend_ver502
)


# Plot the error
comb.df$error <- comb.df$RunSize - comb.df$can.mean

ggplot(comb.df, aes(x = Day, 
                    y = error/1000,
                    color = Version
                    
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
  scale_color_colorblind()+
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
  df <-interval_function(mod = ver502_RR_T, 
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


# Sigma from PSS regression################

ggplot(final.df, aes(x = day, y = median))+
  # geom_ribbon(aes(ymin = median-0.003, ymax = median+0.003, fill = "Median")) +
  geom_line(aes(x = day, y = median))+
  # geom_point(aes(
  #   col = param),
  #   size = 2, 
  #   show.legend = T)+
  geom_ribbon(aes(ymin=low95, ymax=up95, fill = "HDI 95"), alpha = .2,
              show.legend = T) +
  geom_ribbon(aes(ymin=low50,
                  ymax=up50, fill = "HDI 50"), alpha = .5) +
  facet_wrap(~Year)+
  labs(x = "Day",
       y = bquote(sigma^2),
       fill = ""
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
  theme_light(base_size = 15)+
  theme(axis.title    = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90))


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
# List to hold plots
plots_list <- list()

# Loop over outputs to generate plots
for (i in 1:length(ver102_RR_T)){
  
  
  plots_list[[i]] <- density.func(outputList = ver1sstint_RR_T[[i]],
                              CAN_hist = CAN_hist_new)
}

# Vector of names form plot list
test<-names(ver102_RR_T)
names(plots_list)<-test


figure<-ggarrange(plots_list$`2016_238`$Dens.plot + 
                    labs(x = "", y ="")+
                    theme_light(),
                  
                  plots_list$`2017_238`$Dens.plot + 
                    labs(x = "", y ="")+
                    theme_light(),
                  
                  plots_list$`2018_238`$Dens.plot + 
                    labs(x = "", y ="")+
                    theme_light(),
                  
                  plots_list$`2019_238`$Dens.plot + 
                    labs(x = "", y ="")+
                    theme_light(),
                  ncol = 4, 
                  nrow = 1, 
                  common.legend = TRUE,
                  
                  label.x = .4,
                  label.y = 1,
                  labels = c("Aug 26 2016",
                             "Aug 26 2017",
                             "Aug 26 2018",
                             "Aug 26 2019"),
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
