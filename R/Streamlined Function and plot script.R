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

# Define Workflow Paths ============================================
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")

dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

# Call in functions from model_run_funciton.r
source(file = file.path(dir.R,"model_run_function.r"))
source(file = file.path(dir.R,"Retro Function.R"))
######### Import Data ###############
# Historical Canadian EOS reconstructed run
# CAN_hist <- readRDS(file.path(dir.data,"/EOS Reconstructed Canadian.RDS"))
# CAN_hist <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
# This is the reconstructed data from Curry for old reconstructed modeling procedure
CAN_hist <- readRDS(file.path(dir.data,"Can Origin Reconstructed 2Feb22.RDS"))
# PSS historical for Jun 1st = 1, to Aug 31 = 92
# Years 1995 - 2020
# PSScum_hist <- readRDS(file.path(dir.data,"/clean_counts.rds"))

# Read in PSS non-cummulative with date in date format
PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean <- readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
pf_hist <- readRDS(file.path(dir.data,"inv_var_weighted_forcast_v3_Jan282022.RDS"))
# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 27April22.RDS"))
# Control Section ######################
model.version <- "3.3"
# Range of years $$$$ to 2021
# myYear <- 2019

# Range of days 152 -243
# myDay <- 196

# MCMC Parameters
n.chains <- 4
n.iter <- 30000;#5e4
n.thin <- 2

# For runs of multiple days with retrospective testing after
testDays <- seq(from = 153, to = 213, by = 5)
# testDays <- 153
testYears <- c(2007:2010,2013:2021)

# Run the model for a single year and day of interest.
#  Asjust the inputs for the model in the control section above
# model.output <-InSeasonProjection(model.version = model.version,
#                                   myYear = 2013,
#                                   myDay = 165,
#                                   n.chains = n.chains,
#                                   n.iter = n.iter,
#                                   n.thin = n.thin,
#                                   CAN_hist = CAN_hist,
#                                   PSS_hist = PSS_hist,
#                                   # GSI_mean = GSI_mean,
#                                   pf_hist = pf_hist,
#                                   GSI_by_year = GSI_by_year)
# 
# model.output
# Loop to run model over a selected number of days and years ######################
# To run model, uncomment out the following section
# # List to store outputs
outputList<-list()
for(y in c(testYears)){
  for(d in c(testDays)){

    outputList[[paste("",y,"_",d, sep = "")]]<-InSeasonProjection(model.version = model.version,
                                                                  myYear = y,myDay = d,
                                                                  n.chains = n.chains,
                                                                  CAN_hist = CAN_hist,
                                                                  pf_hist = pf_hist,
                                                                  PSS_hist = PSS_hist,
                                                                  n.thin = n.thin,
                                                                  n.iter = n.iter,
                                                                  GSI_by_year = GSI_by_year)
    print(d)
  } #dloop
  print(y)
    } #yloop

outputList$`2007_153`$pars

# Save or read in outputlist #######################################

# Save output
# saveRDS(object = outputList, file = file.path(dir.output, "OutPut_ver21_oldCan_2007_2010_2013_2021_3.3.RDS"))
# outputlist_ver1 <- readRDS(file = file.path(dir.output,"OutPut_ver1_23mar22.RDS"))
# outputlist_ver2 <- readRDS(file = file.path(dir.output,"OutPut_ver2_23mar22.RDS"))
outputlist_ver1_oldCanHist <- readRDS(file = file.path(dir.output,"OutPut_ver1_oldCan_2007_2010_2013_2021_5apr22.RDS"))
outputlist_ver2_oldCanHist <- readRDS(file = file.path(dir.output,"OutPut_ver2_oldCan_2007_2010_2013_2021_5apr22.RDS"))
outputlist_ver2.1_oldCanHist <- readRDS(file = file.path(dir.output,"OutPut_ver21_oldCan_2007_2010_2013_2021_5apr22.RDS"))
outputlist_ver2.2_oldCanHist <- readRDS(file = file.path(dir.output,"OutPut_ver21_oldCan_2007_2010_2013_2021_5apr22.RDS"))
outputlist_ver3.0_oldCanHist <- readRDS(file = file.path(dir.output,"OutPut_ver30_oldCan_2007_2010_2013_2021_9may21.RDS"))
outputlist_ver3.2_oldCanHist <- readRDS(file = file.path(dir.output,"OutPut_ver32_oldCan_2007_2010_2013_2021_9may21.RDS"))
outputlist_ver3.3_oldCanHist <- readRDS(file = file.path(dir.output,"OutPut_ver33_oldCan_2007_2010_2013_2021_9may21.RDS"))
##### Calculations for retrospecitve testing #############################################

RetroList_ver1 <- retrospective.function(outputList = outputlist_ver1_oldCanHist,
                                         testYears = testYears,
                                         testDays = testDays,
                                         CAN_hist = CAN_hist,
                                         pf = FALSE)

RetroList_ver2 <- retrospective.function(outputList = outputlist_ver2_oldCanHist,
                                         testYears = testYears,
                                         testDays = testDays,
                                         CAN_hist = CAN_hist,
                                         pf = FALSE)

RetroList_ver2.1 <- retrospective.function(outputList = outputlist_ver2.1_oldCanHist,
                                            testYears = testYears,
                                            testDays = testDays,
                                            CAN_hist = CAN_hist,
                                            pf = FALSE)
RetroList_ver2.2 <- retrospective.function(outputList = outputlist_ver2.2_oldCanHist,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist,
                                           pf = FALSE)
RetroList_ver3.0 <- retrospective.function(outputList = outputlist_ver3.0_oldCanHist,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist,
                                           pf = FALSE)
RetroList_ver3.2 <- retrospective.function(outputList = outputlist_ver3.2_oldCanHist,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist,
                                           pf = FALSE)
RetroList_ver3.3 <- retrospective.function(outputList = outputlist_ver3.3_oldCanHist,
                                           testYears = testYears,
                                           testDays = testDays,
                                           CAN_hist = CAN_hist,
                                           pf = FALSE)

RetroList_pf <- retrospective.function(outputList = outputlist_ver2.1_oldCanHist,
                                            testYears = testYears,
                                            testDays = testDays,
                                            CAN_hist = CAN_hist,
                                            pf = TRUE)


# RMSE

# Extract RMSE by day into dataframe for each version and posterior of interest
rmseDF_ver1 <- data.frame("Day" = testDays, 
                          "RMSE" = RetroList_ver1$RMSE_by_day_vect,
                          "Version" = "ver1")
rmseDF_ver2 <- data.frame("Day" = testDays, 
                          "RMSE" = RetroList_ver2$RMSE_by_day_vect,
                          "Version" = "ver2")
rmseDF_ver2.1 <- data.frame("Day" = testDays, 
                             "RMSE" = RetroList_ver2.1$RMSE_by_day_vect,
                             "Version" = "ver2.1")
rmseDF_ver2.2 <- data.frame("Day" = testDays, 
                            "RMSE" = RetroList_ver2.2$RMSE_by_day_vect,
                            "Version" = "ver2.2")
rmseDF_ver3.0 <- data.frame("Day" = testDays, 
                            "RMSE" = RetroList_ver3.0$RMSE_by_day_vect,
                            "Version" = "ver3.0")
rmseDF_ver3.2 <- data.frame("Day" = testDays, 
                            "RMSE" = RetroList_ver3.2$RMSE_by_day_vect,
                            "Version" = "ver3.2")
rmseDF_ver3.3 <- data.frame("Day" = testDays, 
                            "RMSE" = RetroList_ver3.3$RMSE_by_day_vect,
                            "Version" = "ver3.3")
rmseDF_pf <- data.frame("Day" = testDays,
                             "RMSE" = RetroList_pf$RMSE_by_day_vect,
                             "Version" = "pf")

full_rmseDF <- rbind(rmseDF_ver1,
                     rmseDF_ver2,
                     rmseDF_ver2.1,
                     rmseDF_ver2.2,
                     rmseDF_ver3.0,
                     rmseDF_ver3.2,
                     rmseDF_ver3.3,
                     rmseDF_pf)

full_rmseDF$Version <- as.factor(full_rmseDF$Version)
str(full_rmseDF)

ggplot(full_rmseDF, aes(x = Day, y = RMSE, fill = Version))+
  geom_col(position = "dodge",width = 3) + 
  labs(title = "RMSE by Day, across Years",
                     fill = "")+
  coord_cartesian(ylim = c(9000,20000))+
  # scale_fill_manual(values = wes_palette("IsleofDogs1"))+
  # labels = c("Ver 1", "Ver 1-PF", "Ver 2", "Ver2-PF"))+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank())


# By years RMSE

rmseYearDF_ver1 <- data.frame("Year" = testYears, 
                              "RMSE" = RetroList_ver1$RMSE_by_year_vect,
                              "Version" = "Ver 1")
rmseYearDF_ver2 <- data.frame("Year" = testYears, 
                              "RMSE" = RetroList_ver2$RMSE_by_year_vect,
                              "Version" = "Ver 2")
rmseYearDF_ver2.1 <- data.frame("Year" = testYears, 
                                 "RMSE" = RetroList_ver2.1$RMSE_by_year_vect,
                                 "Version" = "Ver 2.1")
rmseYearDF_ver2.2 <- data.frame("Year" = testYears, 
                                "RMSE" = RetroList_ver2.2$RMSE_by_year_vect,
                                "Version" = "Ver 2.2")
rmseYearDF_ver3.0 <- data.frame("Year" = testYears, 
                                "RMSE" = RetroList_ver3.0$RMSE_by_year_vect,
                                "Version" = "Ver 3.0")
rmseYearDF_ver3.2 <- data.frame("Year" = testYears, 
                                "RMSE" = RetroList_ver3.2$RMSE_by_year_vect,
                                "Version" = "Ver 3.2")
rmseYearDF_ver3.3 <- data.frame("Year" = testYears, 
                                "RMSE" = RetroList_ver3.3$RMSE_by_year_vect,
                                "Version" = "Ver 3.3")
rmseYearDF_PF <- data.frame("Year" = testYears, 
                                 "RMSE" = RetroList_pf$RMSE_by_year_vect,
                                 "Version" = "PF")

full_rmseYearDF <- rbind(rmseYearDF_ver1,
                     rmseYearDF_ver2,
                     rmseYearDF_ver2.1,
                     rmseYearDF_ver2.2,
                     rmseYearDF_ver3.0,
                     rmseYearDF_ver3.2,
                     rmseYearDF_ver3.3,
                     rmseYearDF_PF)

ggplot(full_rmseYearDF, aes(x = Year, y = RMSE, fill = Version))+
  geom_col(position = "dodge",width = .8) + 
  labs(title = "RMSE by Year, across Days",
       fill = "")+
  coord_cartesian(ylim = c(1000,30000))+
  # scale_fill_manual(values = wes_palette("IsleofDogs1"))+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank())


# Create MAPE dataframe
mapeDF_ver1 <- data.frame("Day" = testDays, 
                          "MAPE" = RetroList_ver1$MAPE_vect, 
                          "Version" = "Ver 1")
mapeDF_ver2 <- data.frame("Day" = testDays, 
                          "MAPE" = RetroList_ver2$MAPE_vect,
                          "Version" = "Ver 2")
mapeDF_ver2.1 <- data.frame("Day" = testDays, 
                             "MAPE" = RetroList_ver2.1$MAPE_vect,
                             "Version" = "Ver 2.1")
mapeDF_ver2.2 <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_ver2.2$MAPE_vect,
                            "Version" = "Ver 2.2")
mapeDF_ver3.0 <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_ver3.0$MAPE_vect,
                            "Version" = "Ver 3.0")
mapeDF_ver3.2 <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_ver3.2$MAPE_vect,
                            "Version" = "Ver 3.2")
mapeDF_ver3.3 <- data.frame("Day" = testDays, 
                            "MAPE" = RetroList_ver3.3$MAPE_vect,
                            "Version" = "Ver 3.3")
mapeDF_PF <- data.frame("Day" = testDays, 
                             "MAPE" = RetroList_pf$MAPE_vect,
                             "Version" = "PF")

full_MAPE.df <- rbind(mapeDF_ver1,
                      mapeDF_ver2,
                      mapeDF_ver2.1,
                      mapeDF_ver2.2,
                      mapeDF_ver3.0,
                      mapeDF_ver3.2,
                      mapeDF_ver3.3,
                      mapeDF_PF)

ggplot(full_MAPE.df, aes(x = Day, y = MAPE*100, fill = Version))+
  geom_col(position = "dodge",width = 3) + 
  labs(title = "MAPE by day, across years",
       fill = "",
       y = "MAPE (%)")+
  # coord_cartesian(ylim = c(1000,30000))+
  # scale_fill_manual(values = wes_palette("IsleofDogs1"))+
                    # labels = c("Ver 1", "Ver 1-PF", "Ver 2", "Ver2-PF"))+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank())



# Plot for PE across year
PE_ver1 <- as.data.frame(RetroList_ver1$PE_mat) %>%
  pivot_longer(cols = starts_with("20"),names_to =  "Year" , values_to = "PE" ) %>%
  as.data.frame()
PE_ver1$version <- "Ver 1"

PE_ver2 <- as.data.frame(RetroList_ver2$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2$version <- "Ver 2"

PE_ver2.1 <- as.data.frame(RetroList_ver2.1$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2.1$version <- "Ver 2.1"

PE_ver2.2 <- as.data.frame(RetroList_ver2.2$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver2.2$version <- "Ver 2.2"

PE_ver3.0 <- as.data.frame(RetroList_ver3.0$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver3.0$version <- "Ver 3.0"

PE_ver3.2 <- as.data.frame(RetroList_ver3.2$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver3.2$version <- "Ver 3.2"

PE_ver3.3 <- as.data.frame(RetroList_ver3.3$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_ver3.3$version <- "Ver 3.3"

PE_pf <- as.data.frame(RetroList_pf$PE_mat) %>% 
  pivot_longer(cols = starts_with("20"), names_to = "Year", values_to = "PE") %>%
  as.data.frame()
PE_pf$version <- "PF"

peDF_total <- rbind(PE_ver1,
                    PE_ver2,
                    PE_ver2.1,
                    PE_ver2.2,
                    PE_ver3.0,
                    PE_ver3.2,
                    PE_ver3.3,
                    PE_pf)



ggplot(peDF_total, aes(x = Year, y = PE, fill = version))+
  geom_boxplot() + 
  geom_hline(yintercept  = 0, linetype = 2)+
  labs(fill = "")+
  scale_fill_colorblind()+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = .5, color = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        panel.background = element_rect(fill = "grey"))



##### Density Plots with plot function #####################################################

# List to hold plots
plots_old_can_ver3.2 <- list()
for (i in 1:length(outputlist_ver3.2_oldCanHist)) {
  

plots_old_can_ver3.2[[i]] <- outPlots(outputList = outputlist_ver3.2_oldCanHist[[i]],
                                   CAN_hist = CAN_hist,GSI = TRUE)
}

plots_old_can_ver3.2[[65]]$PredPlot # Vector of names form plot list
test<-names(outputlist_ver3.2_oldCanHist)


names(plots_old_can_ver3.2)<-test
plots_old_can_ver3.2$`2008_153`$DensPlot


figure<-ggarrange(plots_old_can_ver3.2$`2021_153`$DensPlot,
                  plots_old_can_ver3.2$`2021_178`$DensPlot,
                  plots_old_can_ver3.2$`2021_198`$DensPlot,
                  plots_old_can_ver3.2$`2019_153`$DensPlot,
                  plots_old_can_ver3.2$`2019_178`$DensPlot,
                  plots_old_can_ver3.2$`2019_198`$DensPlot,
                  plots_old_can_ver3.2$`2017_153`$DensPlot,
                  plots_old_can_ver3.2$`2017_178`$DensPlot,
                  plots_old_can_ver3.2$`2017_198`$DensPlot,
                  ncol = 3, 
                  nrow = 3, 
                  common.legend = TRUE,
                  label.x = .5,
                  label.y = 1,
                  labels = c("June 2 2021 ","June 27 2021","July 17 2021","June 2 2019","June 27 2019","July 17 2019",
                              "June 2 2017", "June 27 2017", "July 17 2017"),
                  font.label = list(size = 10, 
                            color = "black",
                            family = "serif"))



annotate_figure(figure, left = textGrob("Probability", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Thousands of Chinook Salmon", gp = gpar(cex = 1.3)))
outp          




plots_old_can_ver3.2[[]]