#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.14.22
# Purpose: This is a working script to do retrospective analysis with
# 
#
#   1) Read in model outputs
#   2) Generate MAPE, PE, MRSE, and other stats
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

library(tidyverse)
library(tidybayes)

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
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
pf_hist <- readRDS(file.path(dir.data,"inv_var_weighted_forcast_v3_Jan282022.RDS"))
# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
# read in  model outputs
outputList <- readRDS(file = file.path(dir.output,"OutPut_ver2_23mar22.RDS"))



# For runs of multiple days with retrospective testing after
testDays <- seq(from = 153, to = 213, by = 5)
# testDays <- 153
testYears <- 2013:2019
# Calculations for retrospecitve testing #############################################

# Create matrix of observed and median outputs
# 
# Create matrix with dimensions of days by years that are in the outputList from model run
median_mat <- matrix(nrow = length(testDays),ncol = length(testYears))

# Give names to matrix
colnames(median_mat) <- c(testYears)
rownames(median_mat) <- c(testDays)

# Make sure it works...
dim(median_mat)
str(median_mat)

# Extract medians from posterior RunSize prediction and place in vector for matrix pop.
# populate matrix with medians from posterior RunSize prediction
counter <- 1
for (y in 1:length(testYears)) {
  for (d in 1:length(testDays)) {
    
    median_mat[d,y]<- median(outputList[[counter]]$pars$RunSize)
    
    counter <- counter+1
  }
}
median_mat

# Matrix with realized reconstructed counts by year and day 

realized_mat <- matrix(nrow = length(testDays),ncol = length(testYears))

# Give names to matrix
colnames(realized_mat) <- c(testYears)
rownames(realized_mat) <- c(testDays)

# Make sure it works...
dim(realized_mat)

# Extract values for realized runsize
realized_vect <- CAN_hist$can.mean[CAN_hist$Year >= 2013]

counter <- 1
# Populate the realized matrix
for (y in 1:length(testYears)) {
  for (d in 1:length(testDays)) {
    
    realized_mat[d,y]<- realized_vect[counter]
    
  }
  counter <- counter+1
}
realized_mat
# Percent Error (PE) matrix
PE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))

# Give names to matrix
colnames(PE_mat) <- c(testYears)
rownames(PE_mat) <- c(testDays)

# Populate the PE matrix
for (y in 1:length(testYears)) {
  for (d in 1:length(testDays)) {
    
    PE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
    
  }
}
# Percent Error (APE) matrix
APE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))

# Give names to matrix
colnames(APE_mat) <- c(testYears)
rownames(APE_mat) <- c(testDays)

# Populate the absolute PE matrix
for (y in 1:length(testYears)) {
  for (d in 1:length(testDays)) {
    
    APE_mat[d,y]<- abs(realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
    
  }
  
}


# Calculate MAPE by Days
MAPE_day_vect <- vector(length = length(PE_mat[,1]))
names(MAPE_day_vect) <- testDays

for (d in 1:length(MAPE_day_vect)) {
  MAPE_day_vect[d] <- mean(APE_mat[d,])
}

# Root Mean Square error for days across years
# Square Error (SE) matrix
SE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))

# Give names to matrix
colnames(SE_mat) <- c(testYears)
rownames(SE_mat) <- c(testDays)

# Populate the SE matrix
for (y in 1:length(testYears)) {
  for (d in 1:length(testDays)) {
    
    SE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])^2
    
  }
  
}

# Calculate RMSE by days
RMSE_day_vect <- vector(length = length(SE_mat[,1]))
names(RMSE_day_vect) <- testDays

for (d in 1:length(RMSE_day_vect)) {
  RMSE_day_vect[d] <- sqrt(mean(SE_mat[d,]))
}

# calculate RMSE by year
RMSE_year_vect <- vector(length = length(SE_mat[1,]))
names(RMSE_year_vect) <- as.character(testYears)

for (y in 1:length(RMSE_year_vect)) {
  RMSE_year_vect[y] <- sqrt(mean(SE_mat[,y]))
}

# How often does the interval capture the true return size? ######
# 95% Interval
true95Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
true80Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
true50Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
# Give names to matrix
colnames(true95Mat) <- c(testYears)
rownames(true95Mat) <- c(testDays)

colnames(true80Mat) <- c(testYears)
rownames(true80Mat) <- c(testDays)

colnames(true50Mat) <- c(testYears)
rownames(true50Mat) <- c(testDays)

count <- 1

for (y in 1:length(testYears)) {
  for (d in 1:length(testDays)) {
    
    # 95%
    test95<-quantile(outputList[[count]]$pars$RunSize, probs = c(.025,.975))
    
    true95Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test95[2] & 
                               CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test95[1],
                             1,
                             0)
    # 80%
    test80<-quantile(outputList[[count]]$pars$RunSize, probs = c(.1,.9))
    
    true80Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test80[2] & 
                               CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test80[1],
                             1,
                             0)
    # 50%
    test50<-quantile(outputList[[count]]$pars$RunSize, probs = c(.25,.75))
    
    true50Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test50[2] & 
                               CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test50[1],
                             1,
                             0)
    count <- count+1
  }
}
true95Mat
true80Mat
true50Mat
# Plot of RMSE

# First by days
rmseDF <- data.frame("Day" = testDays, "RMSE" = RMSE_day_vect)

dev.off()
ggplot(rmseDF, aes(x = Day, y = RMSE))+
  geom_col() + 
  theme_tidybayes()+
  labs(title = paste("RMSE by Day, across Years \n Model Version",outputList[[1]]$version))

# By years

rmseYearDF <- data.frame("Year" = testYears, "RMSE" = RMSE_year_vect)

ggplot(rmseYearDF, aes(x = Year, y = RMSE))+
  geom_col() + 
  theme_tidybayes()+
  labs(title = paste("RMSE by Year, across days \n Model Version",outputList[[1]]$model.version))

# Plot of MAPE

# Create MAPE dataframe
mapeDF <- data.frame("Day" = testDays, "MAPE" = MAPE_day_vect)

# Plot MAPE for days across years
ggplot(mapeDF, aes(x = Day, y = MAPE))+
  geom_col() + 
  theme_tidybayes()+
  labs(title = paste("MAPE \n Model Version",outputList[[1]]$version),
       y = "MAPE")

outputList$`2013_153`$version
# Plot for PE across year
PE_DF <- as.data.frame(PE_mat)

names(PE_DF)

peDF<-as.data.frame(pivot_longer(data = PE_DF,
                                 cols = starts_with("20")))

# Call boxplot for PE across years
ggplot(peDF, aes(x = name, y = value*100))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = 2)+
  labs(x = "Year",
       y = "Percent Error",
       title = paste("PE \n Model Version",outputList[[1]]$version))+
  theme_tidybayes()

# Side by side plot #######################################################

# RMSE by day
# rmseDF_ver1<-rmseDF
rmseDF_ver2 <- rmseDF

rmseDF_ver1$version <- "ver1"
rmseDF_ver2$version <- "ver2"



full_rmseDF <- rbind(rmseDF_ver1,rmseDF_ver2)

rmsePlot_day<-ggplot(full_rmseDF, aes(x = Day, y = RMSE, fill = version))+
  geom_col(position = "dodge") + 
  theme_tidybayes()+
  labs(title = paste("RMSE by Day, across Years"))+
  coord_cartesian(ylim = c(9000,12000))


# RMSE by year
# rmseYearDF_ver1<-rmseYearDF
rmseYearDF_ver2 <- rmseYearDF

rmseYearDF_ver1$version <- "ver1"
rmseYearDF_ver2$version <- "ver2"

full_rmseYearDF <- rbind(rmseYearDF_ver1,rmseYearDF_ver2)

rmsePlot_year <- ggplot(full_rmseYearDF, aes(x = Year, y = RMSE, fill = version))+
  geom_col(position = "dodge") + 
  theme_tidybayes()+
  labs(title = paste("RMSE by Years"))



# MAPE
# mapeDF_ver1 <- mapeDF
# mapeDF_ver2 <- mapeDF
mapeDF_ver1$version <- "ver1"
mapeDF_ver2$version <- "ver2"

mapeDF_long <- rbind(mapeDF_ver1,mapeDF_ver2)
# Plot MAPE for days across years
mapePlot<- ggplot(mapeDF_long, aes(x = Day, y = MAPE, fill = version))+
  geom_col(position = "dodge") + 
  theme_tidybayes()+
  labs(title = paste("MPE \n Realized - median posterior"),
       y = "MPE")


# PE
# peDF_ver1 <- peDF
peDF_ver2 <- peDF

peDF_ver1$version <- "ver1"
peDF_ver2$version <- "ver2"

peDF_long <- rbind(peDF_ver1,peDF_ver2)
# Call boxplot for PE across years
pePlot<-ggplot(peDF_long, aes(x = name, y = value*100, fill = version))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = 2)+
  labs(x = "Year",
       y = "Percent Error",
       title = paste("PE \n Model Version",outputList[[1]]$version))+
  theme_tidybayes()


pdf(file = file.path(dir.figs,"Compare Models first try.pdf"))

mapePlot
pePlot
rmsePlot_day
rmsePlot_year
dev.off()
