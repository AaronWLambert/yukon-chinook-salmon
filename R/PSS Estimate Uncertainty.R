###############################
# Aaron Lambert
# Yukon Chinook Inseason Forecast
# Analysis of PSS passage estimate uncertainty
# 2/15/2023

library(tidyverse)

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
# NOTE!!!! Use "Dataframe preprocess.R" to get  most recent version.
#  If not using Dataframe preprocess.R, uncomment out the following to load relevant PSS historical data
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage Final 2022.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
# GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2007 - 2021)
pf_hist <- readRDS(file.path(dir.data,"temp_PF_insamp_12Dec22.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs) and days not covered by GSI
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 27April22.RDS"))

# Eagle passage
Eagle_hist <- readRDS(file = file.path(dir.data, "Eagle passage 23Nov22.RDS"))

# Variance for daily PSS Chinook Passage
var_pss <- read.csv(file = file.path(dir.data, "Chinook_Variances_Pilot_Station_1995_2022.csv"))

# Retain relevant columns
pssSD <- var_pss[c(4,6,7,8)]

# Break out days of year and year
pssSD  <- pssSD %>% mutate(Year = lubridate::year(Observation.Date),
                              # month = lubridate::month(),
                              # day = lubridate::day(valid),
                              Day = lubridate::yday(Observation.Date))

# Calculate SD
pssSD$sd <- sqrt(pssSD$Variance)

# Calculate CV
pssSD$CV <- (pssSD$sd/pssSD$Fish.Count)*100

# As a factor
pssSD$Size.Category <- as.factor(pssSD$Size.Category)

# Outlier function
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# Label outliers for plotting
pssSD <- na.omit(pssSD) %>%
  group_by(Year) %>%
  mutate(outlier = ifelse(findoutlier(CV), Day, NA))

# Plot Cv by Day
ggplot(pssSD, aes(x = Day, y = CV, color = Size.Category))+
  geom_point()+
  # geom_boxplot()+
  labs(x = "Day of Year",
       y = "CV (%)")+
  facet_wrap(~Year)

# Cv by year
ggplot(pssSD, aes(x = as.factor(Year), y = CV, color = Size.Category))+
  geom_boxplot(na.rm = TRUE ,
               # aes(colour = Size.Category),
               outlier.alpha = .9)+
  geom_text(aes(label = outlier),
            na.rm = TRUE, 
            # nudge_y = 50,
            nudge_x = .5,
            size = 3, 
            alpha = .7)+
  geom_hline(aes(yintercept = 100, 
                 color = "100 %"),
             size = 1,
             linetype = 2)+
  labs(x = "Year",
       y = "CV (%)",
       title = "CV of PSS Chinook passage estimates",
       color = c("Chinook Size"))+
  theme_classic()


na.omit(pssSD[pssSD$CV >= 2,])

pssSD[pssSD$Day>=200,]

sqrt(sum(pssSD$Variance[pssSD$Year == 2019]))

pssVar <- pssSD %>% group_by(Year) %>% 
  summarise(Var = sum(Variance),
            totalChinook = sum(Fish.Count))

pssVar$Sd <- sqrt(pssVar$Var)

pssVar$CV <- pssVar$Sd/pssVar$totalChinook

# CV by year
ggplot(pssVar, aes(x = factor(Year), y = CV*100))+
  geom_point()+
  labs(y = "CV (%)",
       x = "Year",
       title = "CV by Year")+
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 18))

head(pssVar)

# Catagroze cv as high or low
pssVar$HIGH_CV <- ifelse(pssVar$CV>.2, "High", "Low")


head(CAN_hist)

pssVar_DF <- left_join(pssVar, CAN_hist)

# Plot EOS vs EOS PSS passage to identify years with high obs error
ggplot(pssVar_DF, aes(x = totalChinook/1000, y = can.mean/1000))+
  geom_point(aes(col = HIGH_CV, size = CV*100))+
  geom_smooth(method = "lm")+
  geom_text(aes(label = Year),nudge_y = 5)+
  labs(x = "EOS PSS Chinook Passage (1000's)",
       y = "EOS Canadian-origin Chinook (1000's)",
       col = "CV")+
  scale_color_discrete(labels = c("CV > 20%", "CV < 20% "))+
  scale_size(name = "CV")



PSS_sd <- pssSD %>% group_by(Year,Day) %>% summarize(Var = sum(Variance)) %>% as.data.frame()

saveRDS(PSS_sd, file = file.path(dir.data,"PSS SD 1995_2021.RDS"))
