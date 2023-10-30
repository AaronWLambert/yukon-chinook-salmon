#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 01.22.23
# Purpose: Explore environmental covariates and timing
#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 01.22.23
# Purpose: Explore timing differences between Can-orig Chinook and total Chinook


library(bbmle)
library(ggpubr)
library(tidyverse)
library(rstan)
library(bayesplot)
library(ggthemes)
library(viridis)
library(tidybayes)
library(lubridate)

# Parralize for optimum model run time

rstan_options(auto_write = TRUE)

mc.cores = parallel::detectCores()
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

# Emmonak wind direction and speed
emmonak_wind <- read.csv(file = file.path(dir.data, "Emmonak ASOS Data 1995_2022 13Jan23.csv"))

# Midpoints
all.timing.logistic <- read.csv(file = file.path(dir.output,
                                                 "logistic curve parameters All Chinook 1995_2022.csv"))

em <- emmonak_wind[c(2,3,6,7)]

em$tmpf[em$tmpf == "M"] <- NA


emmonak_env  <- em %>% mutate(year = lubridate::year(valid),
                              month = lubridate::month(valid),
                              day = lubridate::day(valid),
                              jday = lubridate::yday(valid))
head(emmonak_env)

str(emmonak_env)

emmonak_env$tmpf <- as.numeric(emmonak_env$tmpf)
# emmonak_env$year <- as.factor(emmonak_env$year)
# emmonak_env$month <- as.factor(emmonak_env$month)
# emmonak_env$day <- as.factor(emmonak_env$day)
# Extract Temperature data and summarise the information


emmonak_env[emmonak_env$year!=1996,]

# Get max temp by month for each year #############################
maxTempDF <- emmonak_env[emmonak_env$year!=1996,] %>% 
  group_by(year,month) %>% 
  summarise(max_month_temp = max(tmpf, na.rm = TRUE))%>% 
  as.data.frame()


head(maxTempDF)

# Replace -INF with NA
maxTempDF$max_month_temp[maxTempDF$max_month_temp == -Inf] <- NA

# Import midpoint 
all.timing.logistic <- read.csv(file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))

midpoint <- all.timing.logistic[,c(2,5)]

hist((midpoint$mid))

midpoint$dev.mid <- mean(midpoint$mid) - midpoint$mid

tempDF_full <- left_join(maxTempDF, midpoint)


# Get max temp for each observed day and get the mean for each month

mean_maxTempDF <- emmonak_env %>% 
  group_by(year,jday,month) %>% 
  drop_na() %>% 
  summarise(max_day_temp = max(tmpf, na.rm = TRUE))%>% 
  as.data.frame()

# mean_maxTempDF$max_day_temp[mean_maxTempDF$max_day_temp == -Inf] <- NA

mean_maxTempDF2 <- mean_maxTempDF %>% 
  group_by(year, month) %>% 
  summarise(mean_month_temp = mean(max_day_temp, na.rm = TRUE))

# Join with full data set containing midpoint and max by month
tempDF_full <- left_join(tempDF_full, mean_maxTempDF2)

head(tempDF_full)

# Pivot longer for plotting
long_tempDF <- pivot_longer(tempDF_full, 
                            cols = c(mean_month_temp, max_month_temp), 
                            names_to = "method",
                            values_to = "temp")

long_tempDF$method <- as.factor(long_tempDF$method)

long_tempDF$month.name <- month.name[long_tempDF$month]

long_tempDF$month.name <- factor(long_tempDF$month.name, levels = month.name)

#
ggplot(long_tempDF, aes(x = temp, y = dev.temp, color = method))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~month.name, scales = "free")+
  scale_color_colorblind()+
  theme_classic()

# Create matrix of temperatures for correlation coeficient calculation

# Max temp by month
temp.mat.max <- matrix(ncol = length(unique(tempDF_full$month)), 
                   nrow = length(unique(tempDF_full$year)))

rownames(temp.mat.max) <- unique(tempDF_full$year)
colnames(temp.mat.max) <- unique(tempDF_full$month)

years <- unique(tempDF_full$year)

for (i in 1:length(unique(tempDF_full$month))) {
  for (j in 1:length(unique(tempDF_full$year))) {
    
  # i <- 1
  # j = 22
    
    y <- years[j]
    
    temp.mat.max[j,i] <- ifelse(length(tempDF_full$max_month_temp[tempDF_full$year == y &
                                            tempDF_full$month == i])==0,
                            NA,
                            tempDF_full$max_month_temp[tempDF_full$year == y & 
                                                   tempDF_full$month == i])
    
  }
}

# Mean temp matrix
temp.mat.mean <- matrix(ncol = length(unique(tempDF_full$month)), 
                       nrow = length(unique(tempDF_full$year)))

rownames(temp.mat.mean) <- unique(tempDF_full$year)
colnames(temp.mat.mean) <- unique(tempDF_full$month)

years <- unique(tempDF_full$year)

for (i in 1:length(unique(tempDF_full$month))) {
  for (j in 1:length(unique(tempDF_full$year))) {
    
    # i <- 1
    # j = 22
    
    y <- years[j]
    
    temp.mat.mean[j,i] <- ifelse(length(tempDF_full$mean_month_temp[tempDF_full$year == y &
                                                              tempDF_full$month == i])==0,
                                NA,
                                tempDF_full$mean_month_temp[tempDF_full$year == y & 
                                                       tempDF_full$month == i])
    
  }
}

mid.mat <- matrix(ncol = length(unique(tempDF_full$month)), 
                   nrow = length(unique(tempDF_full$year)))

# Create matrix of temperatures for correlation coeficient calculation
rownames(mid.mat) <- unique(tempDF_full$year)
colnames(mid.mat) <- unique(tempDF_full$month)

years <- unique(tempDF_full$year)

for (i in 1:length(unique(tempDF_full$month))) {
  for (j in 1:length(unique(tempDF_full$year))) {
    
    # i <- 1
    # j = 22
    
    y <- years[j]
    
    mid.mat[j,i] <- ifelse(length(tempDF_full$dev.mid[tempDF_full$year == y &
                                                          tempDF_full$month == i])==0,
                            NA,
                            tempDF_full$dev.mid[tempDF_full$year == y & 
                                                   tempDF_full$month == i])
    
  }
}

cor.mat <- matrix(nrow = 12,
                  ncol = 3)
cor.mat[,1] <- month.name
colnames(cor.mat) <- c("month", "r_max", "r_avg")

for (i in 1:length(unique(tempDF_full$month))) {
  


 cor.mat[i,2] <- round(as.numeric(unlist(cor.test(x = mid.mat[,i], y = temp.mat.max[,i])[4])),2)

 cor.mat[i,3] <- round(as.numeric(unlist(cor.test(x = mid.mat[,i], y = temp.mat.mean[,i])[4])),2)
 
}

cor.df <- as.data.frame(cor.mat)

cor.df$r_max <- as.numeric(cor.df$r_max)
cor.df$r_avg <- as.numeric(cor.df$r_avg)
cor.df$month <- factor(cor.df$month, levels = month.name)

cor.df <- pivot_longer(cor.df, 
                       names_to = "method", 
                       cols = c(r_max,r_avg), 
                       values_to = "temp")

ggplot(cor.df, aes(x = temp, y = month, color = method))+
  geom_vline(aes(xintercept = 0, color = "0"),
             size = 2,
             linetype = 2)+
  geom_point(size = 5)+
  labs(color = "Method",
       x = "Correlation Coeficient",
       y = "Month")+
  scale_color_colorblind(labels = c("r = 0","Avg Temp", "Max Temp"))+
  theme(text = element_text(size = 18))

# wide_temp<- tempDF_full %>% 
#   pivot_wider( names_from   = month, values_from = max_temp) %>%
#   as.data.frame()

library(Hmisc)
rcorr(as.matrix(wide_temp))


# Wind direction and velocity #####################################################

# Get the wind components and dates
em_wind <- emmonak_env[,-c(1:2)]

em_wind <- emmonak_env[emmonak_env$year != 1996,]
# Replace M's with NA
em_wind$drct[em_wind$drct == "M"] <- NA
em_wind$sknt[em_wind$sknt == "M"] <- NA


em_wind$drct <- as.numeric(em_wind$drct)
em_wind$sknt <- as.numeric(em_wind$sknt)

# Calculate the u and v wind components
em_wind$u.wind <- -1*em_wind$sknt*sin(2*pi*em_wind$drct/360)
em_wind$v.wind <- -1*em_wind$sknt*cos(2*pi*em_wind$drct/360)

# Calculate 
years <- unique(em_wind$year)
 
# matrix to hold avg wind
wind.df <- matrix(ncol = 7,
                  nrow = length(years))

colnames(wind.df) <- c("year","month", "avg.wd", "vect.avg.ws", "scalar.avg.ws", "mean.u", "mean.v")

# month <- 5

# Loop through to get avg wind for each year over a specified month or time interval

for (m in c(3:6)) {
  

  for (y in 1:length(unique(em_wind$year))) {
  
    # m == 3
    # y == 1

  year <- years[y]
  
  # em_wind[em_wind$month >= 4 &
  #         em_wind$month <= 7 &
  #         em_wind$year == year,]



# Calculate avg wind vectors
  mean.u <- mean(em_wind$u.wind[em_wind$year == year &
                                  em_wind$month == m], na.rm = T)
  mean.v <- mean(em_wind$v.wind[em_wind$year == year &
                                  em_wind$month == m], na.rm = T)
  
  # Calculate the resultant vector avg wind direction
  wd.avg <- (atan2(mean.u, mean.v)*360/2/pi)+180
  
  # Calculate the vector avg wind speed
  ws.vect.avg <- ((mean.u^2 + mean.v^2)^0.5)
  
  # Calculate the scalar avg wind speed
  w.scalar.avg <- mean(em_wind$sknt[em_wind$year == year &
                                      em_wind$month == 5], na.rm = T)

  # 
  wind.df[y,1] <- year
  
  wind.df[y,2] <- m
  
  wind.df[y,3] <- wd.avg
  
  wind.df[y,4] <- ws.vect.avg
  
  wind.df[y,5] <- w.scalar.avg
  
  wind.df[y,6] <- mean.u
  
  wind.df[y,7] <- mean.v
  

  }
  if(m == 3){
    
    full_monthly_windDF <- wind.df
  }else{
    
    full_monthly_windDF <- rbind(full_monthly_windDF, wind.df)
  }
}

full_monthly_windDF <- as.data.frame(full_monthly_windDF)
midPoint <- all.timing.logistic[,c(2,5)]


full_windDF <- left_join(full_monthly_windDF, midPoint)

full_windDF$round.mid <- round(full_windDF$mid)

full_windDF$round.mid <- as.factor(full_windDF$round.mid)

for (i in 1:length(full_windDF[,1])) {
  

full_windDF$timing[i] <- if(full_windDF$mid[i] > mean(midPoint$mid)){"Late"}else{if(full_windDF$mid[i] < mean(midPoint$mid)){"Early"}else{"On Time"}}
  

}

full_windDF$timing <- as.factor(full_windDF$timing)

full_windDF$month.name <- month.name[full_windDF$month]

full_windDF$month.name <- factor(full_windDF$month.name, levels = c("March","April","May","June","July"))



ggplot(full_windDF, aes(x = avg.wd , y = mid, size = mid, col = timing))+
  geom_point()+
  scale_x_continuous(breaks=c(0,45,90,135,180,225,270,315),
                     labels=c("N","NE","E","SE","S","SW","W","NW"),
                     limits = c(0,360))+
  # stat_cor(method = "spearman",label.x = 180, label.y = 172)+
  coord_polar()+
  facet_wrap(~month.name)+
  labs(x = "Wind Direction",
       y = "Midpoint",
       col = "Timing",
       size = "Midpoint")+
theme(text = element_text(size = 18))


ggplot(full_windDF, aes(x = avg.wd, y = scalar.avg.ws , fill = timing))+
  geom_col(width = 10)+
  geom_point()+
  scale_x_continuous(breaks=c(0,45,90,135,180,225,270,315),
                     labels=c("N","NE","E","SE","S","SW","W","NW"),
                     limits = c(0,360))+
  # geom_point(aes(x = avg.wd, y = mid/20, color = mid))+
  coord_polar()+
  # stat_cor( )+
  facet_wrap(~month)+
  labs(x = "Wind Direction",
       y = "Scalar Avg Wind Speed knts",
       fill = "Timing")+
  theme(text = element_text(size = 18))



test <- full_windDF %>% drop_na() %>% subset(month == 5) %>% as.data.frame()

cor(test$mid, test$avg.wd)

ggplot(test, aes(x = avg.wd, y = mid))+
  geom_point()+
  geom_smooth()

cor(test$mid, test$vect.avg.ws)

cor(test$mid, test$scalar.avg.ws)

cor(test$mid, test$mean.u)

ggplot(test, aes(x = mean.u, y = mid))+
  geom_point()+
  geom_smooth(method = "lm")+  
  stat_cor(aes(label = after_stat(rr.label)), label.y = 181)

cor(test$mid, test$mean.v)

ggplot(test, aes(x = mean.v, y = mid))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(aes(label = after_stat(rr.label)))


# Polynomial model fitting avg.wd and midpoint
fit <- lm(test$mid~test$avg.wd + I(test$avg.wd^2))

summary(fit)

x <- seq(0,360,by = 1)
f.plot <- coef(fit)[1] + coef(fit)[2]*x + coef(fit)[3]*x^2

plot(f.plot~x,
     ylim = c(166,183), 
     type = "l",
     xlab = "Avg Wind Direction",
     ylab = "Mid Point",
     col = "red",
     lty = 2)
points(x = test$avg.wd, y = test$mid)

# Wind speed by mid point
fit.ws <- lm(test$mid~test$scalar.avg.ws)

summary(fit.ws)

xx <- seq(6.5,10,by = 0.25)
f.plot.ws <- coef(fit.ws)[1] + coef(fit.ws)[2]*xx 

plot(f.plot.ws~xx,
     ylim = c(165,185),
     # xlim = c(0,9),
     type = "l",
     # xlab = "Avg Wind Direction",
     # ylab = "Mid Point",
     col = "red",
     lty = 2)
points(x = test$scalar.avg.ws, y = test$mid)




library(mgcv)

# Simulate  the curve for fitting
sim.parabola <- function(wind, a, h, k){
  
  pred <-  a*(wind-h)^2+k
  
  return(pred)
  
}

# Negative log-likelihood function for mle2
NLL_logistic <- function(a, h, k, obs.mid, sigma){

  # sigma <- exp(ln_sigma)
  
  pred <- sim.parabola(wind = wind,
                         a = a,
                         h = h,
                         k = k)
  
  ll <- dnorm(x = obs.mid, mean = pred, sd = sigma, log = F)
  
  nll <- -1*sum(ll)
  
  return(nll)
}


# pred <- sim.parabola(wind = x, a = .0005, h = 200, k = 165)
# 
# 
# plot(x = test$avg.wd, y = test$mid,
#      ylim = c(165,182))
# lines(x = x, y = pred)


# Use mle2 to fit a line 
fit.parabola <- mle2(NLL_logistic,
                        start = list(
                          a = 0.0005,
                          k = 165,
                          h = 200,
                          sigma = 5
                        ),
                        data = list(wind = test$avg.wd,
                                    obs.mid = test$mid),
                        optimizer = "nlminb",
                        method = "Nelder-Mead",
                        control = list(maxit = 1e6))

# Get line prediction for plotting the results
pred.mid <- sim.parabola(wind = seq(from = 0, to = 360, by = 1), 
                          a = (coef(fit.parabola)[1]),
                          h = (coef(fit.parabola)[2]),
                         k = (coef(fit.parabola)[3]))

plot(x = test$avg.wd, y = test$mid,
     ylim = c(160,200))
lines(x = x, y = pred.mid)



SCRATCH #############################################################

wind.df <- as.data.frame(wind.df)

midPoint <- all.timing.logistic[,c(2,5)]

wind.cor.df <- left_join(wind.df,midPoint)

windy.df <- na.omit(wind.cor.df)

cor(windy.df$mid, windy.df$avg.wd)

# ggplot(windy.df, aes(x = avg.wd , y = mid, size = mid))+
#   geom_point()+
#   scale_x_continuous(breaks=c(0,45,90,135,180,225,270,315),
#                      labels=c("N","NE","E","SE","S","SW","W","NW"),
#                      limits = c(0,360))+
#   coord_polar()

# ggplot(windy.df, aes(x = avg.wd, y = scalar.avg.ws , fill = mid))+
#   geom_col(width = 10,)+
#   scale_x_continuous(breaks=c(0,45,90,135,180,225,270,315),
#                      labels=c("N","NE","E","SE","S","SW","W","NW"),
#                      limits = c(0,360))+
#   geom_point(aes(x = avg.wd, y = mid/20, color = mid))+
#   coord_polar()
  






april.wind <- em_wind[em_wind$month == 6 ,]

aprilwindAVG<-april.wind %>% group_by(year) %>% 
  summarise(mean.u = mean(u.wind, na.rm = T),
            mean.v = mean(v.wind, na.rm = T),
            wd.avg = (atan2(mean.u, mean.v)*360/2/pi)+180,
            month = month)


tail(aprilwindAVG)

midDF <- all.timing.logistic[,c(2,5)]

full_df <- left_join(midDF,aprilwindAVG)

ggplot(full_df, aes(x = wd.avg, y = mid))+
  geom_point()

cor(full_df$mid[full_df$year!=2020& full_df$year !=1995],
    full_df$wd.avg[full_df$year!=2020 &full_df$year!=1995])

ggplot(april.wind, aes(x = jday, y = sknt))+
  geom_line()+
  geom_point()



for (y in 1:12) {
  
  # y <- 2

  windDF <- em_wind[em_wind$month == y ,]

  monthly_windAVG <-  windDF %>% group_by(year) %>% 
    summarise(mean.u = mean(u.wind, na.rm = T),
            mean.v = mean(v.wind, na.rm = T),
            wd.avg = (atan2(mean.u, mean.v)*360/2/pi)+180,
            month = month)


# tail(aprilwindAVG)

  if(y > 1){
    
  full_monthly_df <- rbind(full_monthly_df, monthly_windAVG)
  
  }else{full_monthly_df <- monthly_windAVG}
  # ggplot(full_df, aes(x = wd.avg, y = mid))+
  #   geom_point()

  # cor(full_df$mid[full_df$year!=2020& full_df$year !=1995], full_df$wd.avg[full_df$year!=2020 &full_df$year!=1995])

   } 

full_monthly_df
timeAverage(em_wind, avg.time = "month")
?timeAverage

head(em_wind)

data.em.wind <- em_wind[c(1,3,4)]

data.em.wind$valid <- ymd_hm(data.em.wind$valid)

str(data.em.wind$valid)

clean.wind <- data.em.wind %>% rename(date = valid,
                                      wd = drct,
                                      ws = sknt)


clean.wind  <- clean.wind %>% mutate(year = lubridate::year(date),
                              month = lubridate::month(date),
                              jday = lubridate::yday(date))


april.wind <- timeAverage(clean.wind[clean.wind$year == 2017 &
                                       clean.wind$month == 4,], avg.time = "day" )

windRose(clean.wind,avg.time = , paddle = F, breaks = seq(from = 0,
                                                                 to = 60,
                                                                 by = 5))

ggplot(clean.wind[clean.wind$month == 4 & clean.wind$year == 2019,], aes(y = ws, x = wd, fill = ws))+
  geom_col()+
  coord_polar()

  
ggplot(april.wind, aes(y = ws, x = as.factor(wd), fill = year))+
  geom_col()+
  coord_polar()
