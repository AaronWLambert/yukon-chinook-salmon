#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 06.22.22
# Version: 5.0
# Purpose: To generate predictions from an integrated Bayesian inseason assessment
#           model.
#
#   1) Read in data
#   2) Preprocess the data for model
#   3) Fit a Logistic curve
#   3) Call to Stan model to generate inseason projection
#

#

#=================================================================================
#NOTES:

# This script is the general working script that is used to run single iterations 
#  of the model with any current Stan file.
# library(bbmle)
# library(mvtnorm)
library(rstan)
library(bayesplot)
library(tidyverse)
library(ggthemes)
library(viridis)
library(tidybayes)
library(tidyverse)
library(bbmle)

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
CAN_hist <- readRDS(file.path(dir.data,"Can EOS Abund 3Mar23.RDS"))

# PSS historical for Jun 1st = 152, to Aug 31 = 243
# Data obtained from ADFG website
# https://www.adfg.alaska.gov/index.cfm?adfg=commercialbyareayukon.salmon_escapement
# NOTE!!!! Use "Dataframe preprocess.R" to get  most recent version.
#  If not using Dataframe preprocess.R, uncomment out the following to load relevant PSS historical data
PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage Final 2022.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
# GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - 2021)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
pf_hist <- readRDS(file.path(dir.data,"pf_ver3.1_14Aprl23_eagle+harvest.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs) and days not covered by GSI
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 27April22.RDS"))

# Eagle passage
Eagle_hist <- readRDS(file = file.path(dir.data, "Eagle passage 23Nov22.RDS"))

# Logistic parameter fits from mle2
logistic.all <- read.csv(file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))

# Control Section ######################
# This is where the user can input dates, stan model, and stan model controls.
model.version <- "5.1"

# Range of years 1995 to 2021
myYear <- 2010

# Range of days 152 -243
myDay <- 158


# MCMC Parameters
n.chains <- 4
n.iter <- 2000;#5e4
n.thin <- 2

# Start Years for Predictors
startYearPF <- 2007
startYearGSI <- 2005
startYearPSS <-1997
startYearEagle <- 2005
# Wont typicaly change;
# day 152 = June 1
startDayPSS <- 148
endDayPSS <- 250
#####################################################################################################33

# # Vector containing avg GSI proportions for each day ACROSS ALL Years
# #  (for versions 3.0,3.1,3.2,3.3) ####
# meanGSI_vect <- vector(length = length(148:252))
# names(meanGSI_vect) <- c(148:252)
# sdGSI_vect <- vector(length = length(148:252))
# names(sdGSI_vect) <- c(148:252)
# 
# counter <- 1
# for (d in 148:252) {
#   # d = 175
#   meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
#   sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
#   counter <- counter+1
# }
# # 
# # # No GSI data for day 148 and 149. Use day 150 for these values
# meanGSI_vect[1:2] <- meanGSI_vect[3]
# sdGSI_vect[1:2] <- sdGSI_vect[3]
# # 
# # # No GSI data for days 249:252
# # #  add the last day available GSI information to these days
# meanGSI_vect[is.na(meanGSI_vect)] <- meanGSI_vect["248"]
# sdGSI_vect[is.na(sdGSI_vect)] <- sdGSI_vect["248"]
# 
# # 
# # 
# # # Put meanGSI_vect into dataframe to join with PSS_hist and adjust
# # # counts by GSI proportions for fitting curve in function
# gsiDF <- data.frame("GSI" = meanGSI_vect, "Day" = c(148:252))
# 
# gsiDF$GSI[is.na(gsiDF$GSI)] <- gsiDF$GSI[gsiDF$Day == 248]
# 
# PSS_hist_adj <- left_join(x = PSS_hist, y = gsiDF)
# 
# PSS_hist_adj$adj_count <- PSS_hist_adj$count*PSS_hist_adj$GSI
# 
# temp.df <- PSS_hist_adj %>% group_by(Year) %>% summarize(CumCount = cumsum(count))
# 
# PSS_hist_adj <- cbind(PSS_hist_adj,temp.df)
# 
# # Functions for fitting curve to daily passage ############################################
# extract.data <- function(df = PSS_hist_adj, year) {
#   ### TESTING ###
#   #baywide <- TRUE
#   #dist <- 324
#   # Year <- 2000
#   ###############
#   temp.data <- df[df$Year == year,]
#   min.day <- min(temp.data$Day)
#   max.day <- max(temp.data$Day)
#   days <- min.day:max.day
#   n.days <- length(days)
# 
#   #Calculate standardized date values
#   global.days <- min(df$Day):max(df$Day)
#   std.days <- which(global.days %in% days)
# 
#   rets <- vector(length=n.days)
# 
#   rets <- df$CumCount[df$Year==year]
#   #OUTPUT SECTION
#   output <- NULL
#   # if(baywide==TRUE) { output$dist <- NA }else { output$dist <- dist }
#   output$min.day <- min.day
#   output$max.day <- max.day
#   output$days <- days
#   output$n.days <- n.days
#   output$rets <- rets
#   output$std.days <- std.days
#   return(output)
# }
# 
# # out <- extract.data(PSS_hist = PSS_hist, year = 1995)
# # FUNCTION likelihood model for normal distribution
# #dist options: pois, norm, flynn, ssq
# like.norm.logistic <- function(days, rets, s, m, alpha, sigma, plot, dist='norm') {
#   ### TESTING ###
#   # rets <- out$rets
#   # days <- out$days
#   # sigma <- -2.3
#   # mu <- 5.25
#   # sd <- 12.26
#   # alpha <- 20.2
#   # plot <- TRUE
#   # dist <- 'norm'
#   ###############
#   if(dist!='norm' & dist!='pois' & dist!='flynn' & dist!='ssq')
#   { stop('IMPROPER ERROR DISTRIBUTION SELECTED') }
#   m <- exp(m)
#   s <- exp(s)
#   alpha <- exp(alpha)
#   sigma <- exp(sigma)
# 
#   #days <- out$days
#   #rets <- out$rets
#   n.days <- length(days)
# 
#   pred <- vector(length=n.days)
#   logLike <- vector(length=n.days)
# 
#   # d <- 5
#   # day <- 152
#   for(d in 1:n.days) {
#     day <- days[d]
#     #Predicted data
#     pred[d] <- alpha*(1/(1+exp((-(day-m))/s)))
# 
#     if(is.na(rets[d])) {  #NO DATA
#       logLike[d] <- 0
#     }else {
#       if(dist=='pois') {
#         logLike[d] <- dpois(rets[d],pred[d],log=TRUE)
#       }
#       if(dist=='norm'){
#         logLike[d] <- dnorm(rets[d], pred[d], sigma,log=TRUE)
#       }
#       if(dist=='flynn') {
#         flynn.sigma <- sigma*rets[d] + 0.1
#         logLike[d] <- dnorm(rets[d],pred[d], flynn.sigma, log=TRUE)
#       }
#       if(dist=='ssq') {
#         logLike[d] <- -1*(rets[d] - pred[d])^2
#       }
#     }
#   }#next d
#   #Plotting
#   if(plot==TRUE) {
#     y.lim <- c(0,max(rets,pred, na.rm=TRUE))
#     # plot(rets ~ days, type='l', col='gray', xlab='Day', ylab='Catch + Esc', lty=3, ylim=y.lim)
#     plot(rets ~ days, 
#          type='l', 
#          col='black',
#          xlab='Day', 
#          ylab='Number of Chinook Salmon', 
#          lty=1, 
#          ylim=y.lim, 
#          xlim = c(148,220),
#          lwd=5,
#          cex.axis = 1.5,
#          cex.lab = 1.5,
#          xaxt='n')
#     points(x=days, y=rets, pch=21, bg='gray')
#     #Model data
#     lines(x=days, y=pred, lwd=2, col='red')
#     axis(side=1, at=seq(from=152, to=202, by=10),
#          labels=c("June-1","June-11","June-21","July-1","July-11","July-21"),
#          cex = 1.5)
#     mtext(paste(year), side=3, outer=TRUE, line=-4, font=2)
#   }
# 
#   #Negative Log Likelihood
#   NLL <- -1*sum(logLike)
#   return(NLL)
# }
# # End of functions
# # length(c(1995,1997:2021))
# #NORMAL PRIOR ##############################################################
# norm.prior <- TRUE
# if(norm.prior==TRUE) {
#   years <- unique(PSS_hist_adj$Year[PSS_hist_adj$Year != myYear &
#                                       PSS_hist_adj$Year <= 2021])
#   n.years <- length(years)
# 
#   norm.m <- vector(length=n.years)
#   norm.s <- vector(length=n.years)
#   norm.alpha <- vector(length=n.years)
#   year.df <- vector(length = n.years)
#   # pdf(file = file.path(dir.output,"Meeting Logistic Curve fit 155 29Sept22.pdf"), height=8, width=8)
# 
#   # y <- 24
#   for(y in 1:n.years) {
#     year <- years[y]
#     print(year)
#     # Extract data
#     out <- extract.data(df = PSS_hist_adj, year=year)
#     #Fit normal data
#     fit.norm <-  mle2(like.norm.logistic,
#                       start=list(m=log(170), s=log(5.5), alpha=log(60000)),
#                       fixed=list(plot=FALSE, sigma=log(20)),
#                       data=list(days=out$days, rets=out$rets, dist='ssq'),
#                       method='Nelder-Mead', control=list(maxit=100000, trace=FALSE))
#     norm.s[y] <- exp(coef(fit.norm)[1])
#     norm.m[y] <- exp(coef(fit.norm)[2])
#     norm.alpha[y] <- exp(coef(fit.norm)[3])
#     year.df[y] <- year
#     like.norm.logistic(days=out$days,
#               rets=out$rets,
#               (coef(fit.norm)[1]),
#               coef(fit.norm)[2],
#               coef(fit.norm)[3],
#               coef(fit.norm)[4],
#               plot=TRUE,
#               dist='ssq')
#     # mtext(year, side=3, line=2, font=2)
#   }#next y
# 
#   # dev.off()
#   prior.df.log <- data.frame(norm.m,norm.s,norm.alpha, year.df)
#   names(prior.df.log) <- c('m','s','alpha',"year")
#   # write.csv(prior.df,file = file.path(dir.output,"logisitic curve prior.csv"))
# }
# prior.df<- read.csv(file = file.path(dir.output,"logisitic curve prior.csv"))
# Store each parameter as a vector for model input

ps_m <- logistic.all$mid[logistic.all$year != myYear]
ps_s <- logistic.all$sd[logistic.all$year != myYear]
ps_alpha.log <- logistic.all$alpha[logistic.all$year != myYear]

n_ps_m <- length(ps_m)
n_ps_s <- length(ps_s)
n_ps_alpha_log <- length(ps_alpha.log)

# Preseason Forecast######################

# Current Preseason forecast Can Origin
pf <- log(pf_hist$mean[pf_hist$Year == myYear])

# Vector of historical preseason forecasts for to compute sd for prior in stan model
PF_vect <- pf_hist$mean[pf_hist$Year <= 2021 &
                          pf_hist$Year != myYear]

names(PF_vect) <-  pf_hist$Year[pf_hist$Year <= 2021 &
                                  pf_hist$Year != myYear]


# EOS reconstructed runsize for historic years that have a preseason forecast 
EOS <- CAN_hist$can.mean[CAN_hist$Year >= 2007 & 
                           # CAN_hist$Year != 2011 & 
                           # CAN_hist$Year != 2012
                           CAN_hist$Year != myYear]

names(EOS) <- CAN_hist$Year[CAN_hist$Year >= 2007 & 
                              # CAN_hist$Year != 2011 & 
                              # CAN_hist$Year != 2012 &
                              CAN_hist$Year != myYear]

# SD for preseason forecast prior
pf_sigma <- sd(log(EOS/PF_vect))

# Years in PF
yearPF <- unique(pf_hist$Year[pf_hist$Year != myYear &
                                pf_hist$Year <= 2021])

n_yearsPF <- length(yearPF)

# Historic EOS Reconstructed Can-origin Run Size ############################################

# Vector of run sizes excluding the year of interest
totalEOS <- CAN_hist$can.mean[CAN_hist$Year <= 2021 
                              & CAN_hist$Year != 1996 &
                                CAN_hist$Year != myYear]

# Name the elements for accounting purposes
names(totalEOS) <- CAN_hist$Year[CAN_hist$Year <= 2021
                                 & CAN_hist$Year != 1996 &
                                   CAN_hist$Year != myYear]

# Number of years included in EOS
n_totalEOS <- length(totalEOS)

# Pilot Station Sonar Data ##############################################################

# Vector of historic PSS years excluding myYear
yearPSS <- unique(PSS_hist$Year[PSS_hist$Year != myYear &
                                  PSS_hist$Year <= max(CAN_hist$Year)])

# Number of years used in model
n_yearPSS <- length(yearPSS)

# # PSS days included up to myDay
dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay) & 
                                PSS_hist$Day>= startDayPSS])

# Number of days used
n_dayPSS <- length(dayPSS)

# Number of total days in the season
pss_days_all <- startDayPSS:endDayPSS

n_pss_days_all <- length(pss_days_all)

loc.allDays.myDay <- which(pss_days_all %in% myDay)

# PSS daily passage estimate for days up to myDay for the year of interest (myYear)
curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear &
                              PSS_hist$Day <= myDay &
                              PSS_hist$Day>=startDayPSS])

names(curr_PSS) <- (PSS_hist$Day[PSS_hist$Year == myYear &
                                   PSS_hist$Day <= myDay &
                                   PSS_hist$Day>=startDayPSS])

# Number of days used
n_curr_PSS <- length(curr_PSS)

# Cumulative current PSS counts
cum_curr_PSS <- cumsum(curr_PSS)


# Create matrix with dimensions days x historic Year
PSS_mat <- matrix(nrow = n_dayPSS,
                  #start 1996 to account for missing year 1996
                  ncol = n_yearPSS) 

# Give names to matrix for accounting
colnames(PSS_mat) <- yearPSS

rownames(PSS_mat) <- dayPSS

# Create vector of counts from relevant years and days to input into matrix

(count_vect <- PSS_hist$count[PSS_hist$Year != myYear &
                                PSS_hist$Year <= max(CAN_hist$Year) &
                                PSS_hist$Day <= (myDay) &
                                PSS_hist$Day >= (startDayPSS)])

# The number of observations in count_vector for PSS historic counts
# n_hist_counts <- length(count_vect)

# Use loop to populate matrix with counts by days
# Set counter
counter <-1

for (y in 1:n_yearPSS){
  for (d in 1:n_dayPSS) {
    
    
    PSS_mat[d,y] <- count_vect[counter]
    counter = counter+1
    
  }
}

# Cumulative PSS matrix
cum_PSS_mat <-  apply(PSS_mat, 2, cumsum)





cumPSS_all <- apply(PSS_mat_all, 2, sum)

# Eagle Sonar preprocessing #####################################################################

# Vector of historic Eagle Sonar passage excluding myYear
yearEagle <- unique(Eagle_hist$Year[Eagle_hist$Year != myYear &
                                      Eagle_hist$Year <= 2021])

# Number of Eagle sonar years used in model
n_yearEagle <- length(yearEagle)

# # Eagle days included up to myDay
#     151 is subtracted for use in logistic model
#      #  I.e., June 1 = 1, June 2 = 2 ....
# dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])-147 *******

dayEagle <- unique(Eagle_hist$Day[Eagle_hist$Day <=(myDay) & Eagle_hist$Day>= startDayPSS])

# Number of days used
n_dayEagle <- length(dayEagle)

# Eagle daily passage estimate for days up to myDay for the year of interest (myYear)
# curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay ]) **********
curr_Eagle_vect <- (Eagle_hist$count[Eagle_hist$Year == myYear &
                                       Eagle_hist$Day <= myDay & 
                                       Eagle_hist$Day>=startDayPSS])

curr_Eagle <- array(data = curr_Eagle_vect,
                    dim = length(curr_Eagle_vect))

# Number of days used
# Need this so stan will accept a 0 value for curr_pss
# if(length(curr_Eagle)==1 & curr_Eagle[1]<1){
#   curr_Eagle <- numeric(0)
# }

# Number of days used
n_curr_Eagle <- length(curr_Eagle)

# Cumulative Eagle counts for myYear
cum_curr_Eagle <- cumsum(curr_Eagle)

# Location pointers to index years in common for variance calculations

# 
loc.eagle.years <- which(yearPSS %in% yearEagle)
loc.eagle.years

loc.PF.years.eagle <- which(yearEagle %in% yearPF)

loc.PF.years.PSS <- which(yearPSS %in% yearPF)

# Make sure it works....
# length(loc.eagle.years)
# 
# yearPSS[loc.eagle.years]


# Create matrix with dimensions myday and myYear
Eagle_mat <- matrix(nrow = n_dayPSS,
                    #start 1996 to account for missing year 1996
                    ncol = n_yearEagle) 

# Give names to matrix for accounting
colnames(Eagle_mat) <- yearEagle

rownames(Eagle_mat) <- dayPSS

# Create vector of counts from relevant years and days to input into matrix
# (count_vect <- PSS_hist$count[PSS_hist$Year<= (myYear-1) & *******************************************
#                                 PSS_hist$Day <= (myDay)])
(count_vect_Eagle <- Eagle_hist$count[Eagle_hist$Year != myYear &
                                        Eagle_hist$Year <= 2021 &
                                        Eagle_hist$Day <= (myDay) &
                                        Eagle_hist$Day >= (startDayPSS)])

# The number of observations in count_vector for PSS historic counts
n_hist_Eaglecounts <- length(count_vect_Eagle)

# Use loop to populate matrix with counts by days
# Set counter

# This is to index where to start filling the Eagle matrix 
SE <- n_dayPSS- n_dayEagle + 1

#Counter for loop
counter <-1

# Loop to fill matix by year and day
for (y in 1:n_yearEagle){
  for (d in SE:n_dayPSS) {
    
    
    Eagle_mat[d,y] <- count_vect_Eagle[counter]
    
    counter = counter+1
    
  }
}

# Fill in the rest of the matrix with zeros
Eagle_mat[is.na(Eagle_mat)]<-0

# Cumulative count matrix for historic years
Eagle_cum_hist_mat <- apply(X = Eagle_mat,MARGIN = 2,FUN = cumsum)

# Cumulative Counts by year for plotting
cumEagle <- apply(Eagle_mat,MARGIN = 2,sum)



# Vector containing avg GSI proportions for each day ACROSS ALL Years ##############
#  (for versions 3.0,3.1,3.2,3.3) 

meanGSI_vect <- vector(length = length(startDayPSS:endDayPSS))

names(meanGSI_vect) <- c(startDayPSS:endDayPSS)

sdGSI_vect <- vector(length = length(startDayPSS:endDayPSS))

names(sdGSI_vect) <- c(startDayPSS:endDayPSS)

counter <- 1
for (d in startDayPSS:endDayPSS) {
  # d = 175
  meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  
  sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  
  counter <- counter+1
}

# No GSI data for day 148 and 149. Use day 150 for these values
meanGSI_vect[1:2] <- meanGSI_vect[3]

sdGSI_vect[1:2] <- sdGSI_vect[3]

# No GSI data for days 249:252
#  add the last day available GSI information to these days
meanGSI_vect[is.na(meanGSI_vect)] <- meanGSI_vect["248"]

sdGSI_vect[is.na(sdGSI_vect)] <- sdGSI_vect["248"]

# Mean GSI by mean strata dates. This is used in versions 3.4 and up #############################
meanStartDay <-GSI_by_year %>% group_by(stratum) %>% 
  summarise("meanStartDay" = mean(startday)) %>% 
  as.data.frame()

GSI_mean_by_strata <-GSI_by_year %>% 
  summarize("stratumMean" = c(mean(propCan[stratum == 1]),
                              mean(propCan[stratum == 2]),
                              mean(propCan[stratum == 3 | stratum == 4])),
            "stratumSD" = c(sd(propCan[stratum == 1]),
                            sd(propCan[stratum == 2]),
                            sd(propCan[stratum == 3 | stratum == 4]))) %>%  as.data.frame()

GSI <- cbind(GSI_mean_by_strata,round(meanStartDay[1:3,]))

GSI_avg<-c(rep(GSI$stratumMean[GSI$stratum == 1], 
               times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
           rep(GSI$stratumMean[GSI$stratum == 2], 
               times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
           rep(GSI$stratumMean[GSI$stratum == 3],
               times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))

GSI_avg_vect <- GSI_avg[1:(n_pss_days_all)]

GSI_sd <- c(rep(GSI$stratumSD[GSI$stratum == 1], 
                times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
            rep(GSI$stratumSD[GSI$stratum == 2], 
                times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
            rep(GSI$stratumSD[GSI$stratum == 3],
                times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))

GSI_sd_vect <- GSI_sd[1:(n_pss_days_all)]

# Create matrix with dimensions myday and myYear
PSS_mat_adj <- matrix(nrow = length(startDayPSS:myDay),
                      #start 1996 to account for missing year 1996
                      ncol = n_yearPSS) 

# Give names to matrix
colnames(PSS_mat_adj) <- yearPSS
rownames(PSS_mat_adj) <- c(startDayPSS:(myDay))


for (y in 1:n_yearPSS) {
  for (d in 1:n_dayPSS) {
    # y = 1
    # d = 15
    PSS_mat_adj[d,y] <- PSS_mat[d,y]*meanGSI_vect[d] 
    
  }
  
}

# # Put meanGSI_vect into dataframe to join with PSS_hist and adjust
# # counts by GSI proportions for fitting curve in function
# gsiDF <- data.frame("GSI" = meanGSI_vect, "Day" = c(148:252))
# 
# gsiDF$GSI[is.na(gsiDF$GSI)] <- gsiDF$GSI[gsiDF$Day == 248]
# 
# PSS_hist_adj <- left_join(x = PSS_hist, y = gsiDF)
# 
# PSS_hist_adj$adj_count <- PSS_hist_adj$count*PSS_hist_adj$GSI


adj_curr_PSS <- vector(length = n_dayPSS)

for (d in 1:n_dayPSS) {
  
  
  adj_curr_PSS[d] <- curr_PSS[d]*meanGSI_vect[d]
}


# Create matrix with dimensions myday and myYear
PSS_mat_all_adj <- matrix(nrow = (length(startDayPSS:endDayPSS )),
                          #start 1996 to account for missing year 1996
                          ncol = n_yearPSS) 

# Give names to matrix for accounting
colnames(PSS_mat_all_adj) <- yearPSS

rownames(PSS_mat_all_adj) <- c(startDayPSS:endDayPSS)

# Use loop to populate matrix with counts by days
# Set counter
counter <-1

for (y in 1:n_yearPSS){
  for (d in 1:(length(startDayPSS:endDayPSS))) {
    
    
    PSS_mat_all_adj[d,y] <- PSS_mat_all[d,y] * meanGSI_vect[d]
    
  }
}



# Inits list for ver 5.0
init_fn <-function(){

  list("ps_alpha_curr" = runif(1,1e5,11e4),
       "ps_mid_curr" = runif(1,168,170),
       "ps_shape_curr" = runif(1,2,3),
       "ps_shape_hist" = runif(n_yearPSS,2,3),
       "ps_mid_hist" = runif(n_yearPSS,168,170),
       "ps_alpha_hist"=runif(n_yearPSS,1e5,11e4),
       "sigma" = runif(1,2,3),
       "sigma_hist" = runif(n_yearPSS,2,3),
       "beta" = runif(1,0,0.5),
       "sigma_reg" = runif(1,0.2,0.5),
       "alpha" = runif(1,10,11)
       
)}

# Inits list for ver 6.0
# init_fn <-
# # 
#   list("ps_alpha_curr" = runif(1,50000,200000),
# "ps_mid_curr" = runif(1,175,180),
# "ps_mid_hist" = runif(n_yearPSS,175,180),
# "ps_alpha_hist"=runif(n_yearPSS,50000,200000),
# "alpha_s" = runif(1,14,20),
# "beta_s" = runif(1,-.1,0)
# )
# 
init_ll <- list(init_fn() , init_fn(), init_fn(), init_fn())


# Run stan model ##################################################
fit <- stan(file = file.path(dir.stan,paste("Yukon Inseason Forecast ",
                                            model.version ,".stan", sep = "")),
            data = list("PSS_mat"=PSS_mat,
                        "n_totalEOS"=n_totalEOS,
                        "totalEOS"=totalEOS, 
                        "n_dayPSS"=n_dayPSS,
                        "dayPSS"=dayPSS, 
                        "n_yearPSS"=n_yearPSS,
                        "yearPSS"=yearPSS,
                        "n_curr_PSS"=n_curr_PSS,
                        "curr_PSS"=curr_PSS,
                        # "n_hist_counts"=n_hist_counts,
                        "Pf" = pf, 
                        "Pf_sigma" = pf_sigma,
                        "n_yearsPF" = n_yearsPF,
                        "meanpropCAN" = meanGSI_vect,
                        "sd_meanpropCAN" = sdGSI_vect,
                        "mean_adj_PSS_mat" = PSS_mat_adj,
                        "mean_adj_curr_PSS"= adj_curr_PSS,
                        "cum_PSS_mat" = cum_PSS_mat,
                        "cum_curr_PSS" = cum_curr_PSS,
                        "GSI_mean"= GSI_avg_vect,
                        "GSI_sd" = GSI_sd_vect,
                        # "ps_mu"= ps_mu,
                        # "ps_sd"= ps_sd,
                        "ps_alpha_log"= ps_alpha.log,
                        "pss_days_all" = pss_days_all,
                        "n_pss_days_all" = n_pss_days_all,
                        "PSS_mat_adj" = PSS_mat_adj,
                        "myDay" = myDay,
                        "PSS_mat_all" = PSS_mat_all,
                        "ps_m"= ps_m,
                        "ps_s"= ps_s,
                        "n_ps_m" = n_ps_m,
                        "n_ps_s" = n_ps_s,
                        "n_ps_alpha_log" = n_ps_alpha_log,
                        "loc_pf_years_Eagle" = loc.PF.years.eagle,
                        "loc_pf_years_PSS" = loc.PF.years.PSS,
                        "loc_allDays_myDay" = loc.allDays.myDay ),
            init = init_ll,
            chains = n.chains,
            # chains = 1,
            iter = n.iter, 
            thin = n.thin, 
            # cores = n.chains,
            cores = mc.cores,
            control = list(max_treedepth = 15, adapt_delta = 0.8),
            # verbose = F
            
)


# saveRDS(object = fit, file = file.path(dir.output,"fit5.1_7Nov22_203_2019.RDS"))

# fit <- readRDS(file=file.path(dir.output,"fit6.0_19Sept22_185.RDS"))
# Save to PDF
# pdf(file = file.path(dir.output,"6.0_plots_29sept22_logistic_day1752019.pdf"))

# traceplot(object = fit, c("ps_mid_curr", "ps_shape_curr", "ps_alpha_curr","alpha",
#                           "beta","ln_RunSize", "ps_alpha_hist[4]", "ps_mid_hist[4]", "ps_shape_hist[4]"))
# 
# traceplot(object = fit, c("ps_mid_curr", "pred_shape_curr", "ps_alpha_curr","alpha_s",
#                            "beta_s","ln_RunSize", "ps_alpha_hist[4]", "ps_mid_hist[4]", "pred_shape_hist[4]"))
# traceplot(object = fit, c("ps_mid_curr", "ln_pred_shape_curr", "ps_alpha_curr","alpha",
#                           "beta","ln_RunSize"))
# fit
# shinystan::launch_shinystan(as.shinystan(fit))
# mcmc_pairs(fit, pars = c("ps_mu_curr", "ps_sd_curr", "ps_alpha_curr"))

# Extract parameter estimates
pars <- rstan::extract(fit)


# For checking retro results
# pars <- model.output$pars
# dayPSS <- model.output$dayPSS
# cum_curr_PSS <- model.output$cum_curr_PSS
# myDay <- 180
# myYear <- 2007
# pss_days_all <- model.output$pss_days_all
# n_yearPSS <- 26
# cumPSS <- model.output$cumPSS
# cumPSS_all <- model.output$cumPSS_all
# (totalEOS <- model.output$totalEOS)


# # Limits
x.lim <- c(145,225)
# y.lim <- c(0,max(apply(pars$ps_pred_curr, 2, quantile, probs=0.975), max(cumsum(PSS_hist$count[PSS_hist$Year == myYear]))))
y.lim <- c(0,325000)

# Plot Observed
plot(x=dayPSS,
     y=cum_curr_PSS,
     type='p', pch=21, 
     col='gray',
     xlim=x.lim,
     ylim=y.lim,
     xlab='',
     ylab='',
     # main="PSS: Current Year", 
     xaxt='n',
     cex.axis = 1.45)

axis(side=1, at=seq(from=151, to=206, by=5), labels = FALSE, col='darkgray')
axis(side=1, at=seq(from=161, to=201, by=10),
     labels=c("June-10","June-20","June-30","July-10","July-20"),
     cex.axis = 1.25)
segments(x0=dayPSS, y0=rep(0, length(dayPSS)), x1=dayPSS, y1=cum_curr_PSS, col='black', lwd=2)
points(x=dayPSS, y=cum_curr_PSS, type='p', pch=21, bg='gray')
points(x = c((myDay+1):252), y = cumsum(PSS_hist$count[PSS_hist$Year == myYear])[((myDay+2)-148):(252-147)],
       type = "p",pch = 23, bg = "green")
# Plot Predicted
polygon(x=c(pss_days_all,rev(pss_days_all)), y=c(apply(pars$ps_pred_curr, 2, quantile, probs=0.025),
                                                 rev(apply(pars$ps_pred_curr, 2, quantile, probs=0.975))), border=FALSE, col=rgb(1,0,0,0.2))
polygon(x=c(pss_days_all,rev(pss_days_all)), y=c(apply(pars$ps_pred_curr, 2, quantile, probs=0.25),
                                                 rev(apply(pars$ps_pred_curr, 2, quantile, probs=0.75))), border=FALSE, col=rgb(1,0,0,0.2))
lines(x=pss_days_all, y=apply(pars$ps_pred_curr, 2, median), col=rgb(1,0,0, alpha=0.5), lwd=3)

abline(v=mean(ps_m), col='blue', lty=2, lwd=2)
# abline(h = sum(PSS_hist$count[PSS_hist$Year==myYear]),lty = 4, col = "Green", lwd = 2)
abline(h = mean(ps_alpha.log),lty = 4, col = "Black", lwd = 2)
# abline(h=CAN_hist$can.mean[CAN_hist$Year==myYear])

# Labels
# mtext(paste(" Day,",myDay,'of',myYear), side=3, outer=TRUE, line=-4, font=2)
# mtext('Day of Year', side=1, outer=TRUE, line=0.5, font=2)
# mtext('Cummualtive PSS Chinook Salmon', side=2, outer=TRUE, line=-2, font=2, cex = 1.5)
# legend('topleft', legend=c('Avg. Alpha','Pred.', 'Avg. MidPoint',"Realized "),
#        bg=c(NA,NA,NA,"green"),
#        col=c('black','red','blue','green'), 
#        lty=c(1,1,2,NA),
#        pch = c(NA,NA,NA,19),
#        )

# Historical
# y <- 3
for(y in 1:n_yearPSS) {
  # y <- 1
  
  year <- yearPSS[y]
  
  # Limits
  x.lim <- c(145,255)
  # y.lim <- c(0,max(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.975), cum_PSS_mat[,y]))
  y.lim <- c(0,400000)
  # Plot Observed
  plot(x=dayPSS, y=cum_PSS_mat[,y], type='p', pch=21, col='gray',
       xlim=x.lim,
       ylim=y.lim,
       xlab='', ylab='',
       main=year, xaxt='n')
  
  axis(side=1, at=seq(from=151, to=206, by=5), labels = FALSE, col='darkgray')
  
  axis(side=1, at=seq(from=161, to=201, by=10),
       labels=c("June-10","June-20","June-30","July-10","July-20"))
  
  segments(x0=dayPSS, y0=rep(0, length(dayPSS)), x1=dayPSS, y1=cum_PSS_mat[,y], col='black', lwd=2)
  
  points(x=dayPSS, y=cum_PSS_mat[,y], type='p', pch=21, bg='gray')
  
  # points(x = c((myDay+1):252), y = cumsum(PSS_hist$count[PSS_hist$Year == year])[((myDay+2)-148):(252-147)],
  #        type = "p",pch = 23, bg = "green")
  # 
  # points(x = c((myDay+1):252), y = cumsum(PSS_hist$count[PSS_hist$Year == year])[((myDay+2)-148):(270-myDay)],
  #        type = "p",pch = 23, bg = "green")
  points(x = c((myDay+1):252), y = cumsum(PSS_hist$count[PSS_hist$Year == year])[((myDay+2)-148):(252-147)],
         type = "p",pch = 23, bg = "green")
  # Plot Predicted
  polygon(x=c(pss_days_all,rev(pss_days_all)), 
          y=c(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.025),
                                                   rev(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.975))), 
          border=FALSE, col=rgb(1,0,0,0.2))
  polygon(x=c(pss_days_all,rev(pss_days_all)), y=c(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.25),
                                                   rev(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.75))), border=FALSE, col=rgb(1,0,0,0.2))
  lines(x=pss_days_all, y=apply(pars$ps_pred_hist[,,y], 2, median), col=rgb(1,0,0, alpha=0.5), lwd=3)
  # abline(v=mean(ps_m), col='blue', lty=2, lwd=2)
  abline(h = sum(PSS_hist$count[PSS_hist$Year==year]),lty = 4, col = "Green", lwd = 2)
  # Labels
  mtext(paste("Day",myDay,'of',myYear), side=3, outer=TRUE, line=-5, font=2)
  mtext('Day of Year', side=1, outer=TRUE, line=-2.5, font=2)
  # mtext('Can-origin PSS Daily Chinook Passage', side=2, outer=TRUE, line=-1.7, font=2)
  # legend('bottomright', legend=c('Obs.','Pred.', 'Avg. MidPoint',"Realized PSS Total Passage"),
  #        col=c('black','red','blue',"green"), lty=c(1,1,2,4))
}

# hist((ps_m))
# hist((pars$ps_mid_hist))
# hist((ps_s))
# hist(pars$ps_shape_hist)
# # hist(exp(pars$ln_pred_shape_hist))
# hist(pars$ps_alpha_hist)
# hist(ps_alpha.log)

# Extract estimates for prior, likelihood, and posterior 
df.run <- data.frame("par" = "Runsize", "value" = pars$RunSize)
df.prior <- data.frame("par" = "Prior", "value" = pars$prior_pf)

df.pssPred <- data.frame("par" = "PSS_Pred", "value" = pars$curr_predPSS)

df.postPredPss <- data.frame("par" = "PSS_post_pred", "value" = pars$post_curr_predPSS)
# Bind into one dataframe with 2 columns
# df.comb <- rbind(df.run,df.prior, df.pssPred, df.postPredPss)
df.comb <- rbind(df.prior,df.postPredPss,df.run)
levels(df.comb$par)
summary(df.comb)
# Density plots comparing pf, linear prediction, and posterior estimate
ggplot(df.comb, aes(x = value/1000, fill = par))+
  geom_density(alpha = .5)+
  ggtitle(paste("Density plot for",myYear,", day",myDay,"\n Version", model.version))+
  theme_classic()+
  geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
                 color = "End of Season Runsize"),
             linetype = 2,
             size = 1)+
  # xlim(c(0,300000))+
  coord_cartesian(xlim = c(0,250))+
  # ylim(c(0,500))
  labs(x = "1000's of Chinook Salmon", fill = "Parameters", y = "Probability Density")+
  scale_fill_discrete(name = "Parameters",
                      labels = c( "Preseason Forecast (Prior)",
                                  "PSS Prediction","Runsize"))+
  scale_x_continuous()+
  scale_fill_colorblind(name = "", 
                        labels = c( "Preseason Forecast (Prior)", 
                                    "PSS Prediction","Runsize"))+
  theme_tidybayes()+
  theme(text = element_text(size = 20, family = "serif"),
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")+
  geom_vline(aes(xintercept = median(df.run$value)/1000, color = "Median Runsize"),
             size = 1, linetype = 2)+
  geom_vline(aes(xintercept = median(df.prior$value)/1000, color = "Median PF"),
             size = 1, linetype = 2)+
  geom_vline(aes(xintercept = pf_hist$mean[pf_hist$Year==myYear]/1000, color = "PF Point"),
             size = 1, linetype = 2)

# Get the quantiles for predPSS from model output
quant.predPSS <- apply(X = pars$predPSS, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# quant.predPSS <- quantile(pars$post_curr_predPSS, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

cumPSS <- apply(PSS_mat_all,MARGIN = 2,sum)
# Prder for polygon
ord <- order(cumPSS_all)


adj_cumPSS <- vector(length = n_yearPSS)

for (i in 1:n_yearPSS) {
  adj_cumPSS[i]<-sum(PSS_mat_adj[,i])
}

# Call plot
plot(x = cumPSS_all, y = totalEOS,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = paste("Cumulative PSS Passage Day EOS \n Years 1995-",myYear-1 ),
     ylab = paste("Totoal Canadian-origin Chinook 1995-",myYear
                  -1),
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     xlim = c(0,max(cumPSS_all)),
     ylim = c(min(quant.predPSS), max(quant.predPSS)))
polygon(x = c(cumPSS_all[ord], rev(cumPSS_all[ord])),
        y = c(quant.predPSS[1,ord],rev(quant.predPSS[5,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS_all[ord], rev(cumPSS_all[ord])),
        y = c(quant.predPSS[2,ord],rev(quant.predPSS[4,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS_all[ord], quant.predPSS[3,ord], col = "red", lw = 2)


#
quantile(pars$RunSize, probs = c(.025,.1,.5,.9,.975))

# dev.off()

# Regression plot of shape parameter and mp parameter
# Get the quantiles for predPSS from model output
quant.predPSS_shape <- apply(X = (pars$ln_pred_shape_like), 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# quant.predPSS <- quantile(pars$post_curr_predPSS, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pred.plot.df <- data.frame(ps_m,
                           ps_s,
                           t(quant.predPSS_shape))
names(pred.plot.df) <- c("MidPoint", "Shape","low95","low50","median",
                         "up50","up95")



ggplot(pred.plot.df, aes(x = MidPoint, y = Shape))+
  geom_point() +
  geom_ribbon(aes(ymin=low95, ymax=up95,fill = "HDI 95"), alpha=0.25) +
  geom_ribbon(aes(ymin=low50, ymax=up50, fill = "HDI 50"), alpha=0.25) +
  geom_line(aes(y=median, color = "Median"), show.legend = T) +
  # labs(x = "Cummulative PSS Chinook Passage (1000's)",
  #      y = "Total Reconstructed EOS Chinook Passage (1000's)")+
  theme(text = element_text(size = 20, family = "serif"),
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_fill_colorblind(name = "")+
  scale_color_colorblind(name = "")
# cumPSS <- apply(PSS_mat_all,MARGIN = 2,sum)
# Prder for polygon
ord_m <- order(ps_m)


# adj_cumPSS <- vector(length = n_yearPSS)

# for (i in 1:n_yearPSS) {
#   adj_cumPSS[i]<-sum(PSS_mat_adj[,i])
# }

# Call plot
plot(x = ps_m, y = (ps_s),
     type = "p",
     pch = 21,
     bg = "red")
     # xlab = paste("Cumulative PSS Passage Day ", myDay, "\n Years 1995-",myYear-1 ),
     # ylab = paste("Totoal Canadian-origin Chinook 1995-",myYear
                  # -1),
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     # xlim = c(0,max(ps_s)),
     # ylim = c(min(quant.predPSS_shape), max(quant.predPSS_shape)))
polygon(x = c(ps_m[ord], rev(ps_m[ord])),
        y = c(quant.predPSS_shape[1,ord],rev(quant.predPSS_shape[5,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(ps_m[ord], rev(ps_m[ord])),
        y = c(quant.predPSS_shape[2,ord],rev(quant.predPSS_shape[4,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = ps_m[ord], quant.predPSS_shape[3,ord], col = "red", lw = 2)
boxplot(pars$propCAN)
fit
summary(fit)
mean(pars$sigma_pred)
