#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 01.22.23
# Purpose: Explore timing differences between Can-orig Chinook and total Chinook


library(bbmle)
# library(mvtnorm)
library(tidyverse)
library(rstan)
library(bayesplot)
library(ggthemes)
library(viridis)
library(tidybayes)
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

# Control Section ######################

# Range of years 1995 to 2021
# myYear <- 2017

# # Range of days 152 -243
# myDay <- 180
# 
# # MCMC Parameters
# n.chains <- 4
# n.iter <- 2000;#5e4
# n.thin <- 2

# Start Years for Predictors
startYearPF <- 2007
startYearGSI <- 2005
startYearPSS <-1997
startYearEagle <- 2005

# Wont typicaly change;
# day 152 = June 1
startDayPSS <- 148
# endDayPSS <- 252
endDayPSS <- 225
#####################################################################################################33

# Vector containing avg GSI proportions for each day ACROSS ALL Years
#  (for versions 3.0,3.1,3.2,3.3) ####
meanGSI_vect <- vector(length = length(startDayPSS:endDayPSS))
names(meanGSI_vect) <- c(startDayPSS:endDayPSS)
sdGSI_vect <- vector(length = length(startDayPSS:endDayPSS))
names(sdGSI_vect) <- c(startDayPSS:endDayPSS)

counter <- 1
for (d in startDayPSS:endDayPSS) {
  # d = 175
  meanGSI_vect[counter] <- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  sdGSI_vect[counter] <- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  counter <- counter+1
}
# 
# # No GSI data for day 148 and 149. Use day 150 for these values
meanGSI_vect[1:2] <- meanGSI_vect[3]
sdGSI_vect[1:2] <- sdGSI_vect[3]
# 
# # No GSI data for days 249:252
# #  add the last day available GSI information to these days
meanGSI_vect[is.na(meanGSI_vect)] <- meanGSI_vect["248"]
sdGSI_vect[is.na(sdGSI_vect)] <- sdGSI_vect["248"]

# # # Create matrix with dimensions myday and myYear
# # PSS_mat_adj <- matrix(nrow = (length(148:252 )),
# #                       #start 1996 to account for missing year 1996
# #                       ncol = length(1996:(2021))) 
# # 
# # # Give names to matrix
# # colnames(PSS_mat_adj) <- c("1995", 1997:(2021))
# # rownames(PSS_mat_adj) <- c(148:(252))
# # 
# # 
# # for (y in 1:length(c(1995,1997:2021))) {
# #   for (d in 1:length(148:252)) {
# #     # y = 1
# #     # d = 15
# #     PSS_mat_adj[d,y] <- PSS_mat[d,y]*meanGSI_vect[d] 
# #     
# #   }
# #   
# # }
# 
# 
# # Put meanGSI_vect into dataframe to join with PSS_hist and adjust
# # counts by GSI proportions for fitting curve in function
gsiDF <- data.frame("GSI" = meanGSI_vect, "Day" = c(startDayPSS:endDayPSS))

gsiDF$GSI[is.na(gsiDF$GSI)] <- gsiDF$GSI[gsiDF$Day == endDayPSS]

PSS_hist_adj <- left_join(x = PSS_hist, y = gsiDF)

PSS_hist_adj$adj_count <- PSS_hist_adj$count*PSS_hist_adj$GSI

temp.df <- PSS_hist_adj %>% group_by(Year) %>% summarize(CumCount = cumsum(count),
                                                         CumCount_adj = cumsum(adj_count))

PSS_hist_adj <- cbind(PSS_hist_adj,temp.df)

# # # Functions for fitting curve to daily passage ############################################
# extract.data <- function(df = PSS_hist, year, GSI = TRUE) {
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
#   
#   if(GSI==TRUE){rets <- df$adj_count[df$Year==year]}else{ rets <- df$count[df$Year==year]}
#   
#   
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
# #
# # out <- extract.data(PSS_hist = PSS_hist, year = 1995)
# # FUNCTION likelihood model for normal distribution
# #dist options: pois, norm, flynn, ssq
# like.norm <- function(days, rets, mu, sd, alpha, sigma, plot, dist='norm') {
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
#   mu <- exp(mu)
#   sd <- exp(sd)
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
#     pred[d] <- alpha*dnorm(day,mu,sd)
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
#     #plot(rets ~ days, type='l', col='gray', xlab='Day', ylab='Catch + Esc', lty=3, ylim=y.lim)
#     plot(rets ~ days, 
#          type='h', 
#          col='black', 
#          xlab='Day', 
#          ylab='PSS Chinook Salmon', 
#          lty=1, 
#          ylim=y.lim,
#          xlim = c(148,220),
#          lwd=5,
#          cex.axis = 1.5,
#          cex.lab = 1.5,
#          xaxt='n')
#     #points(x=days, y=rets, pch=21, bg='gray')
#     #Model data
#     lines(x=days, y=pred, lwd=2, col='red')
#     # axis(side=1, at=seq(from=151, to=206, by=5), labels = FALSE, col='darkgray')
#     axis(side=1, at=seq(from=152, to=202, by=10),
#          labels=c("June-1","June-11","June-21","July-1","July-11","July-21"),
#          cex = 1.5)
#   }
#   
#   #Negative Log Likelihood
#   NLL <- -1*sum(logLike)
#   return(NLL)
# }
# # # End of functions
# # length(c(1995,1997:2021))
# 
# # Normal Curve timing Canadian Chinook ####################################3
# norm.prior <- TRUE
# if(norm.prior==TRUE) {
#   # years <- c(1995,1997:(myYear-1))
#   years <- unique(PSS_hist_adj$Year)
#   n.years <- length(years)
#   
#   norm.mu <- vector(length=n.years)
#   norm.sd <- vector(length=n.years)
#   norm.alpha <- vector(length=n.years)
#   # pdf(file = file.path(dir.output,"4.0 Normal MLE Fit 15Jul22"), height=8, width=8)
#   # par(mfrow = c(3,3))
#   # y <- 1
#   for(y in 1:n.years) {
#     year <- years[y]
#     print(year)
#     # Extract data NOTE! Change GSI if needed for model
#     out <- extract.data(df = PSS_hist_adj, year=year, GSI = TRUE)
#     #Fit normal data
#     fit.norm <-  mle2(like.norm,
#                       # start=list(mu=log(170), sd=log(20), alpha=log(4000)),
#                       start=list(mu=log(170), sd=log(20), alpha=log(6000)),
#                       fixed=list(plot=FALSE, sigma=log(20)),
#                       data=list(days=out$days, rets=out$rets, dist='ssq'),
#                       method='Nelder-Mead', control=list(maxit=100000, trace=FALSE))
#     norm.mu[y] <- exp(coef(fit.norm)[1])
#     norm.sd[y] <- exp(coef(fit.norm)[2])
#     norm.alpha[y] <- exp(coef(fit.norm)[3])
#     like.norm(days=out$days,
#               rets=out$rets,
#               (coef(fit.norm)[1]),
#               coef(fit.norm)[2],
#               coef(fit.norm)[3],
#               coef(fit.norm)[4],
#               plot=TRUE,
#               dist='ssq')
#     mtext(year, side=3, line=2, font=2)
#   }#next y
#   
#   # dev.off()
#   prior.df <- data.frame(norm.mu,norm.sd,norm.alpha,years)
#   names(prior.df) <- c('mid','sd','alpha',"year")
#   
#   # write.csv(x = prior.df, file = file.path(dir.output,"normal curve prior GSI.csv"))
# }
# 
# # write.csv(x = prior.df, file = file.path(dir.output,"normal curve parameters CAN Chinook 1995_2022.csv"))
# 
# # Normal Curve timing Canadian Chinook ####################################
# norm.prior <- TRUE
# if(norm.prior==TRUE) {
#   # years <- c(1995,1997:(myYear-1))
#   years <- unique(PSS_hist_adj$Year)
#   n.years <- length(years)
#   
#   norm.mu <- vector(length=n.years)
#   norm.sd <- vector(length=n.years)
#   norm.alpha <- vector(length=n.years)
#   # pdf(file = file.path(dir.output,"4.0 Normal MLE Fit 15Jul22"), height=8, width=8)
#   # par(mfrow = c(3,3))
#   # y <- 1
#   for(y in 1:n.years) {
#     year <- years[y]
#     print(year)
#     # Extract data NOTE! Change GSI if needed for model
#     out <- extract.data(df = PSS_hist_adj, year=year, GSI = FALSE)
#     #Fit normal data
#     fit.norm <-  mle2(like.norm,
#                       # start=list(mu=log(170), sd=log(20), alpha=log(4000)),
#                       start=list(mu=log(170), sd=log(20), alpha=log(6000)),
#                       fixed=list(plot=FALSE, sigma=log(20)),
#                       data=list(days=out$days, rets=out$rets, dist='ssq'),
#                       method='Nelder-Mead', control=list(maxit=100000, trace=FALSE))
#     norm.mu[y] <- exp(coef(fit.norm)[1])
#     norm.sd[y] <- exp(coef(fit.norm)[2])
#     norm.alpha[y] <- exp(coef(fit.norm)[3])
#     like.norm(days=out$days,
#               rets=out$rets,
#               (coef(fit.norm)[1]),
#               coef(fit.norm)[2],
#               coef(fit.norm)[3],
#               coef(fit.norm)[4],
#               plot=TRUE,
#               dist='ssq')
#     mtext(year, side=3, line=2, font=2)
#   }#next y
#   
#   # dev.off()
#   prior.df <- data.frame(norm.mu,norm.sd,norm.alpha,years)
#   names(prior.df) <- c('mid','sd','alpha',"year")
#   
#   # write.csv(x = prior.df, file = file.path(dir.output,"normal curve prior GSI.csv"))
# }
# 
# # write.csv(x = prior.df, file = file.path(dir.output,"normal curve parameters All Chinook 1995_2022.csv"))
# 
# 
# 
# 
# # Logistic Curve #################################################################################
# # Functions for fitting curve to daily passage ############################################
# extract.data.logistic <- function(df = PSS_hist_adj, year, GSI = TRUE) {
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
#   if(GSI==TRUE){rets <- df$CumCount_adj[df$Year==year]}else{ rets <- df$CumCount[df$Year==year]}
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
# 
# # logistic curve parameters (Canadian Chinook) ###########################
# norm.prior <- TRUE
# if(norm.prior==TRUE) {
#   years <- unique(PSS_hist_adj$Year)
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
#     out <- extract.data.logistic(df = PSS_hist_adj, year=year, GSI = TRUE)
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
#                        rets=out$rets,
#                        (coef(fit.norm)[1]),
#                        coef(fit.norm)[2],
#                        coef(fit.norm)[3],
#                        coef(fit.norm)[4],
#                        plot=TRUE,
#                        dist='ssq')
#     # mtext(year, side=3, line=2, font=2)
#   }#next y
#   
#   # dev.off()
#   prior.df.log <- data.frame(norm.m,norm.s,norm.alpha, year.df)
#   names(prior.df.log) <- c('mid','sd','alpha',"year")
#   # write.csv(prior.df,file = file.path(dir.output,"logisitic curve prior.csv"))
# }
# 
# # write.csv(prior.df.log, file = file.path(dir.output,"logistic curve parameters CAN Chinook 1995_2022.csv"))
# 
# # logistic curve parameters (ALL Chinook) ###################################
# norm.prior <- TRUE
# if(norm.prior==TRUE) {
#   years <- unique(PSS_hist_adj$Year)
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
#     out <- extract.data.logistic(df = PSS_hist_adj, year=year, GSI = FALSE)
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
#                        rets=out$rets,
#                        (coef(fit.norm)[1]),
#                        coef(fit.norm)[2],
#                        coef(fit.norm)[3],
#                        coef(fit.norm)[4],
#                        plot=TRUE,
#                        dist='ssq')
#     # mtext(year, side=3, line=2, font=2)
#   }#next y
#   
#   # dev.off()
#   prior.df.log <- data.frame(norm.m,norm.s,norm.alpha, year.df)
#   names(prior.df.log) <- c('mid','sd','alpha',"year")
#   # write.csv(prior.df,file = file.path(dir.output,"logisitic curve prior.csv"))
# }

# write.csv(prior.df.log, file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))

# Import curve parameters ###########################################################################

# Logistic curves
#  All stocks
all.timing.logistic <- read.csv(file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))

#  Yukon Stock only
can.timing.logistic <- read.csv(file = file.path(dir.output,"logistic curve parameters CAN Chinook 1995_2022.csv"))

# Normal Curve
#  All stocks
all.timing.normal <- read.csv(file = file.path(dir.output,"normal curve parameters All Chinook 1995_2022.csv"))

# Yukon Stock only
can.timing.normal <- read.csv(file = file.path(dir.output,"normal curve parameters CAN Chinook 1995_2022.csv"))

# Create a df for plotting
all.timing.logistic$method <- as.factor("Logistic")
can.timing.logistic$method <- as.factor("Logistic")

all.timing.logistic$Chinook <- as.factor("All")
can.timing.logistic$Chinook <- as.factor("Can")

all.timing.normal$method <- as.factor("Normal")
can.timing.normal$method <- as.factor("Normal")

all.timing.normal$Chinook <- as.factor("All")
can.timing.normal$Chinook <- as.factor("Can")

time.df.log <- rbind(all.timing.logistic,
                     can.timing.logistic,
                     all.timing.normal,
                     can.timing.normal)


# Plots ##################################################################

# Pivot wide
test <- pivot_wider(time.df.log, names_from  = c("Chinook","method"), values_from = c("mid", "sd", "alpha","year"))

# Years used for timing comparison
yearPSS <- c(1995,1996:2022)

# Normal vs Logistic All Chinook
plot(all.timing.logistic$mid ,all.timing.normal$mid,
     xlab = "Logistic Curve Midpoint",
     ylab = "Normal Curve Midpoint",
     main = "Logistic vs Normal Midpoint")
abline(0,1, lty = 2)
legend('bottomright', legend=c("1:1 Line"), 
       col=c('gray','red','blue',"green"), lty=c(2),pch = c())
# text(x = all.timing.logistic$mid, y = all.timing.normal$mid+0.5, labels = yearPSS)

# Normal vs Logistic Can Chinook
plot(can.timing.logistic$mid ,can.timing.normal$mid,
     xlab = "Logistic Curve Midpoint",
     ylab = "Normal Curve Midpoint",
     main = "Logistic vs Normal Midpoint Canadian Chinook")
abline(0,1, lty = 2)
legend('bottomright', legend=c("1:1 Line"), 
       col=c('gray','red','blue',"green"), lty=c(2),pch = c())
# text(x = all.timing.logistic$mid, y = all.timing.normal$mid+0.5, labels = yearPSS)

# Can vs All logistic curve timing
plot(x = can.timing.logistic$mid , y = all.timing.logistic$mid,
     xlab = " Canadian Chinook Midpoint",
     ylab = "All Chinook Midpoint",
     main = "All vs Canadian Chinook Midpoint Timing",
     cex= 1.5,
     cex.main = 1.5,
     cex.lab = 1.5, 
     xlim = c(165,180),
     ylim = c(165,180),
     pch = 21,
     bg = "green")
abline(0,1, lty = 2, lwd = 2, col = "red")
legend('bottomright', 
       legend=c("1:1 Line", "Logistic MP", "Normal MP"), 
       col=c("red", "green", "blue"),
       lty=c(2,NA,NA),
       pch = c(NA,21,21),
       bg = c(NA, "green", "blue"))

points(x = can.timing.normal$mid, y = all.timing.normal$mid,
       cex= 1.5,
       pch = 21,
       bg = "Blue")


# Normal vs Logistic  
plot(x = all.timing.normal$mid , y = all.timing.logistic$mid,
     xlab = " Normal Curve Midpoint",
     ylab = "Logistic Curve Midpoint",
     main = "Logistic vs Normal MP",
     cex= 1.5,
     cex.main = 1.5,
     cex.lab = 1.5, 
     xlim = c(165,180),
     ylim = c(165,180),
     pch = 21,
     bg = "Green")
abline(0,1, lty = 2, lwd = 2, col = "red")
legend('bottomright', 
       legend=c("1:1 Line", "All Chinook MP", "Canadian Chinook MP"), 
       col=c("red", "green", "blue"),
       lty=c(2,NA,NA),
       pch = c(NA,21,21),
       bg = c(NA, "green", "blue"))

points(x = can.timing.normal$mid, y = can.timing.logistic$mid,
       cex= 1.5,
       pch = 21,
       bg = "Blue")


#

midPoint<- PSS_hist %>% 
              group_by(Year) %>% 
              summarise(Cumcount = cumsum(count), .groups = ) %>% 
              as.data.frame()

full_mid <- cbind(PSS_hist,midPoint)

full_mid <- full_mid[,-4]

full_mid_no0 <- full_mid[full_mid$Cumcount != 0,]

for (y in c(1995,1997:2021)) {
   
  # y = 1997

      a <- full_mid_no0[full_mid_no0$Year == y &
               full_mid_no0$Cumcount < max(full_mid_no0$Cumcount[full_mid_no0$Year == y]), ]

      b <- full_mid_no0[full_mid_no0$Year == y &
                    full_mid_no0$Cumcount == max(full_mid_no0$Cumcount[full_mid_no0$Year == y]), ][1,]
  
      
      if(y == 1995){temp.mp.df <- rbind(a,b)}else{
      
      temp.mp.df <- rbind(temp.mp.df,a,b)}
  
}


maxYear <- full_mid %>% group_by(Year) %>% summarise(maxcount = max(Cumcount)) %>% as.data.frame()

full_mid <- left_join(temp.mp.df,maxYear)

full_mid$percent <- (full_mid$Cumcount/full_mid$maxcount)*100

head(full_mid)

colnames(all.timing.logistic)[5] <- "Year"

full_mid <- left_join(full_mid, all.timing.logistic)

full_mid$dummy <- 50

full_mid$time.deviation <- full_mid$mid - mean(all.timing.logistic$mid)

full_mid$dev.code <- NA
for (i in 1:length(full_mid$Day)) {
  

full_mid$dev.code[i] <- if(full_mid$time.deviation[i] < 0){"Early"}else{if(full_mid$time.deviation[i] > 0){"Late"}else{"Avg"}}

}

# Plot of cumulative PSS passage and midpoint colored by early or late
ggplot(full_mid[full_mid$Day <= 215,],
       aes(x = Day,
           y = percent))+
  geom_point(aes(x = mid, y = dummy, color = dev.code), size = 3)+
  geom_line()+
  geom_vline(xintercept = (mean(all.timing.logistic$mid)),
             linetype = 5,
             show.legend = T)+
  # geom_vline(aes(xintercept = mid))+
  # geom_point(aes(x = days))+
  facet_wrap(~Year)+
  scale_x_continuous(breaks = c(160,180,200), labels = c(
    "June 9",
    "June 29",
    "July 19"
  ))+
  xlab("Day")+
  ylab("Cummulative PSS Chinook Salmon (%)")+
  theme(text = element_text(size = 17),
        axis.text.x = element_text(angle = 90))


ggplot(full_mid, aes(x = Year, y = time.deviation, color = dev.code))+
  geom_point(size = 3)+
  geom_hline(yintercept = 0,
             linetype = 5,
             color = "red",
             size = 3)+
  scale_y_continuous(name = "Days", limits = c(-7,7), n.breaks = 14)
  
  




curves <-rbind(all.timing.logistic,all.timing.normal)
ggplot(all.timing.logistic, aes(x = year, y = mid, color = method))+
  geom_line()+
  geom_hline(aes(yintercept = mean(mid)))+
  facet_wrap(~method)



##################################################################################################








