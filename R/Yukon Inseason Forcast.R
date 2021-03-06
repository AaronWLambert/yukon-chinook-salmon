#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Version: 1.0
# Purpose: To generate predictions from an integrated Bayesian inseason assessment
#           model.
#
#   1) Read in data
#   2) Use linear model to determine relationship between Canadian PSS passage and
#      end-of-season reconstructed runsize 
#   3) Call to stan model to generate inseason projection
#
#
#=================================================================================
#NOTES:
# This script is the general working script that is used to run single iterations 
#  of the model with any current Stan file.
#
# 
# Next steps: 
# 
#
#=================================================================================
require(rstan)
require(bayesplot)
# require(rstanarm)
require(tidyverse)
require(mgcv)
require(ggthemes)
require(viridis)
require(shinystan)
require(lubridate)
require(reshape2)
require(dplyr)
require(ggthemes)
require(tidybayes)


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
# NOTE!!!! See use "Dataframe preprocess.R" to get  most recent version.
#  If not using Dataframe preprocess.R, uncomment out the following to load relevant PSS historical data
# PSS_hist <- readRDS(file = file.path(dir.data,"pss adfg 27April22.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
# GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - 2021)
# pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
pf_hist <- readRDS(file.path(dir.data,"inv_var_weighted_forcast_v3_June072022.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs) and days not covered by GSI
# GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year unadj 27April22.RDS"))


# Control Section ######################
# This is where the user can input dates, stan model, and stan model controls.
model.version <- "3.2"
# Range of years 1995 to 2021
myYear <- 2022

# Range of days 152 -243
myDay <- 165

# MCMC Parameters
n.chains <- 4
n.iter <- 10000;#5e4
n.thin <- 2

# Start Years for Predictors
startYearPF <- 2007
startYearGSI <- 2005
startYearPSS <-1997

# Wont typicaly change;
# day 152 = June 1
startDayPSS <- 148

# Preseason Forecast######################

# Current Preseason forecast Can Origin
pf <- log(pf_hist$Weighted_Forecast[pf_hist$Year == myYear])

# Vector of historical preseason forecasts for to compute sd for prior in stan model
PF_vect <- pf_hist$Weighted_Forecast[pf_hist$Year<=2021]

# EOS reconstructed runsize for historic years that have a preseason forecast 
EOS <- CAN_hist$can.mean[CAN_hist$Year >= 2007 & 
                          CAN_hist$Year != 2011 & 
                          CAN_hist$Year != 2012]

# SD for preseason forecast prior
pf_sigma <- sd(log(PF_vect/EOS))



# Metadata###############################
# This section is for data pre-processing

# Vector of historic PSS years up to myYear-1
yearPSS <- unique(PSS_hist$Year[PSS_hist$Year<=(myYear-1)])
# Number of years used in model
n_yearPSS <- length(yearPSS)

# # PSS days included up to myDay. 151 is subtracted for use in logistic model
# #  I.e., June 1 = 1, June 2 = 2 ....
dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])-147

# Number of days used
n_dayPSS <- length(dayPSS)


# PSS daily passage estimate for days up to myDay for the year of interest (myYear)
curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay])
# Number of days used
n_curr_PSS <- length(curr_PSS)

# TAKE ME OUT##
cum_curr_PSS <- cumsum(curr_PSS)
# # Historical  Canadian Origin PF (2013 - current)
# histPF <- pf_hist$Weighted_Forecast[pf_hist$Year <= (myYear-1)]
# names(histPF) <- (pf_hist$Year[pf_hist$Year <= (myYear-1)]) 
# n_histPF <-length(histPF)

# End of season Canadian abundance estimates for 1995-2020 excluding 1996 (not included in data) 
#  This is what cum PSS is regressed against in Stan Model and years up to myYear-1 are included
totalEOS <- CAN_hist$can.mean[CAN_hist$Year <= (myYear-1) & CAN_hist$Year != 1996]

# Name the elements for accounting purposes
names(totalEOS) <- CAN_hist$Year[CAN_hist$Year <= (myYear-1)& CAN_hist$Year != 1996]
#Number of 
n_totalEOS <- length(totalEOS)

# Create matrix with PSS daily Chinook passage from day 152 to my day and 1995
# to my year
#  *Note that 1996 is missing*


# Create matrix with dimensions myday and myYear
PSS_mat <- matrix(nrow = (length(startDayPSS:myDay )),
                  #start 1996 to account for missing year 1996
                  ncol = length(1996:(myYear-1))) 

# Give names to matrix for accounting
colnames(PSS_mat) <- c("1995", startYearPSS:(myYear-1))
rownames(PSS_mat) <- c(startDayPSS:(myDay))

# Create vector of counts from relevant years and days to input into matrix
(count_vect <- PSS_hist$count[PSS_hist$Year<= (myYear-1) &
                                PSS_hist$Day <= (myDay)])

# The number of observations in count_vector for PSS historic counts
(n_hist_counts <- length(count_vect))

# Use loop to populate matrix with counts by days
# Set counter
counter <-1

for (y in 1:((myYear-1)-1995)){
  for (d in 1:(length(startDayPSS:myDay))) {
    
    
    PSS_mat[d,y] <- count_vect[counter]
    counter = counter+1
    
  }
}

# Cumulative count matrix for historic years
PSS_cum_hist_mat <- apply(X = PSS_mat,MARGIN = 2,FUN = cumsum)

# Vector for storing cumPSS for plotting
cumPSS <- vector(length = n_yearPSS)

names(cumPSS)<- c(yearPSS)
# Cummulative PSS Counts for plotting
for (i in 1:n_yearPSS) {
  
  cumPSS[i]<- sum(PSS_mat[,i])
}

# Vector containing avg GSI proportions for each day ACROSS ALL YEARS (for versions 3.1,3.2,3.3) ####
meanGSI_vect <- vector(length = length(startDayPSS:myDay))
names(meanGSI_vect) <- c(startDayPSS:myDay)
sdGSI_vect <- vector(length = length(startDayPSS:myDay))
names(sdGSI_vect) <- c(startDayPSS:myDay)

counter <- 1
for (d in startDayPSS:myDay) {
  # d = 175
  meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
  counter <- counter+1
}

# No GSI data for day 148 and 149. Use day 150 for these values
meanGSI_vect[1:2] <- meanGSI_vect[3]
sdGSI_vect[1:2] <- sdGSI_vect[3]

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
GSI_avg_vect <- GSI_avg[1:length(startDayPSS:myDay)]

GSI_sd <- c(rep(GSI$stratumSD[GSI$stratum == 1], 
                times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
            rep(GSI$stratumSD[GSI$stratum == 2], 
                times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
            rep(GSI$stratumSD[GSI$stratum == 3],
                times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))
GSI_sd_vect <- GSI_sd[1:length(startDayPSS:myDay)]

# Create matrix with dimensions myday and myYear
PSS_mat_adj <- matrix(nrow = (length(startDayPSS:myDay )),
                      #start 1996 to account for missing year 1996
                      ncol = length(1996:(myYear-1))) 

# Give names to matrix
colnames(PSS_mat_adj) <- c("1995", startYearPSS:(myYear-1))
rownames(PSS_mat_adj) <- c(startDayPSS:(myDay))


for (y in 1:n_yearPSS) {
  for (d in 1:n_dayPSS) {
    # y = 1
    # d = 15
    PSS_mat_adj[d,y] <- PSS_mat[d,y]*meanGSI_vect[d] 
    
  }
  
}

# Vector for storing cumPSS for plotting
cumPSS_adj <- vector(length = n_yearPSS)

names(cumPSS_adj)<- c(yearPSS)
# Cummulative PSS Counts for plotting
for (i in 1:n_yearPSS) {
  
  cumPSS_adj[i]<- sum(PSS_mat_adj[,i])
}

adj_curr_PSS <- vector(length = n_dayPSS)
# for (d in 1:n_dayPSS) {
#   
# 
# adj_curr_PSS[d] <- curr_PSS[d]*meanGSI_vect[d]
# }
# Stan Model Call ######################################


# inits <- function(){ list(phi = runif(n = 1,0,1),
#               # propCAN_1 = runif(n=1,.1,.9),
#               propCAN_logit = runif(n_dayPSS,0,1),
#               # alpha = 500,
#               # beta = 4,
#               # sigma = 5,
#               ln_predPSS =  runif(23, 2,5))}

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
                        "n_hist_counts"=n_hist_counts,
                        "Pf" = pf, 
                        "Pf_sigma" = pf_sigma,
                        "meanpropCAN" = meanGSI_vect,
                        "sd_meanpropCAN" = sdGSI_vect,
                        "mean_adj_PSS_mat" = PSS_mat_adj,
                        "mean_adj_curr_PSS"= adj_curr_PSS,
                        "PSS_cum_hist_mat" = PSS_cum_hist_mat,
                        "cum_curr_PSS" = cum_curr_PSS,
                        "GSI_mean"= GSI_avg_vect,
                        "GSI_sd" = GSI_sd_vect),
            # init = inits,
            chains = n.chains,
            iter = n.iter, 
            thin = n.thin, 
            # cores = n.chains,
            cores = mc.cores,
            control = list(max_treedepth = 25, adapt_delta = 0.99),
            verbose = F
              
            )



# Figures Section ################################

traceplot(object = fit, c("alpha","beta","sigma", "RunSize",
                          "phi","propCAN[1]","propCAN[4]"))
# shinystan::launch_shinystan(as.shinystan(fit))

mcmc_pairs(fit, pars = c("alpha","sigma", "beta", "phi"))
# Extract parameter estimates
pars <- rstan::extract(fit)
names(pars)

# Plots for looking at predPSS fit

# Get the quantiles for predPSS from model output
quant.predPSS <- apply(X = pars$predPSS, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))

# Prder for polygon
ord <- order(cumPSS)

# Call plot


plot(x = cumPSS, y = totalEOS,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = paste("Cumulative PSS Passage Day ", myDay, "\n Years 1995-",myYear-1 ),
     ylab = paste("Totoal Reconstructed Canadian Passage 1995-",myYear
                  -1),
     # main = paste("Predicted PSS Fit","1995 to", 
                  # myYear-1,"\n Day 152 to",
                  # myDay,"\n Model Version",model.version),
     xlim = c(0,max(cumPSS)),
     ylim = c(min(cumPSS), max(quant.predPSS)))
polygon(x = c(cumPSS[ord], rev(cumPSS[ord])),
        y = c(quant.predPSS[1,ord],rev(quant.predPSS[5,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS[ord], rev(cumPSS[ord])),
        y = c(quant.predPSS[2,ord],rev(quant.predPSS[4,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS[ord], quant.predPSS[3,ord], col = "red", lw = 2)

# zz<-recordPlot()
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
  # theme_classic()+
  # geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
  #            color = "End of Season Runsize"),
  #            linetype = 2,
  #            size = 1)+
  # xlim(c(0,300000))+
  coord_cartesian(xlim = c(0,250))+
  # ylim(c(0,500))
  labs(x = "1000's of Chinook Salmon", fill = "Parameters", y = "Probability Density")+
  # scale_fill_discrete(name = "Parameters", 
                      # labels = c( "Preseason Forecast (Prior)", 
                                  # "PSS Prediction","Runsize"))+
  scale_x_continuous()+
  scale_fill_colorblind(name = "", 
                        labels = c( "Preseason Forecast (Prior)", 
                                    "PSS Prediction","Runsize"))+
  # theme_tidybayes(element_text(size = 20))+
  theme(text = element_text(size = 20, family = "serif"),
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")


pars$propCAN
dim(pars$propCAN)

test <- apply(pars$propCAN , MARGIN = 2, median)
plotDF <- data.frame("median_proportions"=test,"day"= n_dayPSS)


newMat <- pars$propCAN
colnames(newMat)<- 152:myDay
propDF <- as.data.frame(newMat)

longDF <- melt(propDF)


ggplot(longDF, aes(x = variable, y = value) )+
  geom_boxplot(aes(fill = "red"))+
  theme_clean()+
  theme(axis.text.x  = element_text(angle = 90))+
  labs(x = "Day of Year", y = "Estimated Canadian Proportion")


hist(rbeta(n = 1000,pars$paramA[1,2], pars$paramB[1,2]))

pars$
#################################################################################
# Create matrix with dimensions myday and myYear
PSS_mat_adj <- matrix(nrow = (length(startDayPSS:myDay )),
                  #start 1996 to account for missing year 1996
                  ncol = length(1996:(myYear-1))) 

# Give names to matrix
colnames(PSS_mat_adj) <- c("1995", startYearPSS:(myYear-1))
rownames(PSS_mat_adj) <- c(startDayPSS:(myDay))

# Make sure it works...
dim(PSS_mat_adj)


for (y in 1:n_yearPSS) {
  for (d in 1:n_dayPSS) {
    # y = 1
    # d = 15
    PSS_mat_adj[d,y] <- PSS_mat[d,y]*meanGSI_vect[d] 
    
  }
  
}

# Vector for storing cumPSS for plotting
cumPSS_adj <- vector(length = n_yearPSS)

names(cumPSS_adj)<- c(yearPSS)
# Cummulative PSS Counts for plotting
for (i in 1:n_yearPSS) {
  
  cumPSS_adj[i]<- sum(PSS_mat_adj[,i])
}
plot(x = cumPSS_adj, y = totalEOS,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = paste("Cumulative PSS Passage Day ", myDay, "\n Years 1995-",myYear-1 ),
     ylab = paste("Totoal Reconstructed Canadian Passage 1995-",myYear
                  -1),
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     xlim = c(0,max(cumPSS_adj)),
     ylim = c(min(cumPSS_adj), max(quant.predPSS)))
polygon(x = c(cumPSS_adj[ord], rev(cumPSS_adj[ord])),
        y = c(quant.predPSS[1,ord],rev(quant.predPSS[5,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS_adj[ord], rev(cumPSS_adj[ord])),
        y = c(quant.predPSS[2,ord],rev(quant.predPSS[4,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS_adj[ord], quant.predPSS[3,ord], col = "red", lw = 2)

