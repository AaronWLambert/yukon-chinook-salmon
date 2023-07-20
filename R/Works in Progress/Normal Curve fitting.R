#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 06.22.22
# Version: 4.0
# Purpose: To generate predictions from an integrated Bayesian inseason assessment
#           model.
#
#   1) Read in data
#   2) Preprocess the data for model
#   3) Fit a Gaussian curve to past years PSS passage and estimate mu, sd, and alpha 
#       parameters for use in model priors.
#   3) Call to Stan model to generate inseason projection
#
#
#=================================================================================
#NOTES:
# This script is the general working script that is used to run single iterations 
#  of the model with any current Stan file.
library(bbmle)
library(mvtnorm)
library(rstan)
library(bayesplot)
library(tidyverse)
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
# PSS_hist <- readRDS(file = file.path(dir.data,"PSS passage 22Jun22.RDS"))

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
model.version <- "5.0"
# Range of years 1995 to 2021
myYear <- 2017

# Range of days 152 -243
myDay <- 180

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
dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])

# Number of days used
n_dayPSS <- length(dayPSS)


# PSS daily passage estimate for days up to myDay for the year of interest (myYear)
curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay])
# Number of days used
n_curr_PSS <- length(curr_PSS)

# Cumulative current PSS counts
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

# Cumulative PSS matrix
cum_PSS_mat <-  apply(PSS_mat, 2, cumsum)



# PSS martix of complete historic counts ***********
# Create matrix with dimensions myday and myYear
PSS_mat_all <- matrix(nrow = (length(startDayPSS:252 )),
                  #start 1996 to account for missing year 1996
                  ncol = length(1996:(myYear-1))) 

# Give names to matrix for accounting
colnames(PSS_mat_all) <- c("1995", startYearPSS:(myYear-1))
rownames(PSS_mat_all) <- c(startDayPSS:252)

# Create vector of counts from relevant years and days to input into matrix
(count_vect_all <- PSS_hist$count[PSS_hist$Year<= (myYear-1) &
                                PSS_hist$Day <= 252])

# The number of observations in count_vector for PSS historic counts
(n_hist_counts_all <- length(count_vect_all))

# Use loop to populate matrix with counts by days
# Set counter
counter <-1

for (y in 1:((myYear-1)-1995)){
  for (d in 1:(length(startDayPSS:252))) {
    
    
    PSS_mat_all[d,y] <- count_vect_all[counter]
    counter = counter+1
    
  }
}

# Vector containing avg GSI proportions for each day ACROSS ALL Years
#  (for versions 3.0,3.1,3.2,3.3) ####
meanGSI_vect <- vector(length = length(startDayPSS:252))
names(meanGSI_vect) <- c(startDayPSS:252)
sdGSI_vect <- vector(length = length(startDayPSS:252))
names(sdGSI_vect) <- c(startDayPSS:252)

counter <- 1
for (d in startDayPSS:252) {
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



pss_days_all <- 148:252
n_pss_days_all <- length(pss_days_all)



# init_fn <- function(chain_id=1:n.chains){
# 
#   list("ps_alpha_curr" = runif(1,20000,100000),
# "ps_mu_curr" = runif(1,160,186),
# "ps_sd_curr" = runif(1,10,20))
# 
# }
# 
# init_fn(1:4)
# 
# init_ll <- lapply(1, function(id)init_fn(chain_id = id))


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
                        "cum_PSS_mat" = cum_PSS_mat,
                        "cum_curr_PSS" = cum_curr_PSS,
                        "GSI_mean"= GSI_avg_vect,
                        "GSI_sd" = GSI_sd_vect,
                        # "ps_mu"= ps_mu,
                        # "ps_sd"= ps_sd,
                        "ps_alpha"= ps_alpha,
                        "pss_days_all" = pss_days_all,
                        "n_pss_days_all" = n_pss_days_all,
                        "PSS_mat_adj" = PSS_mat_adj,
                        "myDay" = myDay,
                        "PSS_mat_all" = PSS_mat_all,
                        "ps_m"= ps_m,
                        "ps_s"= ps_s),
            # init = init_ll,
            chains = 1,
            iter = 1000, 
            thin = n.thin, 
            # cores = n.chains,
            cores = mc.cores,
            control = list(max_treedepth = 25, adapt_delta = 0.99),
            verbose = T
            
)

# saveRDS(object = fit, file = file.path(dir.output,"fit4.2_13July22.RDS"))

# Save to PDF
# pdf(file = file.path(dir.output,"4.0_plots_13jul22_restrictedMP.pdf"))

traceplot(object = fit, c("ps_mid_curr", "ps_shape_curr", "ps_alpha_curr","alpha", 
                          "beta","ln_RunSize"))

# shinystan::launch_shinystan(as.shinystan(fit))
# mcmc_pairs(fit, pars = c("ps_mu_curr", "ps_sd_curr", "ps_alpha_curr"))

# Extract parameter estimates
pars <- rstan::extract(fit)

# # Limits
x.lim <- c(145,252)
y.lim <- c(0,max(apply(pars$ps_pred_curr, 2, quantile, probs=0.975), curr_PSS))

# Plot Observed
plot(x=dayPSS,
     y=curr_PSS, 
     type='p', pch=21, 
     col='gray',
     xlim=x.lim,
     ylim=y.lim,
     xlab='',
     ylab='',
     main="PSS: Current Year", 
     xaxt='n')

axis(side=1, at=seq(from=151, to=206, by=5), labels = FALSE, col='darkgray')
axis(side=1, at=seq(from=161, to=201, by=10),
     labels=c("June-10","June-20","June-30","July-10","July-20"))
segments(x0=dayPSS, y0=rep(0, length(dayPSS)), x1=dayPSS, y1=curr_PSS, col='black', lwd=2)
points(x=dayPSS, y=curr_PSS, type='p', pch=21, bg='gray')

# Plot Predicted
polygon(x=c(pss_days_all,rev(pss_days_all)), y=c(apply(pars$ps_pred_curr, 2, quantile, probs=0.025),
                                       rev(apply(pars$ps_pred_curr, 2, quantile, probs=0.975))), border=FALSE, col=rgb(1,0,0,0.2))
polygon(x=c(pss_days_all,rev(pss_days_all)), y=c(apply(pars$ps_pred_curr, 2, quantile, probs=0.25),
                                       rev(apply(pars$ps_pred_curr, 2, quantile, probs=0.75))), border=FALSE, col=rgb(1,0,0,0.2))
lines(x=pss_days_all, y=apply(pars$ps_pred_curr, 2, median), col=rgb(1,0,0, alpha=0.5), lwd=3)
abline(v=mean(ps_mu), col='blue', lty=2, lwd=2)

# Labels
mtext(paste(" Day,",myDay,'of',myYear), side=3, outer=TRUE, line=-4, font=2)
mtext('Day of Year', side=1, outer=TRUE, line=0.5, font=2)
mtext('Reconstructed Total Canadian-origin Chinook Salmon', side=2, outer=TRUE, line=-2, font=2)
legend('topright', legend=c('Obs.','Pred.', 'Avg. Peak'), col=c('black','red','blue'), lty=c(1,1,2))

# Historical
# y <- 3
for(y in 1:n_yearPSS) {
  year <- yearPSS[y]
  
  # Limits
  x.lim <- c(145,255)
  y.lim <- c(0,max(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.975), PSS_mat[,y]))
  
  # Plot Observed
  plot(x=dayPSS, y=PSS_mat[,y], type='p', pch=21, col='gray',
       xlim=x.lim,
       ylim=y.lim,
       xlab='', ylab='',
       main=year, xaxt='n')
  
  axis(side=1, at=seq(from=151, to=206, by=5), labels = FALSE, col='darkgray')
  axis(side=1, at=seq(from=161, to=201, by=10),
       labels=c("June-10","June-20","June-30","July-10","July-20"))
  segments(x0=dayPSS, y0=rep(0, length(dayPSS)), x1=dayPSS, y1=PSS_mat[,y], col='black', lwd=2)
  points(x=dayPSS, y=PSS_mat[,y], type='p', pch=21, bg='gray')
  
  # Plot Predicted
  polygon(x=c(pss_days_all,rev(pss_days_all)), y=c(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.025),
                                         rev(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.975))), border=FALSE, col=rgb(1,0,0,0.2))
  polygon(x=c(pss_days_all,rev(pss_days_all)), y=c(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.25),
                                         rev(apply(pars$ps_pred_hist[,,y], 2, quantile, probs=0.75))), border=FALSE, col=rgb(1,0,0,0.2))
  lines(x=pss_days_all, y=apply(pars$ps_pred_hist[,,y], 2, median), col=rgb(1,0,0, alpha=0.5), lwd=3)
  abline(v=mean(ps_mu), col='blue', lty=2, lwd=2)
  
  # Labels
  mtext(paste("Day",myDay,'of',myYear), side=3, outer=TRUE, line=-5, font=2)
  mtext('Day of Year', side=1, outer=TRUE, line=-2.5, font=2)
  mtext('Can-origin PSS Daily Chinook Passage', side=2, outer=TRUE, line=-1.7, font=2)
  legend('topright', legend=c('Obs.','Pred.', 'Avg. Mid-point'), col=c('black','red','blue'), lty=c(1,1,2))
}


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
  scale_color_discrete(name = "")

# Get the quantiles for predPSS from model output
quant.predPSS <- apply(X = pars$predPSS, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# quant.predPSS <- quantile(pars$post_curr_predPSS, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

cumPSS <- apply(PSS_mat_all,MARGIN = 2,sum)
# Prder for polygon
ord <- order(cumPSS)


adj_cumPSS <- vector(length = n_yearPSS)

for (i in 1:n_yearPSS) {
  adj_cumPSS[i]<-sum(PSS_mat_adj[,i])
}

# Call plot
plot(x = cumPSS, y = totalEOS,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = paste("Cumulative PSS Passage Day ", myDay, "\n Years 1995-",myYear-1 ),
     ylab = paste("Totoal Canadian-origin Chinook 1995-",myYear
                  -1),
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     xlim = c(0,max(cumPSS)),
     ylim = c(min(quant.predPSS), max(quant.predPSS)))
polygon(x = c(cumPSS[ord], rev(cumPSS[ord])),
        y = c(quant.predPSS[1,ord],rev(quant.predPSS[5,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS[ord], rev(cumPSS[ord])),
        y = c(quant.predPSS[2,ord],rev(quant.predPSS[4,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS[ord], quant.predPSS[3,ord], col = "red", lw = 2)

#
quantile(pars$RunSize, probs = c(.025,.1,.5,.9,.975))

dev.off()
boxplot(pars$propCAN)
