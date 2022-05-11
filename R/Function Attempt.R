#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 10.14.21
# Version: 1.0
# Purpose: To generate predictions from an integrated Bayesian inseason assessment
#           model.
#
#   
#
#
#=================================================================================
#NOTES:
# This is a script for experimenting with how to create a function for the model
#
# 
# Next steps: 
# 
#
#=================================================================================
require(BEST)
require(rstan)
require(bayesplot)
# require(rstanarm)
require(tidyverse)
require(mgcv)
require(ggthemes)
require(viridis)
require(shinystan)
require(lubridate)
require(beepr)
require(reshape2)
require(dplyr)
require(ggthemes)
require(ggpubr)
require(gridExtra)
require(tidybayes)
require(loo)
require(here)
require(forcats)

# Parralize for optimum model run time

rstan_options(auto_write = TRUE)

mc.cores = parallel::detectCores()

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




InSeasonProjection <- function(model.version = "2.0",
                               myYear ,
                               myDay , 
                               n.chains = 4,
                               n.iter = 30000, 
                               n.thin = 2,
                               CAN_hist,
                               PSS_hist,
                               GSI_mean,
                               pf_hist){
# Control Section ######################
# model.version <- version
# Range of years $$$$ to 2021
# myYear <- my.Year

# Range of days 152 -243
# myDay <- my.Day

# MCMC Parameters
# n.chains <- n.chains
# n.iter <- n.iter;#5e4
# n.thin <- n.thin

# Start Years for Predictors
startYearPF <- 2007
startYearGSI <- 2005
startYearPSS <-1997

# Wont typicaly change;
# day 152 = June 1
startDayPSS <-152
# Preseason Forcast######################

# Current Preseason forecast Can Origin
pf <- log(pf_hist$Weighted_Forecast[pf_hist$Year == myYear])
# pf_sigma <- (0.6)
PF_vect <- pf_hist$Weighted_Forecast[pf_hist$Year<=2019]
EOS <- CAN_hist$can.mean[CAN_hist$Year >= 2007 &
                          CAN_hist$Year != 2011 &
                          CAN_hist$Year != 2012]

pf_sigma <- sd(log(EOS/PF_vect))



# Metadata###############################

# PSS years used
yearPSS <- unique(PSS_hist$Year[PSS_hist$Year<=(myYear-1)])
n_yearPSS <- length(yearPSS)

# PSS Days 
dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])
n_dayPSS <- length(dayPSS)

# PSS Counts for myYear
curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay])
n_curr_PSS <- length(curr_PSS)

# Historic preseason forcast years
# yearPF <- pf_hist$Year[pf_hist$Year <= (myYear-1)]
# n_yearPF <- length(yearPF)

# Historical  Canadian Origin PF (2013 - current)
histPF <- pf_hist$Weighted_Forecast[pf_hist$Year <= (myYear-1)]
names(histPF) <- (pf_hist$Year[pf_hist$Year <= (myYear-1)]) 
n_histPF <-length(histPF)

# End of season Canadian counts 1995-2020 excluding 1996 for comparison to PSS
totalEOS <- CAN_hist$can.mean[CAN_hist$Year <= (myYear-1) &
                               CAN_hist$Year != 1996]
names(totalEOS) <- CAN_hist$Year[CAN_hist$Year <= (myYear-1)
                                          & CAN_hist$Year != 1996]
n_totalEOS <- length(totalEOS)

# Create matrix with PSS daily counts from day 152 to my day and 1995
# to my year
#  *Note that 1996 is missing*


# Create matrix with dimensions myday and myYear
PSS_mat <- matrix(nrow = (length(startDayPSS:myDay )),
                  #start 1996 to account for missing year 1996
                  ncol = length(1996:(myYear-1))) 

# Give names to matrix
colnames(PSS_mat) <- c("1995", startYearPSS:(myYear-1))
rownames(PSS_mat) <- c(startDayPSS:(myDay))

# Make sure it works...
dim(PSS_mat)
str(PSS_hist)

# Create vector of counts from relevant years and days
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

# Vector for storing cumPSS for plotting
cumPSS <- vector(length = n_yearPSS)

names(cumPSS)<- c(yearPSS)
# Cummulative PSS Counts for plotting
for (i in 1:n_yearPSS) {
  
  cumPSS[i]<- sum(PSS_mat[,i])
}

# Stan Model Call ######################################


fit<- stan(file = file.path(dir.stan,paste("Yukon Inseason Forecast ",
                                            model.version ,".stan", sep = "")),
            data = list("PSS_mat"=PSS_mat,
                        "n_histPF" =n_histPF, 
                        "histPF"=histPF,
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
                        "Pf_sigma" = pf_sigma),
            chains = n.chains,
            iter = n.iter, 
            thin = n.thin, 
            cores = n.chains
)





 

# Figures Section ################################


# shinystan::launch_shinystan(as.shinystan(fit))

# Extract parameter estimates
pars <- rstan::extract(fit)
# names(pars)

# mcmc_intervals(fit, pars = c("prior_pf","RunSize", "curr_predPSS"), prob = 0.95)

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
     xlab = "Cumulative PSS Passage",
     ylab = "Totoal Reconstructed Canadian Passage",
     main = paste("Predicted PSS Fit","1995 to", 
                  myYear-1,"\n Day 152 to",
                  myDay,"\n Model Version",model.version),
     xlim = c(0,max(cumPSS)),
     ylim = c(min(cumPSS), max(quant.predPSS)))
polygon(x = c(cumPSS[ord], rev(cumPSS[ord])),
        y = c(quant.predPSS[1,ord],rev(quant.predPSS[5,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(cumPSS[ord], rev(cumPSS[ord])),
        y = c(quant.predPSS[2,ord],rev(quant.predPSS[4,ord])),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = cumPSS[ord], quant.predPSS[3,ord], col = "red", lw = 2)
zz<-recordPlot()

# Extract estimates for prior, likelihood, and posterior 
df.run <- data.frame("par" = "Runsize", "value" = pars$RunSize)
df.prior <- data.frame("par" = "Prior", "value" = pars$prior_pf)
df.pssPred <- data.frame("par" = "PSS_Pred", "value" = pars$curr_predPSS)
df.postPredPss <- data.frame("par" = "PSS_post_pred", "value" = pars$post_curr_predPSS)

# Bind into one dataframe with 2 columns
df.comb <- rbind(df.prior,df.postPredPss,df.run)
# df.comb$par <- levels(df.comb$par, levels = c("Prior","PSS_post_pred","Runsize"))
# Density plots comparing pf, linear prediction, and posterior estimate
testPlot<-ggplot(df.comb, aes(x = value/1000, fill = par))+
  geom_density(alpha = .65)+
  # ggtitle(paste("Year =",myYear,"\n Day =",myDay))+
  
  # theme_classic()+
  geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
             col = "True Runsize"),
             linetype = 2,
             size = 1.5)+
  # ylab("Relative Probability")+
  # xlab("1000's of Chinook Salmon")+
  xlab("")+
  ylab("")+
  coord_cartesian(xlim = c(0,200))+
  # scale_fill_discrete(name = "Parameters", 
                      # labels = c("Preseason Forecast (Prior)", "PSS Prediction","Runsize"))+
  # guides(color = FALSE)+
  scale_fill_colorblind(name = "", 
                        labels = c( "Preseason Forecast (Prior)", 
                                    "PSS Prediction","Runsize Projection"))+
  theme(text = element_text(size = 20, family = "serif"),
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")
#name the list entries
# outputPars <- list(pars)
# outputPlots <- list(testPlot,zz)
outputAll <- list("Pars" = pars,"DensPlot" = testPlot, "PredPlot"= zz)
}

# Plots for paper
output.test.version2<-list()
for(y in 2013){
  for(d in c(153)){
    
    output.test.version2[[paste("Fit Year =",y,"Day =",d)]]<-InSeasonProjection(model.version = "1.0",
                                                                       myYear = y,
                                                                       myDay = d, 
                                                                       n.chains = 4, 
                                                                       CAN_hist = CAN_hist, 
                                                                       pf_hist = pf_hist,
                                                                       PSS_hist = PSS_hist,
                                                                       )
    
  }}

output.test.version2$`Fit Year = 2013 Day = 153`$DensPlot

# Plots for arranging in ggarrange
plotList <- list()
for (i in 1:length(output.test.version2)) {
  plotList[[i]] <- output.test.version2[[i]][[2]]
}
summary(plotList)
ggarrange(plotList[[19]] ,
          plotList[[20]],
          plotList[[21]],
          common.legend = T,
          ncol = 3, nrow = 1,
          labels = c("Jun1-Ver.1","Jun28-Ver.1","Jul15-Ver.1"),
          hjust = -1.5)





plotList[[]]


##################################################### Scratch
output.test[[3]][[2]][[1]]
summary(boo)

ggarrange(boo)
day160Vect<-output.test[[1]][[1]]$RunSize
day175Vect<-output.test[[2]][[1]]$RunSize
day190Vect<-output.test[[3]][[1]]$RunSize

df160 <- data.frame("par" = "Runsize160", "value" = day160Vect)
df175 <- data.frame("par" = "Runsize175", "value" = day175Vect)
df190 <- data.frame("par" = "Runsize190", "value" = day190Vect)

postComb <- rbind(df160, df175, df190)

ggplot(postComb, aes(x = value, fill = par))+
  geom_density(alpha = .5)+
  facet_grid(rows = vars(par))+
  geom_vline(xintercept = CAN_hist$can.mean[CAN_hist$Year == 2017],
             col = "Red",
             linetype = 2,
             size = 5)
?facet_grid

ggarrange(output.test[[1]][[2]],
          output.test[[2]][[2]],
          output.test[[3]][[2]],
          ncol = 1, nrow = 3, 
          common.legend = TRUE,
          # labels = c("June 4 2017 ","July 4 2017","Aug 3 2017","June 4 2018","July 4 2018","Aug 3 2018",
                     # "June 4 2019", "July 4 2017", "Aug 3 2019"),
          font.label = list(size = 12, color = "black", family = "serif"),
          vjust = 1,
          hjust = -1.75)


ggarrange(output.test[[1]][[2]],
          output.test[[2]][[2]],
          output.test[[3]][[2]],
          ncol = 1, nrow = 3, common.legend = TRUE)

par(mfrow = c(3,1))
output.test[[1]][[3]]
output.test[[2]][[3]]
output.test[[3]][[3]]

tryme<-unlist(output.test)

hist(pars$alpha)
hist(pars$beta)
hist(output.test[[1]][[1]]$beta)
?geom_vline
