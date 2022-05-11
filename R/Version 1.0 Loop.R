

# Loop to loop through years and days for stan model version 1.0

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
# Historical Canadian EOS reconstructed passage numbers
# CAN_hist <- readRDS(file.path(dir.data,"/EOS Reconstructed Canadian.RDS"))
CAN_hist <- readRDS(
  file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
# PSS historical for Jun 1st = 1, to Aug 31 = 92
# Years 1995 - 2020
# PSScum_hist <- readRDS(file.path(dir.data,"/clean_counts.rds"))

# Read in PSS non-cummulative with date in date format
PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))

# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))

# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))

# MCMC Parameters
n.chains <- 4
n.iter <- 10000#5e4
n.thin <- 2

# Start Years for Predictors
startYearPF <- 2007
startYearGSI <- 2005
startYearPSS <-1997

# Wont typicaly change;
# day 152 = June 1
startDayPSS <-152


output.version1 <- list()

for(y in 2014:2020){
  for(d in seq(from = 152, to = 222, by = 10)){
    # y = 2015
    # d = 152
    # Preseason Forcast#
    
    # Current Preseason forecast Can Origin
    pf <- log(pf_hist$Mean[pf_hist$Year == y])
    # pf_sigma <- (0.6)
    PF_vect <- pf_hist$Mean[pf_hist$Year <= 2019]
    EOS <- CAN_hist$can.mean[CAN_hist$Year >= 2013]
    
    pf_sigma <- sd(log(EOS/PF_vect))
    
    
    
    # Metadata#
    
    # PSS years used
    yearPSS <- unique(PSS_hist$Year[PSS_hist$Year<=(y-1)])
    n_yearPSS <- length(yearPSS)
    
    # PSS Days 
    dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(d)])
    n_dayPSS <- length(dayPSS)
    
    # PSS Counts for myYear
    curr_PSS <- (PSS_hist$count[PSS_hist$Year == y & PSS_hist$Day <= d])
    n_curr_PSS <- length(curr_PSS)
    
    # Historic preseason forcast years
    # yearPF <- pf_hist$Year[pf_hist$Year <= (myYear-1)]
    # n_yearPF <- length(yearPF)
    
    # Historical  Canadian Origin PF (2013 - current)
    histPF <- pf_hist$Mean[pf_hist$Year <= (y-1)]
    names(histPF) <- (pf_hist$Year[pf_hist$Year <= (y-1)]) 
    n_histPF <-length(histPF)
    
    # End of season Canadian counts 1995-2020
    totalEOS <- CAN_hist$can.mean[CAN_hist$Year <= (y-1)]
    names(totalEOS) <- CAN_hist$Year[CAN_hist$Year <= (y-1)]
    n_totalEOS <- length(totalEOS)
    
    # Create matrix with PSS daily counts from day 152 to my day and 1995
    # to my year
    #  *Note that 1996 is missing*
    
    
    # Create matrix with dimensions myday and myYear
    PSS_mat <- matrix(nrow = (length(startDayPSS:d )),
                      #start 1996 to account for missing year 1996
                      ncol = length(1996:(y-1))) 
    
    # Give names to matrix
    colnames(PSS_mat) <- c("1995", startYearPSS:(y-1))
    rownames(PSS_mat) <- c(startDayPSS:(d))
    
    
    # Create vector of counts from relevant years and days
    (count_vect <- PSS_hist$count[PSS_hist$Year<= (y-1) &
                                  PSS_hist$Day <= (d)])
    
    # The number of observations in count_vector for PSS historic counts
    n_hist_counts <- length(count_vect)
    
    # Use loop to populate matrix with counts by days
    
    # Set counter
    counter <-1
    
    for (i in 1:((y-1)-1995)){
      for (z in 1:(length(startDayPSS:d))) {
        
        
        PSS_mat[z,i] <- count_vect[counter]
        counter = counter+1
        
      }
    }
    
    
    
    
    
    # Stan Model Call #
    
    output.version1[[paste("Fit Year =",y,"Day =",d)]] <-stan(file = file.path(dir.stan,"Yukon Inseason Forcast",model.version,".stan"),
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
    
    
    
    
    # # Figures Section ##
    # 
    # 
    # 
    # 
    # # Extract parameter estimates
    # pars <- rstan::extract(fit)
    # names(pars)
    # 
    # # Extract estimates for prior, likelihood, and posterior 
    # df.run <- data.frame("par" = "Runsize", "value" = pars$RunSize)
    # df.prior <- data.frame("par" = "Prior", "value" = pars$prior_pf)
    # df.pssPred <- data.frame("par" = "PSS_Pred", "value" = pars$curr_predPSS)
    # 
    # # Bind into one dataframe with 2 columns
    # df.comb <- rbind(df.run,df.prior, df.pssPred)
    # 
    # # Density plots comparing pf, linear prediction, and posterior estimate
    # ggplot(df.comb, aes(x = value, fill = par))+
    #   geom_density(alpha = .5)+
    #   ggtitle(paste("Density plot for",myYear,"and day",myDay))+
    #   theme_classic()+
    #   geom_vline(xintercept = CAN_hist$CanCount[CAN_hist$Year == myYear],
    #              col = "Red",
    #              linetype = 2)+
    #   xlab("Number of Chinook Salmon")+
    #   xlim(c(0,250000))
    # 
    # mcmc_intervals(fit, pars = c("alpha","RunSize", "curr_predPSS"), prob = 0.95)
    # 
    
      
      
  }}

   output.version1$`Fit Year = 2016 Day = 212`
pars <- rstan::extract(output.version1$`Fit Year = 2017 Day = 172`)   
