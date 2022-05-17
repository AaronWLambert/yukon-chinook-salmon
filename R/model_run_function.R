#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.14.22
# Version: This function will run any version of the Yukon Inseason Forecast model
# Purpose: Functions that will run the model and return posterior parameter estimates
# 
#
#   1) Read in data
#   2) Use linear model to determine relationship between Canadian passage and
#     PSS and runtiming
#   3) Stan model that estimates total run size
#
#
#=================================================================================
#NOTES:
# This function will be called from another script.
#  This function returns extracted pars, fit, cumPSS (for plotting), and totalEOS 
#   (for plotting)
#
# 
# Next steps: 
# More functions? 
#
#=================================================================================


InSeasonProjection <- function(model.version,
                               myYear,
                               myDay , 
                               n.chains,
                               n.iter, 
                               n.thin ,
                               CAN_hist,
                               PSS_hist,
                               GSI_by_year,
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
  PF_vect <- pf_hist$Weighted_Forecast[pf_hist$Year<=2021]
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
  curr_PSS <- PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay]
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
  count_vect <- PSS_hist$count[PSS_hist$Year<= (myYear-1) &
                                  PSS_hist$Day <= (myDay)]
  
  # The number of observations in count_vector for PSS historic counts
  n_hist_counts <- length(count_vect)
  
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
  
  # Vector containing avg GSI proportions for each day across all years
  meanGSI_vect <- vector(length = length(152:myDay))
  names(meanGSI_vect) <- c(152:myDay)
  sdGSI_vect <- vector(length = length(152:myDay))
  names(sdGSI_vect) <- c(152:myDay)
  
  counter <- 1
  for (d in 152:myDay) {
    meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
    sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
    counter <- counter+1
  }
  
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
  
  adj_curr_PSS <- vector(length = n_dayPSS)
  for (d in 1:n_dayPSS) {
    
    
    adj_curr_PSS[d] <- curr_PSS[d]*meanGSI_vect[d]
  }
  
  
  # Stan Model Call ######################################
  
  
  fit <- stan(file = file.path(dir.stan,paste("Yukon Inseason Forecast ", model.version ,".stan", sep = "")),
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
                         "Pf_sigma" = pf_sigma,
                         "meanpropCAN" = meanGSI_vect,
                         "sd_meanpropCAN" = sdGSI_vect,
                         "mean_adj_PSS_mat"=PSS_mat_adj,
                         "mean_adj_curr_PSS"= adj_curr_PSS,
                         "GSI_mean"= GSI_avg_vect,
                         "GSI_sd" = GSI_sd_vect),
             chains = n.chains,
             iter = n.iter, 
             thin = n.thin, 
             cores = n.chains,
             control = list(max_treedepth = 25,
                            adapt_delta = 0.99)
  )
  
  
  # Extract parameter estimates
  pars <- rstan::extract(fit)
  
  # Extract fit summary for saving
  fit.summary <- summary(fit)
  
  outs <- list("cumPSS" = cumPSS, 
               "cumPSS_adj"= cumPSS_adj,
               "pars" = pars, 
               "totalEOS" = totalEOS, 
               # "fit" = fit, 
               "summary" = fit.summary,
               "myYear" = myYear,
               "version" = model.version,
               "myDay" = myDay)
}

# Plot function for density plots and PSS prediction ribbon plot ########################################################
#'@param outputList A list returned from Inseason function for running the inseason forecast model
outPlots <- function(outputList, CAN_hist, GSI = FALSE){
  outputPlots <- list()


    # i =1
    # Use loop output list from above to call plots and save them to a list ############
    
    # Extract relevant list objects for plotting
    pars <- outputList$pars
    fit <- outputList$fit
    cumPSS <- outputList$cumPSS
    cumPSS_adj <- outputList$cumPSS_adj
    totalEOS <- outputList$totalEOS
    trueCan <- as.double(CAN_hist$can.mean[CAN_hist$Year == outputList$myYear]/1000)
    myD <- outputList$myDay
    myY <- outputList$myYear
    Ver <- outputList$version
    
    # Get the quantiles for predPSS from model output
    quant.predPSS <- apply(X = pars$predPSS, 
                           MARGIN = 2, 
                           FUN = quantile, 
                           probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    
    # If statement to plot with GSI adjusted if TRUE
    if(GSI == FALSE){
    
    pred.plot.df <- data.frame(cumPSS,
                               totalEOS,
                               t(quant.predPSS))
    names(pred.plot.df) <- c("cumPSS", "totalEOS","low95","low50","median",
                             "up50","up95")
    
    
    
    zz <- ggplot(pred.plot.df, aes(x = cumPSS/1000, y = totalEOS/1000 ))+
      geom_point() +
      geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
      geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
      geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
      labs(x = "Cummulative PSS Chinook Passage (1000's)",
           y = "Total Reconstructed EOS Chinook Passage (1000's)")+
      theme(text = element_text(size = 20, family = "serif"),
            panel.grid.major =element_line(),    #strip major gridlines
            panel.grid.minor = element_blank(),    #strip minor gridlines
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"))+
      scale_fill_colorblind(name = "")+
      scale_color_colorblind(name = "")
    } else{
      pred.plot.df <- data.frame(cumPSS_adj,
                                 totalEOS,
                                 t(quant.predPSS))
      names(pred.plot.df) <- c("cumPSS_adj", "totalEOS","low95","low50","median",
                               "up50","up95")
      
      
      
      zz <- ggplot(pred.plot.df, aes(x = cumPSS_adj/1000, y = totalEOS/1000 ))+
        geom_point() +
        geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
        geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
        geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
        labs(x = "Cummulative Canadian-adj PSS Chinook Passage (1000's)",
             y = "Total Reconstructed EOS Chinook Passage (1000's)")+
        theme(text = element_text(size = 20, family = "serif"),
              panel.grid.major =element_line(),    #strip major gridlines
              panel.grid.minor = element_blank(),    #strip minor gridlines
              axis.ticks = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"))+
        scale_fill_colorblind(name = "")+
        scale_color_colorblind(name = "")
      
    }
    
    
    # Extract estimates for prior, likelihood, and posterior 
    df.run <- data.frame("par" = "Runsize", "value" = pars$RunSize)
    df.prior <- data.frame("par" = "Prior", "value" = pars$prior_pf)
    df.pssPred <- data.frame("par" = "PSS_Pred", "value" = pars$curr_predPSS)
    df.postPredPss <- data.frame("par" = "PSS_post_pred", "value" = pars$post_curr_predPSS)
    
    # Bind into one dataframe with 2 columns
    df.comb <- rbind(df.prior,df.postPredPss,df.run)
    # df.comb$par <- levels(df.comb$par, levels = c("Prior","PSS_post_pred","Runsize"))
    # Density plots comparing pf, linear prediction, and posterior estimate
    postDense <-ggplot(df.comb, aes(x = value/1000, fill = par))+
      geom_density(alpha = .65)+
      geom_vline( aes(xintercept = trueCan,
                      col = "True Runsize"),
                  linetype = 2,
                  size = 1.5)+
      # ylab("Relative Probability")+
      ylab("")+
      xlab("")+
      # ggtitle(paste("Version = ",Ver,"\n Day = ",myD,"\n Year = ",myY, sep = ""))+
      # xlab("1000's of Chinook Salmon")+
      coord_cartesian(xlim = c(0,200))+
      scale_fill_colorblind(name = "", 
                            labels = c( "Preseason Forecast (Prior)", 
                                        "PSS Prediction","Runsize Projection"))+
      theme_tidybayes()+
      theme(text = element_text(size = 20, family = "serif"),plot.title = element_text(size = 12))+
      scale_color_discrete(name = "")
    outputPlots <- list("DensPlot" = postDense, "PredPlot"= zz)
  
  return(outputPlots)
}


