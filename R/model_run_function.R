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
#     PSS and run timing
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


# InSeasonProjection <- function(model.version,
#                                myYear,
#                                myDay , 
#                                n.chains,
#                                n.iter, 
#                                n.thin ,
#                                CAN_hist,
#                                PSS_hist,
#                                GSI_by_year,
#                                pf_hist,
#                                Eagle_hist,
#                                # ps_mu,
#                                # ps_sd,
#                                # ps_alpha_norm,
#                                # ps_m,
#                                # ps_s,
#                                # ps_alpha_log,
#                                savefit = FALSE,
#                                logistic = FALSE,
#                                normal = FALSE,
#                                multiplier = 1,
#                                startDayPSS = 148,
#                                startYearPSS = 1995,
#                                startYearEagle = 2005){
#   
#   
#   
#   
#   
# 
#   
#   # Start Years for Predictors
#   startYearPF <- 2007
#   startYearGSI <- 2005
#   # startYearPSS <-1997
#   
#   # Wont typicaly change;
#   # day 152 = June 1
#   # startDayPSS <-148
#   endDayPSS <- 225
#   # Preseason Forcast######################
#   
#   # Current Preseason forecast Can Origin
#   pf <- log(pf_hist$Weighted_Forecast[pf_hist$Year == myYear])
#   # pf_sigma <- (0.6)
#   # PF_vect <- pf_hist$Weighted_Forecast[pf_hist$Year<=2021]#*&*& Change back please
#   PF_vect <- pf_hist$Weighted_Forecast[pf_hist$Year<=2019]
#   EOS <- CAN_hist$can.mean[CAN_hist$Year >= 2007 &
#                              CAN_hist$Year != 2011 &
#                              CAN_hist$Year != 2012]
#   
#   pf_sigma <- sd(log(EOS/PF_vect))
#   
#   
#   
#   # Metadata###############################
#   
#   # PSS years used
#   yearPSS <- unique(PSS_hist$Year[PSS_hist$Year<=(myYear-1)&
#                       PSS_hist$Year >= startYearPSS]) #start year should equal 1995 to include all years data
#   n_yearPSS <- length(yearPSS)
#   
#   # PSS Days 
#   # dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])
#   # n_dayPSS <- length(dayPSS)
#   # dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay) & PSS_hist$Day>= startDayPSS])-147
#   dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay) & PSS_hist$Day>= startDayPSS])
#   # Number of days used
#   n_dayPSS <- length(dayPSS)
#   
#   # Vector of historic PSS years up to myYear-1
#   yearEagle <- unique(Eagle_hist$Year[Eagle_hist$Year<=(myYear-1)])
#   
#   # Number of years used in model
#   n_yearEagle <- length(yearEagle)
#   
#   # # PSS days included up to myDay. 151 is subtracted for use in logistic model
#   # #  I.e., June 1 = 1, June 2 = 2 ....
#   # dayPSS <- unique(PSS_hist$Day[PSS_hist$Day <=(myDay)])-147 *******
#   dayEagle <- unique(Eagle_hist$Day[Eagle_hist$Day <=(myDay) & Eagle_hist$Day>= startDayPSS])-147
#   
#   # Number of days used
#   n_dayEagle <- length(dayEagle)
#   
#   # Eagle daily passage estimate for days up to myDay for the year of interest (myYear)
#   # curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay ]) **********
#   curr_Eagle_vect <- (Eagle_hist$count[Eagle_hist$Year == myYear & Eagle_hist$Day <= myDay & Eagle_hist$Day>=startDayPSS])
#   
#   
#   curr_Eagle <- array(data = curr_Eagle_vect,
#                       dim = length(curr_Eagle_vect))
#   # Number of days used
#   # Need this so stan will accept a 0 value for curr_pss
#   # if(length(curr_Eagle)==1 & curr_Eagle[1]<1){
#   #   curr_Eagle <- numeric(0)
#   # }
#   
#   n_curr_Eagle <- length(curr_Eagle)
#   # TAKE ME OUT##
#   cum_curr_Eagle <- cumsum(curr_Eagle)
#   # PSS Counts for myYear
#   # curr_PSS <- PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay]
#   # n_curr_PSS <- length(curr_PSS)
#   curr_PSS <- (PSS_hist$count[PSS_hist$Year == myYear & PSS_hist$Day <= myDay & PSS_hist$Day>=startDayPSS])
#   # Number of days used
#   n_curr_PSS <- length(curr_PSS)
#   # Historic preseason forcast years
#   # yearPF <- pf_hist$Year[pf_hist$Year <= (myYear-1)]
#   # n_yearPF <- length(yearPF)
#   
#   # # Historical  Canadian Origin PF (2013 - current)
#   # histPF <- pf_hist$Weighted_Forecast[pf_hist$Year <= (myYear-1)]
#   # names(histPF) <- (pf_hist$Year[pf_hist$Year <= (myYear-1)]) 
#   # n_histPF <-length(histPF)
#   
#   # End of season Canadian counts 1995-2020 excluding 1996 for comparison to PSS
#   totalEOS <- CAN_hist$can.mean[CAN_hist$Year <= (myYear-1) &
#                                   CAN_hist$Year >= startYearPSS &
#                                   CAN_hist$Year != 1996]
#   
#   # names(totalEOS) <- CAN_hist$Year[CAN_hist$Year <= (myYear-1) &
#   #                                   CAN_hist$Year >= startYearPSS &
#   #                                   CAN_hist$Year != 1996]
#   
#   n_totalEOS <- length(totalEOS)
#   
#   
#   
#   
#   
#   # Eagle Matrix
#   # Create matrix with dimensions myday and myYear
#   Eagle_mat <- matrix(nrow = (length(startDayPSS:myDay )),
#                       #start 1996 to account for missing year 1996
#                       ncol = length(2005:(myYear-1))) 
#   
#   # Give names to matrix for accounting
#   colnames(Eagle_mat) <- c(2005:(myYear-1))
#   
#   rownames(Eagle_mat) <- c(startDayPSS:(myDay))
#   
#   # Create vector of counts from relevant years and days to input into matrix
#   # (count_vect <- PSS_hist$count[PSS_hist$Year<= (myYear-1) & *******************************************
#   #                                 PSS_hist$Day <= (myDay)])
#   (count_vect_Eagle <- Eagle_hist$count[Eagle_hist$Year<= (myYear-1) &
#                                           Eagle_hist$Day <= (myDay) &
#                                           Eagle_hist$Day >= (startDayPSS)])
#   
#   # The number of observations in count_vector for PSS historic counts
#   n_hist_Eaglecounts <- length(count_vect_Eagle)
#   
#   # Use loop to populate matrix with counts by days
#   # Set counter
#   if(myDay >= 178){
#   counter <-1
#   
#   for (y in 1:length(2005:(myYear-1))){
#     for (d in 31:(length(startDayPSS:myDay))) {
#       
#       
#       Eagle_mat[d,y] <- count_vect_Eagle[counter]
#       counter = counter+1
#       
#     }
#   }
#   } # end of if statement
#   
#   Eagle_mat[is.na(Eagle_mat)]<-0
#   
#   # Cumulative count matrix for historic years
#   Eagle_cum_hist_mat <- apply(X = Eagle_mat,MARGIN = 2,FUN = cumsum)
#   
#   cumEagle <- apply(Eagle_mat,MARGIN = 2,sum)
#   
#   # Create matrix with PSS daily counts from day 152 to my day and 1995
#   # to my year
#   #  *Note that 1996 is missing*
#   
#   
#   # Create matrix with dimensions myday and myYear
#   PSS_mat <- matrix(nrow = (length(startDayPSS:myDay )),
#                     #start 1996 to account for missing year 1996
#                     ncol = length(1996:(myYear-1))) 
#   
#   # Give names to matrix
#   # colnames(PSS_mat) <- c("1995", startYearPSS:(myYear-1))
#   # rownames(PSS_mat) <- c(startDayPSS:(myDay))
#   
#   # Make sure it works...
#   # dim(PSS_mat)
#   # str(PSS_hist)
#   
#   # Create vector of counts from relevant years and days
#   # count_vect <- PSS_hist$count[PSS_hist$Year<= (myYear-1) &
#   #                                 PSS_hist$Day <= (myDay)]
#   (count_vect <- PSS_hist$count[PSS_hist$Year <= (myYear-1) &
#                                   PSS_hist$Year >= startYearPSS &
#                                   PSS_hist$Day <= (myDay) &
#                                   PSS_hist$Day >= (startDayPSS)])
#   # The number of observations in count_vector for PSS historic counts
#   n_hist_counts <- length(count_vect)
#   
#   # Use loop to populate matrix with counts by days
#   
#   # Set counter
#   counter <-1
#   y <-0
#   for (y in 1:length(1996:(myYear-1))){
#     for (d in 1:(length(startDayPSS:myDay))) {
#       
#       
#       PSS_mat[d,y] <- count_vect[counter]
#       counter = counter+1
#       
#     }
#   }
#   
#   #### Double first entry that is non-zero to make up for non-true zero days #############
#   # Set multiplier to 1 in function to use reported values for the first observation in each
#   # year.
#   # for(y in 1:n_yearPSS){
#   #   
#   #   PSS_mat[,y][PSS_mat[,y]>0][1] <- PSS_mat[,y][PSS_mat[,y]>0][1]*multiplier
#   # }
#   # PSS matrix of complete historic counts ***********
#   # Create matrix with dimensions myday and myYear
#   PSS_mat_all <- matrix(nrow = (length(startDayPSS:endDayPSS )),
#                         #start 1996 to account for missing year 1996
#                         ncol = length(1996:(myYear-1))) 
#   
#   # Give names to matrix for accounting
#   # colnames(PSS_mat_all) <- c("1995", startYearPSS:(myYear-1))
#   # rownames(PSS_mat_all) <- c(startDayPSS:endDayPSS)
#   
#   # Create vector of counts from relevant years and days to input into matrix
#   (count_vect_all <- PSS_hist$count[PSS_hist$Year<= (myYear-1) &
#                                       PSS_hist$Year >= startYearPSS &
#                                       PSS_hist$Day <= endDayPSS])
#   
#   # The number of observations in count_vector for PSS historic counts
#   (n_hist_counts_all <- length(count_vect_all))
#   
#   # Use loop to populate matrix with counts by days
#   # Set counter
#   counter <-1
#   
#   for (y in 1:length(1996:(myYear-1))){
#     for (d in 1:(length(startDayPSS:endDayPSS))) {
#       
#       
#       PSS_mat_all[d,y] <- count_vect_all[counter]
#       counter = counter+1
#       
#     }
#   }
#   # Vector for storing cumPSS for plotting
#   cumPSS <- vector(length = n_yearPSS)
#   
#   names(cumPSS)<- c(yearPSS)
#   # Cummulative PSS Counts for plotting
#   for (i in 1:n_yearPSS) {
#     
#     cumPSS[i]<- sum(PSS_mat[,i])
#   }
#   
#   # Vector containing avg GSI proportions for each day across all years
#   meanGSI_vect <- vector(length = length(startDayPSS:myDay))
#   names(meanGSI_vect) <- c(startDayPSS:myDay)
#   sdGSI_vect <- vector(length = length(startDayPSS:myDay))
#   names(sdGSI_vect) <- c(startDayPSS:myDay)
#   # meanGSI_vect <- vector(length = length(startDayPSS:endDayPSS))
#   # names(meanGSI_vect) <- c(startDayPSS:endDayPSS)
#   # sdGSI_vect <- vector(length = length(startDayPSS:endDayPSS))
#   # names(sdGSI_vect) <- c(startDayPSS:endDayPSS)
#   counter <- 1
#   for (d in startDayPSS:myDay) {
#     meanGSI_vect[counter]<- mean(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
#     sdGSI_vect[counter]<- sd(GSI_by_year$propCan[GSI_by_year$startday <= d & GSI_by_year$endday >=d])
#     counter <- counter+1
#   }
#   # No GSI data for day 148 and 149. Use day 150 for these values
#   meanGSI_vect[1:2] <- meanGSI_vect[3]
#   sdGSI_vect[1:2] <- sdGSI_vect[3]
#   
#   meanGSI_vect[is.na(meanGSI_vect)] <- meanGSI_vect["248"]
#   sdGSI_vect[is.na(sdGSI_vect)] <- sdGSI_vect["248"]
#   # Mean GSI by mean strata dates. This is used in versions 3.4 and up #############################
#   meanStartDay <-GSI_by_year %>% group_by(stratum) %>% 
#     summarise("meanStartDay" = mean(startday)) %>% 
#     as.data.frame()
#   
#   GSI_mean_by_strata <-GSI_by_year %>% 
#     summarize("stratumMean" = c(mean(propCan[stratum == 1]),
#                                 mean(propCan[stratum == 2]),
#                                 mean(propCan[stratum == 3 | stratum == 4])),
#               "stratumSD" = c(sd(propCan[stratum == 1]),
#                               sd(propCan[stratum == 2]),
#                               sd(propCan[stratum == 3 | stratum == 4]))) %>%  as.data.frame()
#   
#   GSI <- cbind(GSI_mean_by_strata,round(meanStartDay[1:3,]))
#   
#   GSI_avg<-c(rep(GSI$stratumMean[GSI$stratum == 1], 
#                  times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
#              rep(GSI$stratumMean[GSI$stratum == 2], 
#                  times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
#              rep(GSI$stratumMean[GSI$stratum == 3],
#                  times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))
#   GSI_avg_vect <- GSI_avg[1:length(startDayPSS:myDay)]
#   
#   GSI_sd <- c(rep(GSI$stratumSD[GSI$stratum == 1], 
#                   times = length(startDayPSS:GSI$meanStartDay[GSI$stratum==2]-1)),
#               rep(GSI$stratumSD[GSI$stratum == 2], 
#                   times = length(GSI$meanStartDay[GSI$stratum==2]:GSI$meanStartDay[GSI$stratum == 3]-1)),
#               rep(GSI$stratumSD[GSI$stratum == 3],
#                   times = length(GSI$meanStartDay[GSI$stratum == 3]:max(PSS_hist$Day))))
#   GSI_sd_vect <- GSI_sd[1:length(startDayPSS:myDay)]
#   # Create matrix with dimensions myday and myYear
#   PSS_mat_adj <- matrix(nrow = (length(startDayPSS:myDay )),
#                         #start 1996 to account for missing year 1996
#                         ncol = length(startYearPSS:(myYear-1))) 
#   
#   # Give names to matrix
#   # colnames(PSS_mat_adj) <- c("1995", startYearPSS:(myYear-1))
#   # rownames(PSS_mat_adj) <- c(startDayPSS:(myDay))
#   
#   # Make sure it works...
#   # dim(PSS_mat_adj)
#   
#   
#   for (y in 1:n_yearPSS) {
#     for (d in 1:n_dayPSS) {
#       # y = 1
#       # d = 15
#       PSS_mat_adj[d,y] <- PSS_mat[d,y]*meanGSI_vect[d] 
#       
#     }
#     
#   }
#   
#   # Vector for storing cumPSS for plotting
#   cumPSS_adj <- vector(length = n_yearPSS)
#   
#   names(cumPSS_adj)<- c(yearPSS)
#   # Cummulative PSS Counts for plotting
#   for (i in 1:n_yearPSS) {
#     
#     cumPSS_adj[i]<- sum(PSS_mat_adj[,i])
#   }
#   
#   adj_curr_PSS <- vector(length = n_dayPSS)
#   for (d in 1:n_dayPSS) {
#     
#     
#     adj_curr_PSS[d] <- curr_PSS[d]*meanGSI_vect[d]
#   }
#   
#   # Cumulative PSS matrix
#   cum_PSS_mat <-  apply(PSS_mat, 2, cumsum)
#   
#   # Cumulative current PSS counts
#   cum_curr_PSS <- cumsum(curr_PSS)
#   
#   if(logistic == TRUE){
#     # # counts by GSI proportions for fitting curve in function
#     gsiDF <- data.frame("GSI" = meanGSI_vect, "Day" = c(148:endDayPSS))
#     
#     gsiDF$GSI[is.na(gsiDF$GSI)] <- gsiDF$GSI[gsiDF$Day == 248]
#     
#     PSS_hist_adj <- left_join(x = PSS_hist, y = gsiDF)
#     
#     PSS_hist_adj$adj_count <- PSS_hist_adj$count*PSS_hist_adj$GSI
#     
#     temp.df <- PSS_hist_adj %>% group_by(Year) %>% summarize(count = cumsum(count))
#     
#     PSS_hist_adj <- cbind(PSS_hist_adj,temp.df)
#     
#     norm.prior <- TRUE
#     if(norm.prior==TRUE) {
#       years <- c(1995,1997:(myYear-1))
#       n.years <- length(years)
#       
#       norm.m <- vector(length=n.years)
#       norm.s <- vector(length=n.years)
#       norm.alpha <- vector(length=n.years)
#       year.df <- vector(length = n.years)
#       # pdf(file = file.path(dir.output,"MLE Logistic Curve fit 20Jul22"), height=8, width=8)
#       
#       # y <- 24
#       for(y in 1:n.years) {
#         year <- years[y]
#         print(year)
#         # Extract data
#         out <- extract.data(df = PSS_hist_adj, year=year, GSI = FALSE)
#         #Fit normal data
#         fit.norm <-  mle2(like.norm.logistic,
#                           start=list(m=log(170), s=log(5.5), alpha=log(60000)),
#                           fixed=list(plot=FALSE, sigma=log(20)),
#                           data=list(days=out$days, rets=out$rets, dist='ssq'),
#                           method='Nelder-Mead', control=list(maxit=100000, trace=FALSE))
#         norm.s[y] <- exp(coef(fit.norm)[1])
#         norm.m[y] <- exp(coef(fit.norm)[2])
#         norm.alpha[y] <- exp(coef(fit.norm)[3])
#         year.df[y] <- year
#         like.norm.logistic(days=out$days,
#                   rets=out$rets,
#                   (coef(fit.norm)[1]),
#                   coef(fit.norm)[2],
#                   coef(fit.norm)[3],
#                   coef(fit.norm)[4],
#                   plot= FALSE,
#                   dist='ssq')
#         # mtext(year, side=3, line=2, font=2)
#       }#next y
#       
#       # dev.off()
#       prior.df.log <- data.frame(norm.m,norm.s,norm.alpha, year.df)
#       names(prior.df.log) <- c('m','s','alpha',"year")
#       # write.csv(prior.df,file = file.path(dir.output,"logisitic curve prior.csv"))
#     }
#     ps_m <- prior.df.log$m
#     ps_s <- prior.df.log$s
#     ps_alpha.log <- prior.df.log$alpha
#   }# End if statement for logistic fit function
#   
#   if(normal==TRUE){
#     
#     gsiDF <- data.frame("GSI" = meanGSI_vect, "Day" = c(startDayPSS:endDayPSS))
#     
#     gsiDF$GSI[is.na(gsiDF$GSI)] <- gsiDF$GSI[gsiDF$Day == endDayPSS]
#     
#     PSS_hist_adj <- left_join(x = PSS_hist, y = gsiDF)
#     
#     PSS_hist_adj$adj_count <- PSS_hist_adj$count*PSS_hist_adj$GSI
#     
#     norm.prior <-TRUE
#     if(norm.prior==TRUE) {
#       years <- c(1995,1997:(myYear-1))
#       n.years <- length(years)
#       
#       norm.mu <- vector(length=n.years)
#       norm.sd <- vector(length=n.years)
#       norm.alpha <- vector(length=n.years)
#       # pdf(file = file.path(dir.output,"4.0 Normal MLE Fit 15Jul22"), height=8, width=8)
#       # par(mfrow = c(3,3))
#       # y <- 1
#       for(y in 1:n.years) {
#         year <- years[y]
#         print(year)
#         # Extract data NOTE! Change GSI if needed for model
#         out <- extract.data(df = PSS_hist_adj, year=year, GSI = F)
#         #Fit normal data
#         fit.norm <-  mle2(like.norm,
#                           start=list(mu=log(170), sd=log(20), alpha=log(4000)),
#                           fixed=list(plot=FALSE, sigma=log(20)),
#                           data=list(days=out$days, rets=out$rets, dist='ssq'),
#                           method='Nelder-Mead', control=list(maxit=100000, trace=FALSE))
#         norm.mu[y] <- exp(coef(fit.norm)[1])
#         norm.sd[y] <- exp(coef(fit.norm)[2])
#         norm.alpha[y] <- exp(coef(fit.norm)[3])
#         like.norm(days=out$days,
#                   rets=out$rets,
#                   (coef(fit.norm)[1]),
#                   coef(fit.norm)[2],
#                   coef(fit.norm)[3],
#                   coef(fit.norm)[4],
#                   plot=FALSE,
#                   dist='ssq')
#         # mtext(year, side=3, line=2, font=2)
#       }#next y
#       
#       # dev.off()
#       prior.df <- data.frame(norm.mu,norm.sd,norm.alpha)
#       names(prior.df) <- c('mu','sd','alpha')
#       # write.csv(x = prior.df, file = file.path(dir.output,"normal curve prior GSI.csv"))
#     }
#     ps_mu <- prior.df$mu
#     ps_sd <- prior.df$sd
#     ps_alpha_norm <- prior.df$alpha
#     
#   }
#   
#   pss_days_all <- startDayPSS:endDayPSS
#   n_pss_days_all <- length(pss_days_all)
#   
# 
# 
#   # Stan Model Call ######################################
#   if(logistic==TRUE){
#   
#   fit <- stan(file = file.path(dir.stan,paste("Yukon Inseason Forecast ", model.version ,".stan", sep = "")),
#               data = list("PSS_mat"=PSS_mat,
#                          # "n_histPF" =n_histPF, 
#                          # "histPF"=histPF,
#                          "n_totalEOS"=n_totalEOS,
#                          "totalEOS"=totalEOS, 
#                          "n_dayPSS"=n_dayPSS,
#                          "dayPSS"=dayPSS, 
#                          "n_yearPSS"=n_yearPSS,
#                          "yearPSS"=yearPSS,
#                          "n_curr_PSS"=n_curr_PSS,
#                          "curr_PSS"=curr_PSS,
#                          "n_hist_counts"=n_hist_counts,
#                          "Pf" = pf, 
#                          "Pf_sigma" = pf_sigma,
#                          "meanpropCAN" = meanGSI_vect,
#                          "sd_meanpropCAN" = sdGSI_vect,
#                          "mean_adj_PSS_mat"=PSS_mat_adj,
#                          "mean_adj_curr_PSS"= adj_curr_PSS,
#                          "GSI_mean"= GSI_avg_vect,
#                          "GSI_sd" = GSI_sd_vect,
#                          # "ps_mu"= ps_mu,
#                          # "ps_sd"= ps_sd,
#                          # "ps_alpha"= ps_alpha_norm,
#                          "ps_m"= ps_m,
#                          "ps_s"= ps_s,
#                          "ps_alpha_log" = ps_alpha.log,
#                          "pss_days_all" = pss_days_all,
#                          "n_pss_days_all" = n_pss_days_all,
#                          # "PSS_mat_adj" = PSS_mat_adj,
#                          "myDay" = myDay,
#                          "PSS_mat_all" = PSS_mat_all,
#                          "cum_PSS_mat" = cum_PSS_mat,
#                          "cum_curr_PSS" = cum_curr_PSS,
#                          "Eagle_mat" = Eagle_mat, 
#                          "n_dayEagle" = n_dayEagle,
#                          "dayEagle" = dayEagle, 
#                          "n_yearEagle" = n_yearEagle,
#                          "yearEagle" = yearEagle,
#                          "n_curr_Eagle" = n_curr_Eagle,
#                          "curr_Eagle" = curr_Eagle),
#              # init = init_ll,
#              chains = n.chains,
#              iter = n.iter, 
#              thin = n.thin, 
#              cores = n.chains,
#              control = list(max_treedepth = 25,
#                             adapt_delta = 0.99)
#   )
#   
#   }else{if(normal == TRUE){
#     
#     fit <- stan(file = file.path(dir.stan,paste("Yukon Inseason Forecast ", model.version ,".stan", sep = "")),
#                 data = list("PSS_mat"=PSS_mat,
#                             # "n_histPF" =n_histPF, 
#                             # "histPF"=histPF,
#                             "n_totalEOS"=n_totalEOS,
#                             "totalEOS"=totalEOS, 
#                             "n_dayPSS"=n_dayPSS,
#                             "dayPSS"=dayPSS, 
#                             "n_yearPSS"=n_yearPSS,
#                             "yearPSS"=yearPSS,
#                             "n_curr_PSS"=n_curr_PSS,
#                             "curr_PSS"=curr_PSS,
#                             "n_hist_counts"=n_hist_counts,
#                             "Pf" = pf, 
#                             "Pf_sigma" = pf_sigma,
#                             "meanpropCAN" = meanGSI_vect,
#                             "sd_meanpropCAN" = sdGSI_vect,
#                             "mean_adj_PSS_mat"=PSS_mat_adj,
#                             "mean_adj_curr_PSS"= adj_curr_PSS,
#                             "GSI_mean"= GSI_avg_vect,
#                             "GSI_sd" = GSI_sd_vect,
#                             "ps_mu"= ps_mu,
#                             "ps_sd"= ps_sd,
#                             "ps_alpha"= ps_alpha_norm,
#                             # "ps_m"= ps_m,
#                             # "ps_s"= ps_s,
#                             # "ps_alpha_log" = ps_alpha_log,
#                             "pss_days_all" = pss_days_all,
#                             "n_pss_days_all" = n_pss_days_all,
#                             # "PSS_mat_adj" = PSS_mat_adj,
#                             "myDay" = myDay,
#                             "PSS_mat_all" = PSS_mat_all,
#                             "cum_PSS_mat" = cum_PSS_mat,
#                             "cum_curr_PSS" = cum_curr_PSS,
#                             "Eagle_mat" = Eagle_mat, 
#                             "n_dayEagle" = n_dayEagle,
#                             "dayEagle" = dayEagle, 
#                             "n_yearEagle" = n_yearEagle,
#                             "yearEagle" = yearEagle,
#                             "n_curr_Eagle" = n_curr_Eagle,
#                             "curr_Eagle" = curr_Eagle),
#                 # init = init_ll,
#                 chains = n.chains,
#                 iter = n.iter, 
#                 thin = n.thin, 
#                 cores = n.chains,
#                 control = list(max_treedepth = 25,
#                                adapt_delta = 0.99)
#     )
#     
#   }
#     
#   }
#   if(logistic == FALSE & normal ==FALSE){
#     
#     
#     # init_fn <-  function(){
#     #   list("ps_alpha_curr" = runif(1,200,1000),
#     #        "ps_mu_curr" = runif(1,172,181),
#     #        "ps_sd_curr" = runif(1,10,20),
#     #        "ps_alpha_hist" = runif(n_yearPSS,200,1000),
#     #        "ps_mu_hist" = runif(n_yearPSS,172,181),
#     #        "ps_sd_hist" = runif(n_yearPSS,10,20))
#     # }
#     # 
#     # 
#     # 
#     # # 
#     # init_ll <- list(init_fn(), init_fn(), init_fn(),init_fn())
#     
#     fit <- stan(file = file.path(dir.stan,paste("Yukon Inseason Forecast ", model.version ,".stan", sep = "")),
#                 data = list("PSS_mat"=PSS_mat,
#                             # "n_histPF" =n_histPF, 
#                             # "histPF"=histPF,
#                             "n_totalEOS"=n_totalEOS,
#                             "totalEOS"=totalEOS, 
#                             "n_dayPSS"=n_dayPSS,
#                             "dayPSS"=dayPSS, 
#                             "n_yearPSS"=n_yearPSS,
#                             "yearPSS"=yearPSS,
#                             "n_curr_PSS"=n_curr_PSS,
#                             "curr_PSS"=curr_PSS,
#                             "n_hist_counts"=n_hist_counts,
#                             "Pf" = pf, 
#                             "Pf_sigma" = pf_sigma,
#                             "meanpropCAN" = meanGSI_vect,
#                             "sd_meanpropCAN" = sdGSI_vect,
#                             "mean_adj_PSS_mat"=PSS_mat_adj,
#                             "mean_adj_curr_PSS"= adj_curr_PSS,
#                             "GSI_mean"= GSI_avg_vect,
#                             "GSI_sd" = GSI_sd_vect,
#                             # "ps_mu"= ps_mu,
#                             # "ps_sd"= ps_sd,
#                             # "ps_alpha"= ps_alpha_norm,
#                             # "ps_m"= ps_m,
#                             # "ps_s"= ps_s,
#                             # "ps_alpha_log" = ps_alpha_log,
#                             "pss_days_all" = pss_days_all,
#                             "n_pss_days_all" = n_pss_days_all,
#                             # "PSS_mat_adj" = PSS_mat_adj,
#                             "myDay" = myDay,
#                             "PSS_mat_all" = PSS_mat_all,
#                             "cum_PSS_mat" = cum_PSS_mat,
#                             "cum_curr_PSS" = cum_curr_PSS,
#                             "Eagle_mat" = Eagle_mat,
#                             "n_dayEagle" = n_dayEagle,
#                             "dayEagle" = dayEagle, 
#                             "n_yearEagle" = n_yearEagle,
#                             "yearEagle" = yearEagle,
#                             "n_curr_Eagle" = n_curr_Eagle,
#                             "curr_Eagle" = curr_Eagle),
#                 # init = init_ll,
#                 chains = n.chains,
#                 iter = n.iter, 
#                 thin = n.thin, 
#                 cores = n.chains,
#                 control = list(max_treedepth = 25,
#                                adapt_delta = 0.99)
#     )
#     }
#   # Extract parameter estimates
#   pars <- rstan::extract(fit)
#   
#   # Extract fit summary for saving
#   fit.summary <- summary(fit)
#   if(savefit == FALSE){
#   outs <- list("cumPSS" = cumPSS, 
#                "cumPSS_adj"= cumPSS_adj,
#                "cumEagle" = cumEagle,
#                "pars" = pars, 
#                "totalEOS" = totalEOS,
#                "summary" = fit.summary,
#                "myYear" = myYear,
#                "version" = model.version,
#                "myDay" = myDay)
#   }else{
#     outs <- list("cumPSS" = cumPSS, 
#                  "cumPSS_adj"= cumPSS_adj,
#                  "cumEagle" = cumEagle,
#                  "pars" = pars, 
#                  "totalEOS" = totalEOS,
#                  "summary" = fit.summary,
#                  "myYear" = myYear,
#                  "version" = model.version,
#                  "myDay" = myDay,
#                  "fit" = fit)
#   }#Savefit for false
# }

# Plot function for density plots and PSS prediction ribbon plot ########################################################
#'@param outputList A list returned from Inseason function for running the inseason forecast model




outPlots <- function(outputList, CAN_hist, GSI = FALSE, Retrospective = FALSE, eagle = FALSE){
  outputPlots <- list()


    # i =1
    # Use loop output list from above to call plots and save them to a list ############
    
    # Extract relevant list objects for plotting
    pars <- outputList$pars
    fit <- outputList$fit
    cumPSS <- outputList$cumPSS
    cumPSS_adj <- outputList$cumPSS_adj
    cumEagle <- outputList$cumEagle
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
    
    
    # If statement to plot without retrospective testing
    #  When Retrospective = FALSE, then the realized Canadian-origin line is omitted
    if(Retrospective == TRUE){
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
      coord_cartesian(xlim = c(0,200),
                      ylim = c(0,.0315))+
      scale_fill_colorblind(name = "", 
                            labels = c( "Preseason Forecast (Prior)", 
                                        "PSS Prediction","Runsize Projection"))+
      theme_tidybayes()+
      theme(text = element_text(size = 20, family = "serif"),plot.title = element_text(size = 12))+
      scale_color_discrete(name = "")
    
    }else{ # Else statement for Retrospective argument
      
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
        geom_vline( aes(xintercept = median(pars$RunSize)/1000,
                        col = "Median Projection"),
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
    }
    
    
    
    # Eagle regression plot
    if(eagle == TRUE){
      
      quant.predEagle <- apply(X = pars$predEagle, 
                             MARGIN = 2, 
                             FUN = quantile, 
                             probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
      
      eagle.pred.plot.df <- data.frame(cumEagle,
                                 totalEOS[10:max(length(totalEOS))],
                                 t(quant.predEagle))
      names(eagle.pred.plot.df) <- c("cumEagle", "totalEOS","low95","low50","median",
                               "up50","up95")
      
      
      
      zz.Eagle <- ggplot(eagle.pred.plot.df, aes(x = cumEagle/1000, y = totalEOS/1000 ))+
        geom_point() +
        geom_ribbon(aes(ymin=low95/1000, ymax=up95/1000,fill = "HDI 95"), alpha=0.25) +
        geom_ribbon(aes(ymin=low50/1000, ymax=up50/1000, fill = "HDI 50"), alpha=0.25) +
        geom_line(aes(y=median/1000, color = "Median"), show.legend = T) +
        labs(x = "",
             y = "")+
        theme(text = element_text(size = 20, family = "serif"),
              panel.grid.major =element_line(),    #strip major gridlines
              panel.grid.minor = element_blank(),    #strip minor gridlines
              axis.ticks = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"))+
        scale_fill_colorblind(name = "")+
        scale_color_colorblind(name = "")
      
      
      
      
      # Posterior EOS Can-orig run size estimate
      df.run <- data.frame("par" = "Runsize",
                           "value" = pars$RunSize)
      
      # PF (prior)
      df.prior <- data.frame("par" = "Prior", 
                             "value" = pars$prior_pf)
      
      # PSS prediction density
      df.postPredPss <- data.frame("par" = "PSS_post_pred", 
                                   "value" = pars$post_curr_predPSS)
      
      # Eagle prediction density
      df.postPredEagle <- data.frame("par" = "Eagle_post_pred", 
                                     "value" = pars$post_curr_predEagle)
      
      # Bind into one data frame with 2 columns for plotting
      df.comb <- rbind(df.prior,
                       df.postPredPss,
                       df.postPredEagle,
                       df.run)
      str(df.comb)
      df.comb$par <- as.factor(df.comb$par)
      
      df.comb$par <- relevel(df.comb$par,"Prior")
      
      # Density plots comparing pf, linear prediction, and posterior estimate
      Eagle.Dense<- ggplot(df.comb, aes(x = value/1000, fill = par))+
        geom_density(alpha = .5)+
        # ggtitle(paste("Density plot for",myY,", day",myD,"\n Version", model.version))+
        theme_classic()+
        geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myY]/1000,
                       color = "End of Season Runsize"),
                   linetype = 2,
                   size = 1)+
        
        # xlim(c(0,300000))+
        coord_cartesian(xlim = c(0,200),
                        ylim = c(0,.04))+
        # ylim(c(0,500))
        labs(x = "", 
             fill = "Parameters", 
             y = "")+
        # scale_fill_discrete(name = "Parameters", 
        # labels = c( "Preseason Forecast (Prior)", 
        # "PSS Prediction","Runsize"))+
        scale_x_continuous()+
        scale_fill_colorblind(name = "", 
                              labels = c( "Preseason Forecast (Prior)", 
                                          "Eagle Prediction",
                                          "PSS Prediction",
                                          "Runsize"))+
        # theme_tidybayes(element_text(size = 20))+
        theme(text = element_text(size = 20, family = "serif"),
              panel.grid.major =element_line(),    #strip major gridlines
              panel.grid.minor = element_blank(),    #strip minor gridlines
              axis.ticks = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"))+
        scale_color_discrete(name = "")
      
      
    }
      
      
      
    if(eagle == FALSE){
    outputPlots <- list("DensPlot" = postDense, "PredPlot"= zz)
    }else{
      outputPlots <- list("DensPlot" = Eagle.Dense, "PredPlot"= zz, "EaglePlot" = zz.Eagle)
    }
    
    
  
  return(outputPlots)
}


# Function for easily getting day of year for this project
myDay_func <- function(Month, Day){
  # Function to get doy for the Inseason Bayesian Projection Model
  # Enter Month as a number
  #Testing
  # Month <- 6 ; Day <- 25
  
  # Required package for function
  library(lubridate)
  
  # Df for storing dates
  date.df <- data.frame("dayofyear" = c(148:260), 
                        "date" = seq(from = as.Date("2022-05-28"), 
                                     to = as.Date("2022-9-17"), by=1) )
  # Extract month and day
  date.df$month <-  month(date.df$date)
  date.df$day <- day(date.df$date)
  
  MYDAY <- as.double(date.df$dayofyear[date.df$month == Month & date.df$day == Day])
  
  return(MYDAY)
}

# Extract Data Function ##########################################################################
# Used for curve fitting with mle2

extract.data <- function(df = PSS_hist_adj, year, GSI) {
  ### TESTING ###
  #baywide <- TRUE
  #dist <- 324
  # Year <- 2000
  ###############
  temp.data <- df[df$Year == year,]
  min.day <- min(temp.data$Day)
  max.day <- max(temp.data$Day)
  days <- min.day:max.day
  n.days <- length(days)
  
  #Calculate standardized date values
  global.days <- min(df$Day):max(df$Day)
  std.days <- which(global.days %in% days)
  
  rets <- vector(length=n.days)
  
  # rets <- df$CumCount[df$Year==year]
  if(GSI==TRUE){rets <- df$adj_count[df$Year==year]}else{ rets <- df$count[df$Year==year]}
  #OUTPUT SECTION
  output <- NULL
  # if(baywide==TRUE) { output$dist <- NA }else { output$dist <- dist }
  output$min.day <- min.day
  output$max.day <- max.day
  output$days <- days
  output$n.days <- n.days
  output$rets <- rets
  output$std.days <- std.days
  return(output)
}

# out <- extract.data(PSS_hist = PSS_hist, year = 1995)
# FUNCTION likelihood model for normal distribution
#dist options: pois, norm, flynn, ssq
like.norm.logistic <- function(days, rets, s, m, alpha, sigma, plot, dist='norm') {
  ### TESTING ###
  # rets <- out$rets
  # days <- out$days
  # sigma <- -2.3
  # mu <- 5.25
  # sd <- 12.26
  # alpha <- 20.2
  # plot <- TRUE
  # dist <- 'norm'
  ###############
  if(dist!='norm' & dist!='pois' & dist!='flynn' & dist!='ssq')
  { stop('IMPROPER ERROR DISTRIBUTION SELECTED') }
  m <- exp(m)
  s <- exp(s)
  alpha <- exp(alpha)
  sigma <- exp(sigma)
  
  #days <- out$days
  #rets <- out$rets
  n.days <- length(days)
  
  pred <- vector(length=n.days)
  logLike <- vector(length=n.days)
  
  # d <- 5
  # day <- 152
  for(d in 1:n.days) {
    day <- days[d]
    #Predicted data
    pred[d] <- alpha*(1/(1+exp((-(day-m))/s)))
    
    if(is.na(rets[d])) {  #NO DATA
      logLike[d] <- 0
    }else {
      if(dist=='pois') {
        logLike[d] <- dpois(rets[d],pred[d],log=TRUE)
      }
      if(dist=='norm'){
        logLike[d] <- dnorm(rets[d], pred[d], sigma,log=TRUE)
      }
      if(dist=='flynn') {
        flynn.sigma <- sigma*rets[d] + 0.1
        logLike[d] <- dnorm(rets[d],pred[d], flynn.sigma, log=TRUE)
      }
      if(dist=='ssq') {
        logLike[d] <- -1*(rets[d] - pred[d])^2
      }
    }
  }#next d
  #Plotting
  if(plot==TRUE) {
    y.lim <- c(0,max(rets,pred, na.rm=TRUE))
    # plot(rets ~ days, type='l', col='gray', xlab='Day', ylab='Catch + Esc', lty=3, ylim=y.lim)
    plot(rets ~ days, type='l', col='black', xlab='Day', ylab='Catch + Esc', lty=1, ylim=y.lim, lwd=5)
    points(x=days, y=rets, pch=21, bg='gray')
    #Model data
    lines(x=days, y=pred, lwd=2, col='red')
    mtext(paste('Year = ',year), side=3, outer=TRUE, line=-4, font=2)
  }
  
  #Negative Log Likelihood
  NLL <- -1*sum(logLike)
  return(NLL)
}


####### Normal Curve fitting function #######################
# FUNCTION likelihood model for normal distribution
#dist options: pois, norm, flynn, ssq
like.norm <- function(days, rets, mu, sd, alpha, sigma, plot, dist='norm') {
  ### TESTING ###
  # rets <- out$rets
  # days <- out$days
  # sigma <- -2.3
  # mu <- 5.25
  # sd <- 12.26
  # alpha <- 20.2
  # plot <- TRUE
  # dist <- 'norm'
  ###############
  if(dist!='norm' & dist!='pois' & dist!='flynn' & dist!='ssq')
  { stop('IMPROPER ERROR DISTRIBUTION SELECTED') }
  mu <- exp(mu)
  sd <- exp(sd)
  alpha <- exp(alpha)
  sigma <- exp(sigma)
  
  #days <- out$days
  #rets <- out$rets
  n.days <- length(days)
  
  pred <- vector(length=n.days)
  logLike <- vector(length=n.days)
  
  # d <- 5
  # day <- 152
  for(d in 1:n.days) {
    day <- days[d]
    #Predicted data
    pred[d] <- alpha*dnorm(day,mu,sd)
    
    if(is.na(rets[d])) {  #NO DATA
      logLike[d] <- 0
    }else {
      if(dist=='pois') {
        logLike[d] <- dpois(rets[d],pred[d],log=TRUE)
      }
      if(dist=='norm'){
        logLike[d] <- dnorm(rets[d], pred[d], sigma,log=TRUE)
      }
      if(dist=='flynn') {
        flynn.sigma <- sigma*rets[d] + 0.1
        logLike[d] <- dnorm(rets[d],pred[d], flynn.sigma, log=TRUE)
      }
      if(dist=='ssq') {
        logLike[d] <- -1*(rets[d] - pred[d])^2
      }
    }
  }#next d
  #Plotting
  if(plot==TRUE) {
    y.lim <- c(0,max(rets,pred, na.rm=TRUE))
    #plot(rets ~ days, type='l', col='gray', xlab='Day', ylab='Catch + Esc', lty=3, ylim=y.lim)
    plot(rets ~ days, type='h', col='black', xlab='Day', ylab='Catch + Esc', lty=1, ylim=y.lim, lwd=5)
    #points(x=days, y=rets, pch=21, bg='gray')
    #Model data
    lines(x=days, y=pred, lwd=2, col='red')
  }
  
  #Negative Log Likelihood
  NLL <- -1*sum(logLike)
  return(NLL)
}