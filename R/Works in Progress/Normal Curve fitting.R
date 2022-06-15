
library(bbmle)
library(mvtnorm)

extract.data <- function(PSS_hist = PSS_hist, year) {
  ### TESTING ###
  #baywide <- TRUE
  #dist <- 324
  # Year <- 2000
  ###############
    temp.data <- PSS_hist[PSS_hist$Year == year,]
  min.day <- min(temp.data$Day)
  max.day <- max(temp.data$Day)
  days <- min.day:max.day
  n.days <- length(days)
  
  #Calculate standardized date values
  global.days <- min(PSS_hist$Day):max(PSS_hist$Day)
  std.days <- which(global.days %in% days)
  
  rets <- vector(length=n.days)
  
  # d <- 1
  # for(d in 1:n.days) {
  #   day <- days[d]
  #   if(length(temp.data$num[temp.data$jdate==day]) > 0) {
  #     rets[d] <- sum(temp.data$num[temp.data$jdate==day], na.rm=TRUE)*1000
  #   }else {
  #     rets[d] <- NA
  #   }
  # }#next d
  rets <- PSS_hist$count[PSS_hist$Year==year]
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
  if(dist!='norm' & dist!='pois' & dist!='flynn' & dist!='ssq') { stop('IMPROPER ERROR DISTRIBUTION SELECTED') }
  mu <- exp(mu)
  sd <- exp(sd)
  alpha <- exp(alpha)
  sigma <- exp(sigma)
  
  #days <- out$days
  #rets <- out$rets
  n.days <- length(days)
  
  pred <- vector(length=n.days)
  logLike <- vector(length=n.days)
  
  d <- 5
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


#NORMAL PRIOR
norm.prior <- TRUE
if(norm.prior==TRUE) {
  years <- c(1995,1997:2021)
  n.years <- length(years)
  
  norm.mu <- vector(length=n.years)
  norm.sd <- vector(length=n.years)
  norm.alpha <- vector(length=n.years)
  # pdf('Plots/Norm Dist - ssq.pdf', height=8, width=8)
  
  y <- 1
  for(y in 1:n.years) {
    year <- years[y]
    print(year)
    # Extract data
    out <- extract.data(PSS_hist = PSS_hist, year=year)
    #Fit normal data
    fit.norm <-  mle2(like.norm,
                      start=list(mu=log(170), sd=log(20), alpha=log(4000)),
                      fixed=list(plot=FALSE, sigma=log(20)), 
                      data=list(days=out$days, rets=out$rets, dist='ssq'), method='Nelder-Mead', control=list(maxit=100000, trace=FALSE))
    norm.mu[y] <- exp(coef(fit.norm)[1])
    norm.sd[y] <- exp(coef(fit.norm)[2])
    norm.alpha[y] <- exp(coef(fit.norm)[3])
    like.norm(days=out$days, rets=out$rets, (coef(fit.norm)[1]), coef(fit.norm)[2], coef(fit.norm)[3], coef(fit.norm)[4], plot=TRUE, dist='ssq')
    mtext(year, side=3, line=2, font=2)
  }#next y
  
  dev.off()
  prior.df <- data.frame(norm.mu,norm.sd,norm.alpha)
  names(prior.df) <- c('mu','sd','alpha')
  # write.csv(prior.df,'Prior Data/Prior Norm Dist - ssq.csv')
}
prior.df
mean(prior.df$mu)

tt <- vector(length = length(152:220))
d <- 152
for(d in 152:220){
 tt[d] <- norm.alpha[26]*dnorm(d,norm.mu[26],norm.sd[26])}
sum(tt, na.rm = TRUE)
