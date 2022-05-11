library(doParallel)
registerDoParallel(cores = 6)
# outputList<-list()
posterior<- foreach(y = testYears)%:% 
  foreach( d = testDays, .packages = "rstan") %dopar%
  {
   InSeasonProjection(model.version = model.version, 
                      myYear = y,
                      myDay = d,
                      n.chains = n.chains,
                      CAN_hist = CAN_hist,
                      pf_hist = pf_hist,
                      PSS_hist = PSS_hist,
                      n.thin = n.thin,
                      n.iter = n.iter )
  }



plot(posterior[[3]][[2]]$pars$RunSize)



outPlots(posterior[[1]][[1]], CAN_hist = CAN_hist)








