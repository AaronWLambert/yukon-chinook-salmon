#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.14.22
# Version: This function will run any version of the Yukon Inseason Forecast model outputlist
# Purpose: Functions that will calculate retrospective testing stats
# 
#
#
#
#
#=================================================================================
#NOTES:
# This function will be called from another script.
# 
# Next steps: 
# More functions? 
#
#=================================================================================

# Function for MAPE, PE, RMSE, and other things listed in the function description

retrospective.function<- function(outputList, 
                                  testYears, 
                                  testDays, 
                                  CAN_hist,
                                  pf = FALSE){
  # OutputList <- a list of outputs from the Stan model runs for days and years
  # testYears <- vector of days included in OutputList
  # testDays <- vector of years included in OutputList
  # CAN_hist <- dataframe of historic EOS reconstructed run sizes
  # Returns a list including:
  #     Matrix of medians from posterior RunSize projections
  #     Matrix of historic EOS Reconstructed Can-orig Runsize
  #     Percent error matrix
  #     APE matrix
  #     MAPE for days across years vector
  #     RMSE for days across years vector
  #     RMSE for years across days vector
  #     Matrices for how often the "true" EOS can-orig runsize falls in the 
  #      50, 80, and 95% intervals.
  #
  # Calculations for retrospecitve testing #############################################
  if(pf == FALSE){
  # Create matrix of observed and median outputs
  # 
  # Create matrix with dimensions of days by years that are in the outputList from model run
  median_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
  
  # Give names to matrix
  colnames(median_mat) <- c(testYears)
  rownames(median_mat) <- c(testDays)
  
  # Make sure it works...
  dim(median_mat)
  str(median_mat)
  
  # Extract medians from posterior RunSize prediction and place in vector for matrix pop.
  # populate matrix with medians from posterior RunSize prediction
  counter <- 1
  for (y in 1:length(testYears)) {
    for (d in 1:length(testDays)) {
      
      median_mat[d,y]<- median(outputList[[counter]]$pars$RunSize)
      
      counter <- counter+1
    }
  }
  median_mat
  
  # Matrix with realized reconstructed counts by year and day 
  
  realized_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
  
  # Give names to matrix
  colnames(realized_mat) <- c(testYears)
  rownames(realized_mat) <- c(testDays)
  
  # Make sure it works...
  dim(realized_mat)
  
  # Extract values for realized runsize
  realized_vect <- CAN_hist$can.mean[CAN_hist$Year >= 2007&
                                       CAN_hist$Year != 2011 & 
                                       CAN_hist$Year != 2012]
  
  counter <- 1
  # Populate the realized matrix
  for (y in 1:length(testYears)) {
    for (d in 1:length(testDays)) {
      
      realized_mat[d,y]<- realized_vect[counter]
      
    }
    counter <- counter+1
  }
  realized_mat
  # Percent Error (PE) matrix
  PE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
  
  # Give names to matrix
  colnames(PE_mat) <- c(testYears)
  rownames(PE_mat) <- c(testDays)
  
  # Populate the PE matrix
  for (y in 1:length(testYears)) {
    for (d in 1:length(testDays)) {
      
      PE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
      
    }
  }
  # Percent Error (APE) matrix
  APE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
  
  # Give names to matrix
  colnames(APE_mat) <- c(testYears)
  rownames(APE_mat) <- c(testDays)
  
  # Populate the absolute PE matrix
  for (y in 1:length(testYears)) {
    for (d in 1:length(testDays)) {
      
      APE_mat[d,y]<- abs(realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
      
    }
    
  }
  
  
  # Calculate MAPE by Days
  MAPE_day_vect <- vector(length = length(PE_mat[,1]))
  names(MAPE_day_vect) <- testDays
  
  for (d in 1:length(MAPE_day_vect)) {
    MAPE_day_vect[d] <- mean(APE_mat[d,])
  }
  
  # Root Mean Square error for days across years
  # Square Error (SE) matrix
  SE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
  
  # Give names to matrix
  colnames(SE_mat) <- c(testYears)
  rownames(SE_mat) <- c(testDays)
  
  # Populate the SE matrix
  for (y in 1:length(testYears)) {
    for (d in 1:length(testDays)) {
      
      SE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])^2
      
    }
    
  }
  
  # Calculate RMSE by days
  RMSE_day_vect <- vector(length = length(SE_mat[,1]))
  names(RMSE_day_vect) <- testDays
  
  for (d in 1:length(RMSE_day_vect)) {
    RMSE_day_vect[d] <- sqrt(mean(SE_mat[d,]))
  }
  
  # calculate RMSE by year
  RMSE_year_vect <- vector(length = length(SE_mat[1,]))
  names(RMSE_year_vect) <- as.character(testYears)
  
  for (y in 1:length(RMSE_year_vect)) {
    RMSE_year_vect[y] <- sqrt(mean(SE_mat[,y]))
  }
  
  # How often does the interval capture the true return size? ######
  # 95% Interval
  true95Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
  true80Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
  true50Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
  # Give names to matrix
  colnames(true95Mat) <- c(testYears)
  rownames(true95Mat) <- c(testDays)
  
  colnames(true80Mat) <- c(testYears)
  rownames(true80Mat) <- c(testDays)
  
  colnames(true50Mat) <- c(testYears)
  rownames(true50Mat) <- c(testDays)
  
  count <- 1
  
  for (y in 1:length(testYears)) {
    for (d in 1:length(testDays)) {
      
      # 95%
      test95<-quantile(outputList[[count]]$pars$RunSize, probs = c(.025,.975))
      
      true95Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test95[2] & 
                                 CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test95[1],
                               1,
                               0)
      # 80%
      test80<-quantile(outputList[[count]]$pars$RunSize, probs = c(.1,.9))
      
      true80Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test80[2] & 
                                 CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test80[1],
                               1,
                               0)
      # 50%
      test50<-quantile(outputList[[count]]$pars$RunSize, probs = c(.25,.75))
      
      true50Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test50[2] & 
                                 CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test50[1],
                               1,
                               0)
      count <- count+1
    }
  }
  
  }else{
    # Create matrix with dimensions of days by years that are in the outputList from model run
    median_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
    
    # Give names to matrix
    colnames(median_mat) <- c(testYears)
    rownames(median_mat) <- c(testDays)
    
    # Make sure it works...
    dim(median_mat)
    str(median_mat)
    
    # Extract medians from posterior RunSize prediction and place in vector for matrix pop.
    # populate matrix with medians from posterior RunSize prediction
    counter <- 1
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        median_mat[d,y]<- median(outputList[[counter]]$pars$prior_pf )
        
        counter <- counter+1
      }
    }
    median_mat
    
    # Matrix with realized reconstructed counts by year and day 
    
    realized_mat <- matrix(nrow = length(testDays),ncol = length(testYears))
    
    # Give names to matrix
    colnames(realized_mat) <- c(testYears)
    rownames(realized_mat) <- c(testDays)
    
    # Make sure it works...
    dim(realized_mat)
    
    # Extract values for realized runsize
    realized_vect <- CAN_hist$can.mean[CAN_hist$Year >= 2007 &
                                         CAN_hist$Year != 2011 & 
                                         CAN_hist$Year != 2012]
    
    counter <- 1
    # Populate the realized matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        realized_mat[d,y]<- realized_vect[counter]
        
      }
      counter <- counter+1
    }
    realized_mat
    # Percent Error (PE) matrix
    PE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(PE_mat) <- c(testYears)
    rownames(PE_mat) <- c(testDays)
    
    # Populate the PE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        PE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
        
      }
    }
    # Percent Error (APE) matrix
    APE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(APE_mat) <- c(testYears)
    rownames(APE_mat) <- c(testDays)
    
    # Populate the absolute PE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        APE_mat[d,y]<- abs(realized_mat[d,y]- median_mat[d,y])/realized_mat[d,y]
        
      }
      
    }
    
    
    # Calculate MAPE by Days
    MAPE_day_vect <- vector(length = length(PE_mat[,1]))
    names(MAPE_day_vect) <- testDays
    
    for (d in 1:length(MAPE_day_vect)) {
      MAPE_day_vect[d] <- mean(APE_mat[d,])
    }
    
    # Root Mean Square error for days across years
    # Square Error (SE) matrix
    SE_mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    
    # Give names to matrix
    colnames(SE_mat) <- c(testYears)
    rownames(SE_mat) <- c(testDays)
    
    # Populate the SE matrix
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        SE_mat[d,y]<- (realized_mat[d,y]- median_mat[d,y])^2
        
      }
      
    }
    
    # Calculate RMSE by days
    RMSE_day_vect <- vector(length = length(SE_mat[,1]))
    names(RMSE_day_vect) <- testDays
    
    for (d in 1:length(RMSE_day_vect)) {
      RMSE_day_vect[d] <- sqrt(mean(SE_mat[d,]))
    }
    
    # calculate RMSE by year
    RMSE_year_vect <- vector(length = length(SE_mat[1,]))
    names(RMSE_year_vect) <- as.character(testYears)
    
    for (y in 1:length(RMSE_year_vect)) {
      RMSE_year_vect[y] <- sqrt(mean(SE_mat[,y]))
    }
    
    # How often does the interval capture the true return size? ######
    # 95% Interval
    true95Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    true80Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    true50Mat <- matrix(nrow = length(testDays), ncol = length(testYears))
    # Give names to matrix
    colnames(true95Mat) <- c(testYears)
    rownames(true95Mat) <- c(testDays)
    
    colnames(true80Mat) <- c(testYears)
    rownames(true80Mat) <- c(testDays)
    
    colnames(true50Mat) <- c(testYears)
    rownames(true50Mat) <- c(testDays)
    
    count <- 1
    
    for (y in 1:length(testYears)) {
      for (d in 1:length(testDays)) {
        
        # 95%
        test95<-quantile(outputList[[count]]$pars$prior_pf, probs = c(.025,.975))
        
        true95Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test95[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test95[1],
                                 1,
                                 0)
        # 80%
        test80<-quantile(outputList[[count]]$pars$prior_pf, probs = c(.1,.9))
        
        true80Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test80[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test80[1],
                                 1,
                                 0)
        # 50%
        test50<-quantile(outputList[[count]]$pars$prior_pf, probs = c(.25,.75))
        
        true50Mat[d,y] <-ifelse( CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] < test50[2] & 
                                   CAN_hist$can.mean[CAN_hist$Year == outputList[[count]]$myYear] > test50[1],
                                 1,
                                 0)
        count <- count+1
      }
    }
  } 
  
  RetroList <- list("Med_mat" = median_mat, 
                    "Realized_mat" = realized_mat,
                    "PE_mat" = PE_mat,
                    "APE_mat" = APE_mat,
                    "MAPE_vect" = MAPE_day_vect,
                    "RMSE_by_day_vect" = RMSE_day_vect,
                    "RMSE_by_year_vect" = RMSE_year_vect,
                    "True95" = true95Mat,
                    "True80" = true80Mat,
                    "True50" = true50Mat,
                    "Version" = outputList$`2013_153`$version,
                    "TestDaysIncluded" <- testDays,
                    "TestYearsIncluded" <- testYears)
  
  return(RetroList)
}

# Retro Plots function ###########################################################3

