// Yukon King Inseason Forecast Version 4.2
// Notes: fitting a normal curve to PSS passage to predict EOS Can-origin 
//  Chinook abundance w/ GSI information incorporated


// Components:
// 1) Pre-season forecast  (Prior)
//  
//
//  



data {
  
    // This years Preseason forecast
  real<lower=0> Pf;
  real<lower=0> Pf_sigma;
  int <lower=0> n_yearsPF;

  // EOS Canadian counts by year up to myyear -1
  int<lower=0> n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int n_dayPSS ;
  vector <lower=0> [n_dayPSS]dayPSS;
  
  // Historic PSS years up to myYear -1
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // All the days in PSS season
  int <lower=0> n_pss_days_all;
  vector <lower=0> [n_pss_days_all]pss_days_all;
  
  // Matrix of historic counts by days & years
  // int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Matrix of historic counts by days & years adjusted by meanGSI
  // matrix <lower=0> [n_dayPSS, n_yearPSS]PSS_mat_adj;
  
  // Matrix of historic PSS passage for complete historic years
  matrix <lower=0> [n_pss_days_all, n_yearPSS]PSS_mat_all;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // GSI proportions
  real <lower=0> meanpropCAN_all[n_pss_days_all];
  real <lower=0> sd_meanpropCAN_all[n_pss_days_all];

  //Current years passage adjusted by mean GSI
  real <lower=0> mean_adj_curr_PSS[n_dayPSS];
  
  // Estimated Mean, sd, and alphas for normal curve
  int  <lower=0> n_ps_mu;
  int  <lower=0> n_ps_sd;
  int  <lower=0> n_ps_alpha;
  
  real <lower=0> ps_mu[n_ps_mu]; 
  real <lower=0> ps_sd[n_ps_sd];
  real <lower=0> ps_alpha_norm[n_ps_alpha];
  
  // Pointer vector for PF years in PSS
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  int <lower=0> loc_allDays_myDay; 
  
  int <lower=0> myDay;
  
}

transformed data{
  
   // Parameters for GSI beta likelihood
  real paramA[n_pss_days_all];
  
  real paramB[n_pss_days_all];

   
   // Calculate the A an B parameter for beta dist
  for(x in 1:n_pss_days_all){
    paramA[x] = ((1-meanpropCAN_all[x])/(sd_meanpropCAN_all[x]^2)-(1/meanpropCAN_all[x]))*meanpropCAN_all[x]^2;

    paramB[x] = paramA[x] * ((1/meanpropCAN_all[x])-1);}
}

// 
parameters {
   
  // Normal curve paramaters for current year
  real <lower=0> ps_alpha_curr;
  real <lower=0> ps_mu_curr;
  real <lower=0> ps_sd_curr;
  real <lower=0> sigma;
  
  // Normal curve parameters for histric years
  real <lower=0> ps_alpha_hist[n_yearPSS];
  real <lower=0> ps_mu_hist[n_yearPSS];
  real <lower=0> ps_sd_hist[n_yearPSS];
  real <lower=0> sigma_hist[n_yearPSS];
  
  // Regression parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma_reg;
  
  // // GSI parameters
  real ln_phi;
  // 
  //   // Logit canadian proportions
  real propCAN_logit[n_pss_days_all];
  
  
  // Log runsize projection
  real ln_RunSize;
  
}

transformed parameters{
  
  // Vector for holding this years daily predicted passage
  vector [n_pss_days_all]ps_pred_curr;
  
  // Vector for fitting to all observed days
  vector [n_dayPSS]ps_pred_curr_obs;
  
  // // Matrix for fitting curve to complete historic passage/days
  matrix [n_pss_days_all,n_yearPSS]ps_pred_hist;
  // 
  // // Matrix for fitting curve to historic years up to myDay
  matrix [n_dayPSS, n_yearPSS]ps_pred_hist_obs;
  
  // Vector to hold cumulative sums for historic passage up to myDay 
   //Empty vector for cum sum of counts for each year to myyear -1
  vector [n_yearPSS]cumHistPSS;
  
    // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]ln_predPSS;
  
  // Vector of cumulative complete PSS passage up to myYear-1
  vector [n_yearPSS]cumHistPSS_all;
  
  // PSS prediction from current year
  real ln_curr_predPSS;
  
  // PSS historic log prediction for emirical sd for update
  vector [n_yearPSS]ln_predPSS_hist;
  
  // PSS historic prediction for emiricap sd for update
  // vector [n_yearPSS]predPSS_hist;
  
  // Empirical sd for update
  real sigma_pred;
  
  real phi;
  
  // // Canadian proportions after inv-logit transformation
  vector<lower=0, upper=1> [n_pss_days_all]propCAN;
 
 
  // Exponentiate the variance term in the GSI random walk
  phi = exp(ln_phi);
  
  //  Bring propCAN_logit into beta space
  for(t in 1:n_pss_days_all){
    propCAN[t] = exp(propCAN_logit[t])/(1+exp(propCAN_logit[t]));}
    
  // Curve fit to entire current season/year of interest
  ps_pred_curr[1:n_pss_days_all] =  ps_alpha_curr * (1/(ps_sd_curr*sqrt(2*pi())))*exp(-0.5*square((pss_days_all - ps_mu_curr)/ps_sd_curr));
  
  
  // Curve fit to days up to myDay for the current year of interest
  ps_pred_curr_obs[1:n_dayPSS] = ps_alpha_curr * (1/(ps_sd_curr*sqrt(2*pi())))*exp(-0.5*square((dayPSS - ps_mu_curr)/ps_sd_curr));
  

  for(y in 1:n_yearPSS){
    for(d in 1:n_pss_days_all)
    ps_pred_hist[d,y] =  ps_alpha_hist[y] * (1/(ps_sd_hist[y]*sqrt(2*pi())))*exp(-0.5*square((pss_days_all[d] - ps_mu_hist[y])/ps_sd_hist[y]));

    for(d in 1:n_dayPSS){
      ps_pred_hist_obs[d,y] = ps_alpha_hist[y] * (1/(ps_sd_hist[y]*sqrt(2*pi())))*exp(-0.5*square((dayPSS[d] - ps_mu_hist[y])/ps_sd_hist[y]));
    }

  }
   
   // Sum up prediction from Normal curve fitting for days myDay 
   // Used to calculate empirical sd in update
     // Loop to get cumulative counts
    for (i in 1:n_yearPSS){
   
       cumHistPSS[i] = 0;
   
         for(j in 1:n_dayPSS){
     
             cumHistPSS[i] += (PSS_mat[j,i] * propCAN[j]);
                  }
                      }
     // Loop to get cumulative counts
    for (i in 1:n_yearPSS){
   
        for(j in (n_dayPSS + 1):n_pss_days_all){
     
            cumHistPSS[i] += (ps_pred_hist[j,i] * propCAN[j]);
                }
                }
   
     // Loop to get cumulative counts
   // for (i in 1:n_yearPSS){
   //   cumHistPSS[i] = 0;
   // 
   //  // for(j in 1:n_dayPSS){
   //  //   cumHistPSS[i] += ((PSS_mat[j,i] * propCAN[j]) + (ps_pred_hist[j+(225-myDay),i]*propCAN[j+(225-myDay)]));
   //  // }
   //  //   }
   //  for(j in 1:n_dayPSS){
   //    cumHistPSS[i] += ((PSS_mat[j,i] * propCAN[j]) + (ps_pred_hist[j+(225-myDay),i]*propCAN[j+(225-myDay)]));
   //  }
   //    }
    
    
    
    // 
    // cum_predicted_histPSS[1:n_yearPSS] = cumHistPSS + sum(ps_pred_hist[(myDay-147):n_pss_days_all,]);
    
   
    // Sum curr_PSS
    // cum_current_PSS = sum(curr_PSS)+sum(ps_pred_curr[myDay-147:n_pss_days_all]);
    
    cum_current_PSS = 0;
    for(j in 1:n_dayPSS){
      cum_current_PSS += (curr_PSS[j] * propCAN[j]);
      }
    
    for(j in (n_dayPSS + 1):n_pss_days_all){
    cum_current_PSS += ps_pred_curr[j]*propCAN[j];}
      
     
    // Sum historic PSS observed after estimated adjustment for regression relationship
       for (i in 1:n_yearPSS){
     cumHistPSS_all[i] = 0;
   
    for(j in 1:n_pss_days_all){
      cumHistPSS_all[i] += (PSS_mat_all[j,i] * propCAN[j]);
    }
      }
    
    
    //  Regression fitting for EOS vs Run size
    for (i in 1:n_yearPSS){
      ln_predPSS[i] = alpha + beta * cumHistPSS_all[i];}
      
     // Historic prediction for empirical sd caluclation
    ln_predPSS_hist[1:n_yearPSS] = alpha + beta * cumHistPSS;
    
    // predPSS_hist = exp(ln_predPSS_hist);
    
    // Empirical SD
    sigma_pred = sd(log(totalEOS[loc_pf_years_PSS]) - ln_predPSS_hist[loc_pf_years_PSS]);
    
    // Get current years pss prediction for update
    ln_curr_predPSS = alpha + beta * cum_current_PSS;
 

   
   
 // Print statements used to debug script
 // print("ps_mu=",ps_mu);
 // print("ps_sd=",ps_sd);
 // print("alpha=",ps_alpha);
 // print("ps_pred_curr=",ps_pred_curr);
 // print("ps_pred_curr_obs=",ps_pred_curr_obs);
 // print("ps_pred_hist=",ps_pred_hist);
 // print("ps_pred_curr_alpha=",ps_pred_curr_alpha);
  // print("cum_predicted_histPSS=",cum_predicted_histPSS)
  // print("cum_current_PSS:", cum_current_PSS);
  // print("predPSS:",predPSS);
  // print("sigma_pred:",sigma_pred);
}

model {
  //priors
  ps_alpha_curr ~ normal(mean(ps_alpha_norm), sd(ps_alpha_norm)) ;
  // ps_alpha_curr ~ uniform(19000,112000);
  ps_mu_curr ~ normal(mean(ps_mu), sd(ps_mu));
  ps_sd_curr ~ normal(mean(ps_sd), sd(ps_sd));
  // sigma_day ~ normal(0,10);
  sigma ~ normal(0,2);
  
   ps_alpha_hist ~ normal(mean(ps_alpha_norm), sd(ps_alpha_norm)) ;
  // ps_alpha_curr ~ uniform(19000,112000);
  ps_mu_hist ~ normal(mean(ps_mu), sd(ps_mu));
  ps_sd_hist ~ normal(mean(ps_sd), sd(ps_sd));
  // sigma_day_hist~ normal(0,10);
  sigma_hist ~ normal(0,2);

  //prior
  alpha ~ normal(0,1e6);
  beta ~ normal(0,10);
  sigma_reg ~ normal(0,5);
  
  // Prior for GSI proportion
  ln_phi ~ normal(-1,1);
  // Draw a GSI value for day 1 of passage
  propCAN_logit[1] ~ normal(0,1);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf - Pf_sigma^2/2, Pf_sigma);
  
  // Likelihood for fitting normal curve to observed passage
  for(d in 1:n_dayPSS){
    if(curr_PSS[d]>0){
  log(curr_PSS[d] ) ~ normal(log(ps_pred_curr_obs[d]), sigma);}
  
    for(y in 1:n_yearPSS){
      if(PSS_mat[d,y]>0){
      log(PSS_mat[d,y]) ~normal(log(ps_pred_hist_obs[d,y]), sigma_hist[y]);}
    }

  }
  
  
  // Random walk
  for(d in 2:n_pss_days_all){
   propCAN_logit[d] ~ normal(propCAN_logit[d-1], phi);
  }


  // Likelihood for regression
  for(i in 1:n_yearPSS){
    log(totalEOS)[i] ~ normal(ln_predPSS[i], sigma_reg);
  }
  

  // Canadian Proportion likelihood
  for(i in 1:n_pss_days_all){
  propCAN[i] ~ beta(paramA[i], paramB[i]);}
  
  // Update runsize for the posterior
  target += normal_lpdf(ln_RunSize | ln_curr_predPSS - sigma_pred^2/2, sigma_pred);

}

generated quantities{
  // //
  //
  real RunSize;
  real ln_prior_pf;
  real prior_pf;
  
  real ln_post_curr_predPSS;
  real post_curr_predPSS;
  
  // Exponentiate the projection
  RunSize = exp(ln_RunSize);
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf - Pf_sigma^2/2, Pf_sigma);

  // Exponential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);


  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(ln_curr_predPSS, sigma_pred);

  post_curr_predPSS = exp(ln_post_curr_predPSS);

}






