// Yukon King Inseason Forecast Version 4.1
//Notes: fitting a normal curve to PSS passage to predict log EOS Can-origin Chinook border passage
//        w/o GSI adjustment


//Components:
//  1) Pre-season forecast  (Prior)
//  
//
//  



data {
  
  // This years Preseason forecast
  real<lower=0> Pf;
  real<lower=0> Pf_sigma;
  

  // EOS Canadian counts by year up to myyear -1
  int<lower=0> n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  // real <lower=0> totalEOS[n_totalEOS];
  
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
  
  // // Matrix of historic counts by days & years adjusted by meanGSI
  // matrix <lower=0> [n_dayPSS, n_yearPSS]PSS_mat_adj;
  
  // Matrix of historic PSS passage for complete historic years
  matrix <lower=0> [n_pss_days_all, n_yearPSS]PSS_mat_all;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // // GSI proportions
  // real <lower=0> meanpropCAN[n_dayPSS];
  // real <lower=0> sd_meanpropCAN[n_dayPSS];
  
  // Estimated Mean, sd, and alphas for normal curve
  real <lower=0> ps_mu[n_yearPSS]; 
  real <lower=0> ps_sd[n_yearPSS];
  real <lower=0> ps_alpha_norm[n_yearPSS];
  
  // //Current years passage adjusted by mean GSI
  // real <lower=0> mean_adj_curr_PSS[n_dayPSS];
  
  
  int <lower=0> myDay;
}

// 
parameters {
    // Normal curve paramaters for current year
  // real <lower=10000, upper=1e6> ps_alpha_curr;
  real <lower=0> ps_alpha_curr;
  real <lower=169, upper=181> ps_mu_curr;
  real <lower=0> ps_sd_curr;
  real <lower=0> sigma;
  
  // Normal curve parameters for histric years
  // real <lower=10000, upper=1e6> ps_alpha_hist[n_yearPSS];
  real <lower=0> ps_alpha_hist[n_yearPSS];
  real <lower=169, upper=181> ps_mu_hist[n_yearPSS];
  real <lower=0> ps_sd_hist[n_yearPSS];
  real <lower=0> sigma_hist[n_yearPSS];
  

  // Regression parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma_reg;
  
  
  
  // real <lower=0> ps_pred_curr;
  // // GSI parameters
  // real <lower=0> phi;
  // // real <lower=0, upper=1> propCAN_1;
  // // real<lower=0, upper=1> propCAN[n_dayPSS];
  // 
  //   // Logit canadian proportions
  // real propCAN_logit[n_dayPSS];
  
    // Empirical sd for update
  // real sigma_pred;
  //
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
  
  // Historic cumulative curve predictions to add to cumHistPSS
  
  vector [n_yearPSS]cum_predicted_histPSS;
  
    // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector <lower=0> [n_yearPSS]ln_predPSS;
  
  // Vector of cumulative complete PSS passage up to myYear-1
  vector [n_yearPSS]cumHistPSS_all;
  
  // // Projection
  // real RunSize;
  
  // PSS prediction from current year
  real ln_curr_predPSS;
  
  // PSS historic prediction for emiricap sd for update
  vector [n_yearPSS]ln_predPSS_hist;
  
  // vector [n_yearPSS]predPSS_hist;
  
  // // Empirical sd for update
  real sigma_pred;
  
  
  // Curve fit to entire current season/year of interest
  ps_pred_curr[1:n_pss_days_all] = ps_alpha_curr * (1/(ps_sd_curr*sqrt(2*pi())))*exp(-0.5*square((pss_days_all - ps_mu_curr)/ps_sd_curr));
  
  
  // Curve fit to days up to myDat for the current year of interest
  ps_pred_curr_obs[1:n_dayPSS] = ps_alpha_curr * (1/(ps_sd_curr*sqrt(2*pi())))*exp(-0.5*square((dayPSS - ps_mu_curr)/ps_sd_curr));
  
  // Curves fit to complete and partial historic passage by year
  for(y in 1:n_yearPSS){
    for(d in 1:n_pss_days_all)
    ps_pred_hist[d,y] = ps_alpha_hist[y] * (1/(ps_sd_hist[y]*sqrt(2*pi())))*exp(-0.5*square((pss_days_all[d] - ps_mu_hist[y])/ps_sd_hist[y]));

    for(d in 1:n_dayPSS){
      ps_pred_hist_obs[d,y] = ps_alpha_hist[y] * (1/(ps_sd_hist[y]*sqrt(2*pi())))*exp(-0.5*square((dayPSS[d] - ps_mu_hist[y])/ps_sd_hist[y]));
    }

  }
  
   
   // Sum up prediction from Normal curve fitting for days myDay +1 to 252
   // Used to calculate empirical sd in update
    for( i in 1:n_yearPSS){
      
      cumHistPSS[i] = sum(PSS_mat[,i]);
    }
    
    cum_predicted_histPSS[1:n_yearPSS] = cumHistPSS + sum(ps_pred_hist[(myDay-147):n_pss_days_all,]);
   
    
    // Sum curr_PSS
    cum_current_PSS = sum(curr_PSS)+sum(ps_pred_curr[myDay-147:n_pss_days_all]);
 
    // Sum up complete PSS histoic days for use in regression
     for( i in 1:n_yearPSS){
      cumHistPSS_all[i] = sum(PSS_mat_all[,i]);
    }
    
    //  Calculate predPSS from aplha, beta and cumHistPSS for model section
    // for (i in 1:n_yearPSS){
    //   ln_predPSS[i] = alpha + beta * cumHistPSS_all[i];}
    ln_predPSS = alpha + beta * cumHistPSS_all;
    
    // Historic prediction
    ln_predPSS_hist = alpha + beta * cum_predicted_histPSS;
    
    // predPSS_hist = exp(ln_predPSS_hist);
    
    // Empirical SD for update  
    sigma_pred = sd(log(totalEOS) - (ln_predPSS_hist));
      // sigma_pred = sd(log(totalEOS ./ cum_predicted_histPSS));
      
      
    // Get current years pss prediction for update
    ln_curr_predPSS = alpha + beta * cum_current_PSS;
 
    // Exponentiate the projection
    // RunSize = exp(ln_RunSize);
   
   
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
  
  // print("ln_predPSS_hist:",ln_predPSS_hist);
  // print("ln_curr_predPSS:",ln_curr_predPSS);
  // print("sigma_pred:",sigma_pred);
}

model {
  //priors
  ps_alpha_curr ~ normal(mean(ps_alpha_norm), sd(ps_alpha_norm)) ;
  // ps_alpha_curr ~ uniform(1000,1e6);
  ps_mu_curr ~ normal(mean(ps_mu), sd(ps_mu)*2);
  ps_sd_curr ~ normal(mean(ps_sd), sd(ps_sd));
  // sigma_day ~ normal(0,10);
  sigma ~ normal(0,5);
  
   ps_alpha_hist ~ normal(mean(ps_alpha_norm), sd(ps_alpha_norm)) ;
  // ps_alpha_curr ~ uniform(1000,1e6);
  ps_mu_hist ~ normal(mean(ps_mu), sd(ps_mu)*2);
  ps_sd_hist ~ normal(mean(ps_sd), sd(ps_sd));
  // sigma_day_hist~ normal(0,10);
  sigma_hist ~ normal(0,5);

  //prior
  alpha ~ normal(0,1e6);
  beta ~ normal(0,10);
  sigma_reg ~ normal(0,5);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  
  for(d in 1:n_dayPSS){
    if(curr_PSS[d]>0){
  log(curr_PSS[d] ) ~ normal(log(ps_pred_curr_obs[d]), sigma);}
  
    for(y in 1:n_yearPSS){
      if(PSS_mat[d,y]>0){
      log(PSS_mat[d,y] ) ~normal(log(ps_pred_hist_obs[d,y]), sigma_hist[y]);}
    }

  }
  
  // Likelihood function for regression
  // log(totalEOS[1:n_yearPSS])~normal(ln_predPSS[1:n_yearPSS], sigma_reg);
  for(i in 1:n_yearPSS){
    
    log(totalEOS[i]) ~ normal((ln_predPSS[i]),sigma_reg);
   
  }

  //Canadian Proportion likelihood
  // // propCAN ~ beta(((1 - meanpropCAN)/(sd_meanpropCAN) - (1/meanpropCAN)) * meanpropCAN^2, ((1/meanpropCAN)-1));
  // for(i in 1:n_dayPSS){
  // propCAN[i] ~ beta(paramA[i], paramB[i]);}
  // // target += beta_lpdf(propCAN|paramA, paramB);
  // 
  
  // // Update for the posterior
  target += normal_lpdf(ln_RunSize | ln_curr_predPSS, sigma_pred);
  // target += normal_lpdf(ln_RunSize | ln_curr_predPSS, sd(log(totalEOS) ./ exp(ln_predPSS_hist)));

}

generated quantities{
  // //
  //
    // Projection
  real RunSize;
  
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real curr_predPSS;
  
  vector [n_yearPSS]predPSS;
  
  // Exponentiate the projection
  RunSize = exp(ln_RunSize);
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);

  // Exponential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);

  predPSS = exp(ln_predPSS);

  curr_predPSS = exp(ln_curr_predPSS);

  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(ln_curr_predPSS, sigma_pred);

  post_curr_predPSS = exp(ln_post_curr_predPSS);



}
