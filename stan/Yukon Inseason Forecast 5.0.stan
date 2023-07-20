// Yukon King Inseason Forecast Version 5.0
//Notes: fitting a logistiv curve to PSS passage to predict EOS Can-origin Chinook border passage


//Components:
//  1) Pre-season forecast  (Prior)
//  2) PSS log curve and EOS regression
//  3) Eagle proportion estimator
//  

// Next to address:
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
  
    // Matrix of cumulative historic PSS passage by days and years
  matrix <lower=0> [n_dayPSS, n_yearPSS]cum_PSS_mat;
  
  // Matrix of historic counts by days & years adjusted by meanGSI
  // matrix <lower=0> [n_dayPSS, n_yearPSS]PSS_mat_adj;
  
  // Matrix of historic PSS passage for complete historic years
  matrix <lower=0> [n_pss_days_all, n_yearPSS]PSS_mat_all;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
    // Cumulative current year counts up to myDay
  real<lower=0> cum_curr_PSS[n_dayPSS];
  
  // Estimated Mean, sd, and alphas for normal curve
    int  <lower=0> n_ps_m;
  int  <lower=0> n_ps_s;
  int  <lower=0> n_ps_alpha_log;
  
  real <lower=0> ps_m[n_ps_m]; 
  real <lower=0> ps_s[n_ps_s];
  real <lower=0> ps_alpha_log[n_ps_alpha_log];
  
  // Historic PSS years excluding myYear
  int<lower=0> n_yearEagle;
  
  // Matrix of historic counts by days & years
  matrix<lower=0> [n_dayPSS, n_yearEagle]Eagle_mat;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_Eagle;
  vector<lower=0> [n_curr_Eagle]curr_Eagle;
  
  // Location vector of integers
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  int <lower=0> loc_allDays_myDay; 
  
  // Propotion of observed fish at Eagle
  real mean_propEagle;
  
  int <lower=0> myDay;
}

// 
parameters {
   
  // logistic curve paramaters for current year
  // real <lower=10000, upper=1e6> ps_alpha_curr;
  real <lower=0> ps_alpha_curr;
  real <lower=0> ps_mid_curr;
  // real<lower=0> ps_mid_curr;
  real <lower=0> ps_shape_curr;
  real <lower=0> sigma;
  
  // Logistic curve parameters for histric years
  // real <lower=10000, upper=1e6> ps_alpha_hist[n_yearPSS];
  vector <lower=0> [n_yearPSS]ps_alpha_hist;
  real <lower=0>ps_mid_hist[n_yearPSS];
  // real<lower=0> ps_mid_hist[n_yearPSS];
  real <lower=0> ps_shape_hist[n_yearPSS];
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
  
  // // Historic cumulative curve predictions to add to cumHistPSS
  // 
  vector [n_yearPSS]cum_predicted_histPSS;
  
    // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]predPSS;
  
  // Vector of cumulative complete PSS passage up to myYear-1
  vector [n_yearPSS]cumHistPSS_all;
  
  // Projection
  // real RunSize;
  
  // PSS prediction from current year
  real curr_predPSS;
  
  // PSS historic prediction for emirical sd for update
  vector [n_yearPSS]predPSS_hist;
  
  // Empirical sd for update
  real sigma_pred;
  
  
  // // Curve fit to entire current season/year of interest
  // ps_pred_curr[1:n_pss_days_all] = ps_alpha_curr * (1/(ps_sd_curr*sqrt(2*pi())))*exp(-0.5*square((pss_days_all - ps_mu_curr)/ps_sd_curr));
  // 
  for(d in 1:n_pss_days_all){
   ps_pred_curr[d] = ps_alpha_curr *(1/(1 + exp((-(pss_days_all[d] - ps_mid_curr)) / ps_shape_curr)));}
   // ps_pred_curr = ps_alpha_curr *(1/(1 + exp((-(pss_days_all - ps_mid_curr)) / ps_shape_curr)));
  // Curve fit to days up to myDat for the current year of interest
  for(d in 1:n_dayPSS){
  ps_pred_curr_obs[d] = ps_alpha_curr * (1/(1 + exp((-(dayPSS[d] - ps_mid_curr))/ps_shape_curr)));}
  
  // Curves fit to complete and partial historic passage by year
  for(y in 1:n_yearPSS){
    for(d in 1:n_pss_days_all)
    ps_pred_hist[d,y] = ps_alpha_hist[y] * (1/(1 + exp((-(pss_days_all[d] - ps_mid_hist[y]))/ps_shape_hist[y])));

    for(d in 1:n_dayPSS){
      ps_pred_hist_obs[d,y] = ps_alpha_hist[y] * (1/(1 + exp((-(dayPSS[d] - ps_mid_hist[y]))/ps_shape_hist[y])));
    }

  }
  
 
    // Sum up complete PSS histoic days for use in regression
     for( i in 1:n_yearPSS){
      cumHistPSS_all[i] = sum(PSS_mat_all[,i]);
    }
    

    predPSS = alpha +beta * cumHistPSS_all;
    
    // // Historic prediction
    // predPSS_hist = alpha + beta * ps_alpha_hist;
    
        for( i in 1:n_yearPSS){
      
      cumHistPSS[i] = sum(PSS_mat[,i]);
    }
    
    for(y in 1:n_yearPSS){
    cum_predicted_histPSS[y] = cumHistPSS[y] + (ps_pred_hist[n_pss_days_all,y]-ps_pred_hist[(n_dayPSS+1),y]);}
    
    // Historical prediction for empirial sd
    predPSS_hist = alpha + beta * cum_predicted_histPSS;
    
    // Empirical SD for update  
    sigma_pred = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS_hist[loc_pf_years_PSS]));
      // sigma_pred = sd(log(totalEOS ./ cum_predicted_histPSS));
      // sigma_pred = sd(log(totalEOS) - (ln_predPSS_hist));
      
    // Sum curr_PSS
    cum_current_PSS = sum(curr_PSS) + (ps_pred_curr[n_pss_days_all] - ps_pred_curr[n_dayPSS+1]);  
      
    // Get current years pss prediction for update
    curr_predPSS = alpha + beta * cum_current_PSS;
 
    // // Exponentiate the projection
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
  // print("predPSS:",predPSS);
  // print("sigma_pred:",sigma_pred);
}

model {
  //priors
  ps_alpha_curr ~ normal(mean(ps_alpha_log), sd(ps_alpha_log)*2);
  // ps_alpha_curr ~ uniform(30000,400000);
  ps_mid_curr ~ normal(mean(ps_m), sd(ps_m)/2);
  // ps_mid_curr ~uniform(167,182);
  ps_shape_curr ~ normal(mean(ps_s), sd(ps_s)*2);
  sigma ~ normal(0,2);
  
   ps_alpha_hist ~ normal(mean(ps_alpha_log), sd(ps_alpha_log)*2);
  // ps_alpha_hist ~ uniform(30000,400000);
  ps_mid_hist ~ normal(mean(ps_m), sd(ps_m)/2);
  // ps_mid_hist ~uniform(167,182);
  ps_shape_hist ~ normal(mean(ps_s), sd(ps_s)*2);
  sigma_hist ~ normal(0,2);

  //prior
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma_reg ~ normal(0,10);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  
  for(d in 1:n_dayPSS){ 
    if(cum_curr_PSS[d] > 0){
      
  log(cum_curr_PSS[d] ) ~ normal(log(ps_pred_curr_obs[d]), sigma);}

    for(y in 1:n_yearPSS){
    if(cum_PSS_mat[d,y] > 0 ){
      
      log(cum_PSS_mat[d,y]) ~normal(log(ps_pred_hist_obs[d,y]), sigma_hist[y]);}
    }

  }
  
  // Likelihood function for regression
  log(totalEOS[1:n_yearPSS])~normal(log(predPSS[1:n_yearPSS]), sigma_reg);

  //Canadian Proportion likelihood
  // // propCAN ~ beta(((1 - meanpropCAN)/(sd_meanpropCAN) - (1/meanpropCAN)) * meanpropCAN^2, ((1/meanpropCAN)-1));
  // for(i in 1:n_dayPSS){
  // propCAN[i] ~ beta(paramA[i], paramB[i]);}
  // // target += beta_lpdf(propCAN|paramA, paramB);
  // 
  // // Update for the posterior
  target += normal_lpdf(ln_RunSize | log(curr_predPSS), sigma_pred);

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
  
  // Exponentiate the projection
  RunSize = exp(ln_RunSize);
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);

  // Exponential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);

  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_pred);

  post_curr_predPSS = exp(ln_post_curr_predPSS);


}









