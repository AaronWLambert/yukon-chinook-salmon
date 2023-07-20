// Yukon King Inseason Forecast Version 2.c.0.2.sst
//Notes: 


//Components:
//  1) Pre-season forecast  (Prior)
//  2) log(EOS-Can) ~ PSS 
//
//  



data {
  
  // This years Preseason forecast
  real<lower=0> Pf;
  real<lower=0> Pf_sigma;
  int <lower=0> n_yearsPF;

  // EOS Canadian counts excluding myYear
  int<lower=0> n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int n_dayPSS ;
  int<lower=0> dayPSS[n_dayPSS];
  
  // Historic PSS years excluding myYear
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // Matrix of historic counts by days & years
  // int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // Historic PSS years excluding myYear
  int<lower=0> n_yearEagle;
  
  // Matrix of historic counts by days & years
  // int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearEagle]Eagle_mat;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_Eagle;
  vector<lower=0> [n_curr_Eagle]curr_Eagle;
  // Pointer vector
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  
  // April air mean temp
  vector [n_yearPSS]may_sst_hist;
  real  may_sst_curr;
  
  // Propotion of observed Chinook at Eagle
  real mean_propEagle;
  
}

// 
parameters {
   
  
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  real ln_RunSize;
  
}

transformed parameters{
  
   //Empty vector for cum sum of counts for each year to myyear -1
  real cumHistPSS[n_yearPSS];
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]ln_predPSS;
  
  // Predicted Runsize
  real RunSize;
  
  //
  real cum_current_PSS;
  
  //
  real ln_curr_predPSS;
  
  //
  real sigma_predPSS;
  
  //Empty vector for cum sum of counts for historic Eagle passage
  vector [n_yearEagle]cumHistEagle;
  
  // Log prediction of EOS Eagle passage
  vector [n_yearEagle]predEagle;
  
  // Cummulative current Eagle passage
  real cum_current_Eagle;
  
  // log current eagle prediction
  real curr_predEagle;
  
  // Empirical sd for eagle posterior update
  real sigma_predEagle;
  
  
  // Loop to get cumulative counts
 for (i in 1:n_yearPSS){
   cumHistPSS[i] = sum(PSS_mat[,i]);
 }
 
 //  Calculate predPSS from aplha, beta and cumHistPSS for model section
 for (i in 1:n_yearPSS){
 ln_predPSS[i] = alpha + beta * cumHistPSS[i]*may_sst_hist[i];}
 
   // Empirical SD for update  
 sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS]) - ln_predPSS[loc_pf_years_PSS]);
 
 // Bring it back to reality
 RunSize = exp(ln_RunSize);
 
 // Sum curr_PSS
 cum_current_PSS = sum(curr_PSS);
 
 //
 ln_curr_predPSS = alpha + beta * cum_current_PSS * may_sst_curr ;
 
 // Eagle Sonar calculations ////////////////////////////////
 
   // Loop to get cumulative Eagle counts
 for (i in 1:n_yearEagle){
   cumHistEagle[i] = sum(Eagle_mat[,i]);
    }
 
   // Summ current year Eagle passage
 cum_current_Eagle = sum(curr_Eagle);

 // Eagle Component using observed propotions
  
  if(mean_propEagle > 0){
  curr_predEagle = (cum_current_Eagle + 1)/mean_propEagle;
  
  predEagle =  (cumHistEagle+1)/mean_propEagle;
 

  // Empirical SD for update  
 sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS] ./ predEagle[loc_pf_years_Eagle]));}
 
 
 // print("ln_predPSS",ln_predPSS);
}

model {
  //prior
  alpha ~ normal(0,1e6);
  beta ~ normal(0,10);
  sigma ~ normal(0,5);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf-(Pf_sigma^2/2), Pf_sigma);
  
    //Likelihood Function
  for(i in 1:n_yearPSS){
    log(totalEOS[i])~normal(ln_predPSS[i],sigma);
  }
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);

  // Eagle update to the posterior run size projection
  if(mean_propEagle > 0){  
    target += normal_lpdf(ln_RunSize | log(curr_predEagle), sigma_predEagle);}
}

generated quantities{
  //
  real ln_prior_pf;
  
  real prior_pf;
  // 
  vector [n_yearPSS]predPSS;
  //
  real curr_predPSS;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real ln_post_curr_predEagle;
  
  real post_curr_predEagle;
  
  ln_prior_pf = normal_rng(Pf - Pf_sigma^2/2, Pf_sigma);
  
  prior_pf = exp(ln_prior_pf);
  //
  predPSS = exp(ln_predPSS);
  //
  curr_predPSS = exp(ln_curr_predPSS);
  
  ln_post_curr_predPSS = normal_rng(ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  if(mean_propEagle > 0){
  ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);}
  
}






