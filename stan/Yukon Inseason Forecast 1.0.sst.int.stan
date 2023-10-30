// Yukon King Inseason Forecast Version 1.0 with sst interaction only
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

  // EOS Canadian counts by year up to myyear -1
  int<lower=0> n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int n_dayPSS ;
  int<lower=0> dayPSS[n_dayPSS];
  
  // Historic PSS years up to myYear -1
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // Matrix of historic counts by days & years
  // int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
      // Historic PSS years up to myYear -1
  int<lower=0> n_yearEagle;
  
  // Matrix of historic counts by days & years
  // int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearEagle]Eagle_mat;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_Eagle;
  vector<lower=0> [n_curr_Eagle]curr_Eagle;
  
  // Location vector of integers
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  int <lower=0> loc_devMP_allYears[n_yearEagle];
  
  // Propotion of observed fish at Eagle
  // vector [n_yearEagle]propEagle;
  real mean_propEagle;
  
  // April air mean temp
  int <lower=0> n_yearSST;
  vector [n_yearSST]may_sst_hist;
  real  may_sst_curr;
}

// 
parameters {
   
  
  real <lower=0> alpha;
  real  beta;
  real  beta2;
  real <lower=0> sigma;
  
  real ln_RunSize;
  
}

transformed parameters{
  
   //Empty vector for cum sum of counts for each year to myyear -1
  real cumHistPSS[n_yearPSS];
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]predPSS;
  
  //
  real cum_current_PSS;
  
  //
  real curr_predPSS;
  
  //
  real sigma_predPSS;
    
  //Empty vector for cum sum of counts for each year to myyear -1
  vector [n_yearEagle]cumHistEagle;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearEagle]predEagle;  
  real cum_current_Eagle;
  
  //
  real curr_predEagle;
  
  real sigma_predEagle;
  
  vector [n_yearPSS]may_sst_hist2;
  
  // Loop to get cumulative counts
 for (i in 1:n_yearPSS){
   cumHistPSS[i] = sum(PSS_mat[,i]);
 }
 
 // Extract the years of sst observations that are used (2005 onward, excluding myYear)
 may_sst_hist2 = may_sst_hist[loc_devMP_allYears];
 
 //  Calculate predPSS from aplha, beta and cumHistPSS for model section
 for (i in 1:n_yearPSS){
 predPSS[i] = alpha + beta * cumHistPSS[i] + beta2 * cumHistPSS[i] * may_sst_hist2[i];}
 
   // Empirical SD for update  
 sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS[loc_pf_years_PSS]));
 
 // Sum curr_PSS
 cum_current_PSS = sum(curr_PSS);
 
 //
 curr_predPSS = alpha + beta * cum_current_PSS * may_sst_curr ;
 
 // Loop to get cumulative Eagle counts
  for (i in 1:n_yearEagle){
     cumHistEagle[i] = sum(Eagle_mat[,i]);
    }
 
  // Summ current year Eagle passage
   cum_current_Eagle = sum(curr_Eagle);
 
   // Eagle Component using observed propotions
  if(mean_propEagle > 0 && cum_current_Eagle > 0){
  curr_predEagle = (cum_current_Eagle+1)/mean_propEagle;
  
  predEagle =  (cumHistEagle+1)/mean_propEagle;
 

  // Empirical SD for update  
 sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS] ./ predEagle[loc_pf_years_Eagle]));}
 
 // print("predPSS:",predPSS);
}

model {
  //prior
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  beta2 ~normal(0,1e15);
  sigma ~ normal(0,10);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
    //Likelihood Function
  for(i in 1:n_yearPSS){
    
    log(totalEOS[i]) ~ normal(log(predPSS[i]),sigma);
  }
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | log(curr_predPSS), (sigma_predPSS));
  
  // Update likelihood for the posterior prediction
  if(mean_propEagle > 0 && cum_current_Eagle > 0){
  target += normal_lpdf(ln_RunSize | log(curr_predEagle), sigma_predEagle);}

}

generated quantities{
  //
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real ln_post_curr_predEagle;
  
  real post_curr_predEagle;
  
  // Predicted Runsize
  real RunSize;
  
  ln_prior_pf = normal_rng(Pf, Pf_sigma);
  
  prior_pf = exp(ln_prior_pf);
  
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  if(mean_propEagle > 0 && cum_current_Eagle > 0){
  ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);}
  
  // Bring it back to reality
  RunSize = exp(ln_RunSize);
  
}






