// Yukon King Inseason Forecast version 1.0.2
//Notes: Regression in normal space with eagle added into likelihood


//Components:
//  1) Pre-season forecast  (Prior)
//  2) EOS-Can ~ PSS in normal space
//  3) Eagle Sonar regression added to likelihood
//  



data {
  
  // This years Preseason forecast
  real <lower=0>Pf;
  real <lower=0>Pf_sigma;
  int <lower=0> n_yearsPF;
  

  // EOS Canadian counts by year up to myyear -1
  int <lower=0>n_totalEOS ;
  vector<lower=0> [n_totalEOS]totalEOS;
  
  // PSS days up to myDay
  int<lower=0> n_dayPSS ;
  int<lower=0> dayPSS[n_dayPSS];
  
  // Historic PSS years up to myYear -1
  int<lower=0> n_yearPSS;
  int <lower=0>yearPSS[n_yearPSS];
  
  // Matrix of historic counts by days & years
  // int <lower=0>n_hist_counts;
  matrix <lower=0>[n_dayPSS, n_yearPSS]PSS_mat;
  
  // Current year counts up to myDay
  int <lower=0> n_curr_PSS;
  real <lower=0> curr_PSS[n_curr_PSS];
  
  // Mean propotion of PSS of EOS can passage
  real mean_propPSS;
  
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
  
  // Propotion of observed fish at Eagle
  // vector [n_yearEagle]propEagle;
  real mean_propEagle;
  
}

// 
parameters {
  
  real ln_RunSize;
  
}

transformed parameters{
  
   //Empty vector for cum sum of counts for each year to myyear -1
  vector [n_yearPSS]cumHistPSS;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]predPSS;
  
  // //Empty vector for cum sum of counts for each year to myyear -1
  // vector [n_yearEagle]cumHistEagle;
  // 
  // // Empty vector for predicted PSS counts in transformed parameters
  // vector [n_yearEagle]predEagle;
  
  // Variable for posterior abundance
  // real RunSize;
  
  // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  //
  real curr_predPSS;
  
  //
  // real cum_current_Eagle;
  // 
  // //
  // real curr_predEagle;
  // 
  real sigma_predPSS;
  // 
  // real sigma_predEagle;
  
  // Loop to get cumulative counts
  for (i in 1:n_yearPSS){
     cumHistPSS[i] = sum(PSS_mat[,i]);
    }
 
  // // Loop to get cumulative Eagle counts
  // for (i in 1:n_yearEagle){
  //    cumHistEagle[i] = sum(Eagle_mat[,i]);
  //   }
 
  // Summ current year Eagle passage
   // cum_current_Eagle = sum(curr_Eagle);
   
  //  Calculate predPSS from aplha, beta and cumHistPSS for model section
  if(cumHistPSS[n_yearPSS]>0){
  for (i in 1:n_yearPSS){
      
    predPSS[i] = (cumHistPSS[i] + 1) / mean_propPSS;

    }
  
    sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS[loc_pf_years_PSS]));
 
  // //
  //   RunSize = exp(ln_RunSize);
 
  // Sum curr_PSS
    cum_current_PSS = sum(curr_PSS);
 
  //
    curr_predPSS = cum_current_PSS/mean_propPSS;}
 
  // Eagle Component using observed propotions
 //  if(mean_propEagle > 0){
 //  curr_predEagle = (cum_current_Eagle+1)/mean_propEagle;
 //  
 //  predEagle =  (cumHistEagle+1)/mean_propEagle;
 // 
 // 
 //  // Empirical SD for update  
 // sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS] ./ predEagle[loc_pf_years_Eagle]));}
 
 // print("predEagle", predEagle);
 // print("sigma_predEagle", sigma_predEagle);
 
}

model {
  // PSS regression priors
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | log(curr_predPSS),sigma_predPSS);
  
  // // Update likelihood for the posterior prediction
  // if(mean_propEagle > 0){
  // target += normal_lpdf(ln_RunSize | log(curr_predEagle), sigma_predEagle);}

}

generated quantities{
  //
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  // real ln_post_curr_predEagle;
  // 
  // real post_curr_predEagle;
  
  // Variable for posterior abundance
  real RunSize;
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);
  
  // Expenential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);
  
  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  // if(mean_propEagle > 0){
  // ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  // 
  // post_curr_predEagle = exp(ln_post_curr_predEagle);}
  
  //
  RunSize = exp(ln_RunSize);
  
}

