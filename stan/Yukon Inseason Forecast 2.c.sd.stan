// Yukon King Inseason Forecast Version 2.0
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
  
  // Pointer vector
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  
  // SD for daily passage estimates
  // matrix<lower=0> [n_dayPSS, n_yearPSS]SD_mat;
  
  vector <lower=0> [n_yearPSS]PSS_year_sd;
  
  real <lower=0> curr_PSS_year_sd;

  
  
  
}

// 
parameters {
   
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  real ln_RunSize;
  
  vector <lower=0> [n_yearPSS]cumPSS_est_vect;
  
  real <lower=0> cumPSS_curr_est;
  
}

transformed parameters{
  
   //Empty vector for cum sum of counts for each year to myyear -1
  vector<lower = 0> [n_yearPSS]cumHistPSS;
  
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
  
  // Loop to get cumulative counts
 for (i in 1:n_yearPSS){
   cumHistPSS[i] = sum(PSS_mat[,i]);
 }
 
 //  Calculate predPSS from aplha, beta and cumHistPSS for model section
 for (i in 1:n_yearPSS){
 ln_predPSS[i] = alpha + beta * cumPSS_est_vect[i];}
 
   // Empirical SD for update  
 sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS]) - ln_predPSS[loc_pf_years_PSS]);
 
 // Bring it back to reality
 RunSize = exp(ln_RunSize);
 
 // Sum curr_PSS
 cum_current_PSS = sum(curr_PSS);
 
 //
 ln_curr_predPSS = alpha + beta * cumPSS_curr_est ;
 
 // print("ln_predPSS",ln_predPSS);
  // print("cumPSS_est_vect:",cumPSS_est_vect);

}

model {
  //prior
  alpha ~ normal(0,1e6);
  beta ~ normal(0,10);
  sigma ~ normal(0,5);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf-(Pf_sigma^2/2), Pf_sigma);
  
  // // //Likelihood Function
  // for(y in 1:n_yearPSS){
  //   for(d in 1:n_dayPSS)
  //   
  //   // if(PSS_mat[d,y]>0){
  //     
  //   PSS_est_mat[d,y] ~ normal(PSS_mat[d,y],SD_mat[d,y]);
  //   
  // }
  
  
  for(y in 1:n_yearPSS){
    if(PSS_year_sd[y] > 0){
    cumPSS_est_vect[y] ~ normal(cumHistPSS[y], PSS_year_sd[y]);
  }}
  
  if(curr_PSS_year_sd >0)
   cumPSS_curr_est ~normal(cum_current_PSS,curr_PSS_year_sd);
  
  
    //Likelihood Function
  for(i in 1:n_yearPSS){
    // log(totalEOS[i]) ~ normal(log(predPSS), sigma);
    // target += normal_lpdf(totalEOS|predPSS,sigma);
    // totalEOS[i]~normal(predPSS[i], sigma);
    
    log(totalEOS[i])~normal(ln_predPSS[i],sigma);
    // log(totalEOS[i])~normal(log(predPSS[i]),sigma);
  }
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);

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
  
  ln_prior_pf = normal_rng(Pf - Pf_sigma^2/2, Pf_sigma);
  
  prior_pf = exp(ln_prior_pf);
  //
  predPSS = exp(ln_predPSS);
  //
  curr_predPSS = exp(ln_curr_predPSS);
  
  ln_post_curr_predPSS = normal_rng(ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
}









