// Yukon King Inseason Forecast
//Notes:


//Components:
//  1) Pre-season forecast  (Prior)
//  2) PSS current counts
//  3) Regression for EOS Canadian
//  4) Logsistic curve for EOS Canadian
//  5) Steepness of curve (s) perhaps mean and sd of past years? As prior?
//
//  



data {
  
  // This years Preseason forecast
  real <lower=0>Pf;
  real <lower=0>Pf_sigma;
  

  // EOS Canadian counts by year up to myyear -1
  int <lower=0>n_totalEOS ;
  real<lower=0> totalEOS[n_totalEOS];
  
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
  real predPSS[n_yearPSS];
  
  // Variable for posterior abundance
  real RunSize;
  
  // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  //
  real curr_predPSS;
  
  // Loop to get cumulative counts
 for (i in 1:n_yearPSS){
   cumHistPSS[i] = sum(PSS_mat[,i]);
 }
 
 //  Calculate predPSS from aplha, beta and cumHistPSS for model section
 for (i in 1:n_yearPSS){
 predPSS[i] = alpha + beta * cumHistPSS[i];}
 
 //
 RunSize = exp(ln_RunSize);
 
 // Sum curr_PSS
 cum_current_PSS = sum(curr_PSS);
 
 //
 curr_predPSS = alpha + beta * cum_current_PSS;
}

model {
  //prior
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma ~ normal(0,10);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  //Likelihood Function
  for(i in 1:n_yearPSS){
    // log(totalEOS[i]) ~ normal(log(predPSS), sigma);
    // target += normal_lpdf(totalEOS|predPSS,sigma);
    // totalEOS[i]~normal(predPSS[i], sigma);
    // Change here
    log(totalEOS[i])~normal(log(predPSS[i]),sigma);
    // log(totalEOS[i])~normal(log(predPSS[i]),sigma);
  }
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | log(curr_predPSS),sigma);

}

generated quantities{
  //
  real ln_prior_pf;
  real prior_pf;
  
  real ln_post_curr_predPSS;
  real post_curr_predPSS;
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);
  // Expenential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);
  
  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma);
  post_curr_predPSS = exp(ln_post_curr_predPSS);
}



