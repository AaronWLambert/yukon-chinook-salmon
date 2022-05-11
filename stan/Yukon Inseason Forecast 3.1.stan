// Yukon King Inseason Forecast Version 3.0
//Notes: log(EOS)~cumPSS
//       Average across days GSI proportions 


//Components:
//  1) Pre-season forecast  (Prior)
//  2) PSS current counts
//  3) Regression for EOS Canadian
//
//  



data {
  
  // This years Preseason forecast
  real<lower=0> Pf;
  real<lower=0> Pf_sigma;
  

  // EOS Canadian counts by year up to myyear -1
  int<lower=0> n_totalEOS ;
  real<lower=0> totalEOS[n_totalEOS];
  
  // PSS days up to myDay
  int n_dayPSS ;
  int<lower=0> dayPSS[n_dayPSS];
  
  // Historic PSS years up to myYear -1
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // Matrix of historic counts by days & years
  int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // GSI proportions
  real <lower=0> meanpropCAN[n_dayPSS];
  real <lower=0> sd_meanpropCAN[n_dayPSS];
  
}

// 
parameters {
   
  // Regression Parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  // GSI parameters

  // real <lower=0, upper=1> propCAN_1;
  real<lower=.2, upper=1> propCAN[n_dayPSS];
  
    real <lower=0, upper=0.1> phi;
  //
  real ln_RunSize;
  
}

transformed parameters{
  
   //Empty vector for cum sum of counts for each year to myyear -1
  real cumHistPSS[n_yearPSS];
  
  // Empty vector for predicted PSS counts in transformed parameters
  real ln_predPSS[n_yearPSS];
  
  // Predicted Runsize
  real RunSize;
  
  // Cummualative PSS for the year of interest
  real cum_current_PSS;
  
  // log current prediction for PSS
  real ln_curr_predPSS;
  
  // // draw for GSI proportion
  // real<lower=0, upper=1> propCAN[n_dayPSS];
  
  // Parameters for GSI beta likelihood
  real paramA[n_dayPSS];
  
  real paramB[n_dayPSS];
  
  // First days proportion
  // propCAN[1] = propCAN_1;
  
  // Loop to get cumulative counts
 for (i in 1:n_yearPSS){
   cumHistPSS[i] = 0;
   
   for(j in 1:n_dayPSS){
   cumHistPSS[i] += (PSS_mat[j,i] * propCAN[j]);
 }
 }
 
 //  Calculate predPSS from aplha, beta and cumHistPSS for model section
 for (i in 1:n_yearPSS){
 ln_predPSS[i] = alpha + beta * cumHistPSS[i];}
 
 // Bring it back to reality
 RunSize = exp(ln_RunSize);
 
 // Sum curr_PSS
 cum_current_PSS = sum(curr_PSS);
 
 //
 ln_curr_predPSS = alpha + beta * cum_current_PSS;
 
//  // Calculate the A an B parameter for beta dist
for(x in 1:n_dayPSS){
 paramA[x] = ((1-meanpropCAN[x])/(sd_meanpropCAN[x]^2)-(1/meanpropCAN[x]))*meanpropCAN[x]^2;

 paramB[x] = paramA[x] * ((1/meanpropCAN[x])-1);}
//  
  // print("cumhist:",cumHistPSS);
  // print("propCAN:", propCAN);
}

model {
  //prior
  alpha ~ normal(0,1e6);
  beta ~ normal(0,10);
  sigma ~ normal(0,5);
  // phi ~ normal(0,1);
  phi~beta(1,5);
  // propCAN[1] ~ uniform(.5,1);
  propCAN[1]~ beta(2,2)T[0.2,1];
  
  for(d in 2:n_dayPSS){
  propCAN[d] ~ beta((((1-propCAN[(d-1)])/phi)-(1/propCAN[d-1]))*propCAN[d-1]^2,
  ((((1-propCAN[(d-1)])/phi)-(1/propCAN[d-1]))*propCAN[d-1]^2)*((1/propCAN[d-1])-1));}
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  //Likelihood Function
  for(i in 1:n_yearPSS){
    // log(totalEOS[i]) ~ normal(log(predPSS), sigma);
    // target += normal_lpdf(totalEOS|predPSS,sigma);
    // totalEOS[i]~normal(predPSS[i], sigma);
    // Change here
    log(totalEOS[i])~normal((ln_predPSS[i]),sigma);
    // log(totalEOS[i])~normal(log(predPSS[i]),sigma);
  }
  
  //Canadian Proportion likelihood
  // propCAN ~ beta(((1 - meanpropCAN)/(sd_meanpropCAN) - (1/meanpropCAN)) * meanpropCAN^2, ((1/meanpropCAN)-1));
  propCAN ~ beta(paramA, paramB);
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | ln_curr_predPSS,sigma);

}

generated quantities{
  //
  real ln_prior_pf;
  
  real prior_pf;
  // 
  real predPSS[n_yearPSS];
  //
  real curr_predPSS;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  ln_prior_pf = normal_rng(Pf,Pf_sigma);
  
  prior_pf = exp(ln_prior_pf);
  //
  predPSS = exp(ln_predPSS);
  //
  curr_predPSS = exp(ln_curr_predPSS);
  
  ln_post_curr_predPSS = normal_rng(ln_curr_predPSS, sigma);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
}






