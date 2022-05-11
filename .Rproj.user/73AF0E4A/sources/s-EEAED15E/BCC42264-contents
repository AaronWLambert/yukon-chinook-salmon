// Yukon King Inseason Forecast Version 4.0
//Notes: CDF version
//       Average across days GSI proportions with logit transformation
//       and beta liklihood for propCAN


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
  
  // Data for CDF
  matrix <lower=0> [n_dayPSS, n_yearPSS]PSS_cum_hist_mat;
  real <lower=0> cum_curr_PSS[n_dayPSS];
  
}

// 
parameters {
   
  // // Regression Parameters
  // real <lower=0> alpha;
  // real <lower=0> beta;
  real <lower=0> sigma;
  
  // CDF parameters
  real midPoint;
  real shape;
  
  // GSI parameters
  real <lower=0> phi;
  // real <lower=0, upper=1> propCAN_1;
  // real<lower=0, upper=1> propCAN[n_dayPSS];
  
    // First days proportion
  real propCAN_logit[n_dayPSS];
  
  
  //
  real ln_RunSize;
  
}

transformed parameters{
  
   //Empty vector for cum sum of counts for each year to myyear -1
  matrix [n_dayPSS, n_yearPSS]cumHistPSS;
  
  matrix [n_dayPSS, n_yearPSS]cummulative_histPSS;
  
  // Empty vector for predicted PSS counts in transformed parameters
  matrix [n_yearPSS, n_yearPSS]predPSS;
  
  // Predicted Runsize
  real RunSize;
  
  // Cummualative PSS for the year of interest
  real cum_current_PSS;
  
  // current prediction for PSS
  real curr_predPSS;
  
  // Parameters for GSI beta likelihood
  real paramA[n_dayPSS];
  
  real paramB[n_dayPSS];
  
  real delta[n_dayPSS];
  
  matrix[n_dayPSS, n_yearPSS]changeDelta;
  
  // // First days proportion
  // real propCAN_logit[n_dayPSS];
  vector<lower=0, upper=1> [n_dayPSS]propCAN;
  
  // real <lower=0> cumCurrentPSS[n_dayPSS];
   
   //  // Calculate the A an B parameter for beta dist
for(x in 1:n_dayPSS){
 paramA[x] = ((1-meanpropCAN[x])/(sd_meanpropCAN[x]^2)-(1/meanpropCAN[x]))*meanpropCAN[x]^2;

 paramB[x] = paramA[x] * ((1/meanpropCAN[x])-1);}
 
 
 //  Bring propCAN_logit into beta space
for(t in 1:n_dayPSS){
 propCAN[t] = exp(propCAN_logit[t])/(1+exp(propCAN_logit[t]));
}
 
  // Loop to get cumulative counts
 for (i in 1:n_yearPSS){
   for(j in 1:n_dayPSS){
     cumHistPSS[j,i] = 0;
   cumHistPSS[j,i] += (PSS_mat[j,i] * propCAN[j]);
 }
 }
 for(y in 1:n_yearPSS){
 cummulative_histPSS[,y] = cumulative_sum(cumHistPSS[,y]);}
 
 //
 
 
   for(d in 1:n_dayPSS){
 delta[d] = 1/(1+exp(-(dayPSS[d] - midPoint)/shape));}
 
 for(i in 1: n_yearPSS){
 changeDelta[1,i] = 0;}
 
 for(y in 2:n_yearPSS){
   for(d in 2:n_dayPSS)
 changeDelta[d,y] = delta[d]-delta[d-1];}
 //

 
 for(y in 1:n_yearPSS){
   for(d in 1:n_dayPSS){
 
 
 predPSS[d,y] = PSS_cum_hist_mat[d,y]/delta[d];}}
 
 
 
 // //  Calculate predPSS from aplha, beta and cumHistPSS for model section
 // for (i in 1:n_yearPSS){
 // ln_predPSS[i] = alpha + beta * cumHistPSS[i];}
 
 // // Bring it back to reality
 // RunSize = exp(ln_RunSize);
 
 // Sum curr_PSS after ajusting by propCAN
 
 cum_current_PSS = 0;
 for(j in 1:n_dayPSS){
 cum_current_PSS += (curr_PSS[j]*propCAN[j]);
 }
 //

 
 // ln_curr_predPSS = alpha + beta * cum_current_PSS;
 
  curr_predPSS = cum_current_PSS/delta[n_dayPSS];


  // print("ln_predPSS: ",ln_predPSS[]);
  // print("sigma:", sigma);
  // print("alpha:",alpha);
  // print("beta:",beta);
  // print("cumhist:",cumHistPSS);
  // print("propCAN:", propCAN);
  // print("cum_current_PSS:", cum_current_PSS);
  // print("ln_curr_predPSS:", ln_curr_predPSS);
  // print("logitCAn", propCAN_logit);
  // print("pss_mat",PSS_mat );

}

model {
  //prior
  // alpha ~ normal(0,1e6);
  // beta ~ normal(0,10);
  midPoint ~ normal(28,2);
  shape ~ normal(5,2);
  sigma ~ normal(0,5);
  phi ~ normal(0,1);
  propCAN_logit[1] ~ normal(0,1);
  
  
  // for(d in 2:n_dayPSS){
  // propCAN[d] ~ beta((((1-propCAN[(d-1)])/phi)-(1/propCAN[d-1]))*propCAN[d-1]^2,
  // ((((1-propCAN[(d-1)])/phi)-(1/propCAN[d-1]))*propCAN[d-1]^2)*((1/propCAN[d-1])-1));
// }
   
  for(d in 2:n_dayPSS){
   propCAN_logit[d] ~ normal(propCAN_logit[d-1],phi);
  }


  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  // //Likelihood Function
  // for(i in 1:n_yearPSS){
  //   // log(totalEOS[i]) ~ normal(log(predPSS), sigma);
  //   // target += normal_lpdf(totalEOS|predPSS,sigma);
  //   // totalEOS[i]~normal(predPSS[i], sigma);
  //   // Change here
  //   log(totalEOS[i])~normal((predPSS[i]),sigma);
  //   // log(totalEOS[i])~normal(log(predPSS[i]),sigma);
  // }
  
  // CDF Likelihood
  for(y in 1:n_yearPSS){
    for(d in 1:n_dayPSS){
  
    log(cumHistPSS[d,y])~normal(log(changeDelta[d,y]*predPSS[d,y]),sigma);}}
  
  

  //Canadian Proportion likelihood
  // propCAN ~ beta(((1 - meanpropCAN)/(sd_meanpropCAN) - (1/meanpropCAN)) * meanpropCAN^2, ((1/meanpropCAN)-1));
  for(i in 1:n_dayPSS){
  propCAN[i] ~ beta(paramA[i], paramB[i]);}
  // target += beta_lpdf(propCAN|paramA, paramB);
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | curr_predPSS,sigma);

}

// generated quantities{
//   //
//   real ln_prior_pf;
//   
//   real prior_pf;
//   // 
//   matrix [n_dayPSS,n_yearPSS]predPSS;
//   //
//   real curr_predPSS;
//   
//   real ln_post_curr_predPSS;
//   
//   real post_curr_predPSS;
//   
//   ln_prior_pf = normal_rng(Pf,Pf_sigma);
//   
//   prior_pf = exp(ln_prior_pf);
//   //
//   predPSS = exp(ln_predPSS);
//   //
//   curr_predPSS = exp(ln_curr_predPSS);
//   
//   ln_post_curr_predPSS = normal_rng(ln_curr_predPSS, sigma);
//   
//   post_curr_predPSS = exp(ln_post_curr_predPSS);
// }










