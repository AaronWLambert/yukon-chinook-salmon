// The model to be estimated. We model the output// Yukon King Inseason Forecast Version 2.c.1.1
//Notes: 


//Components:
//  1) Pre-season forecast  (Prior)
//  2) log(EOS-Can) ~ PSS
//  3) GSI incorporated (Beta likelihood)
//  4) Eagle sonar regression



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
  
  // GSI proportions
  real <lower=0> meanpropCAN[n_dayPSS];
  real <lower=0> sd_meanpropCAN[n_dayPSS];
  
  // Pointer vector
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  int <lower=0> loc_pf_years_PSS[n_yearsPF];  
  
  // GSI likelihod beta parameters
  // real <lower=0> paramA[n_dayPSS];
  
  // real <lower=0> paramB[n_dayPSS];
  
  
  
}

transformed data{
  
    // Parameters for GSI beta likelihood
  real <lower=0> paramA[n_dayPSS];
  // 
  real <lower=0> paramB[n_dayPSS];
  
  
  
    // Calculate the A an B parameter for beta dist
  for(x in 1:n_dayPSS){

      paramA[x] = ((1-meanpropCAN[x])/(sd_meanpropCAN[x]^2)-(1/meanpropCAN[x]))*meanpropCAN[x]^2;

      paramB[x] = paramA[x] * ((1/meanpropCAN[x])-1);}
  
  
  
}

// 
parameters {
   
  // Regression Parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  // Eagle regression paramerters
  real <lower=0> alpha_eagle;
  real <lower=0> beta_eagle;
  real <lower=0> sigma_eagle;
  
  // GSI parameters
  real ln_phi;
  
  // Logit canadian proportion draws
  real propCAN_logit[n_dayPSS];
  
  
  // Log response variable
  real ln_RunSize;
  
}

transformed parameters{
  
  // Empty vector for cum sum of counts for each year to myyear -1
  real <lower=0> cumHistPSS[n_yearPSS];
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]ln_predPSS;
  
  //Empty vector for cum sum of counts for historic Eagle passage
  real cumHistEagle[n_yearEagle];
  
  // Log prediction of EOS Eagle passage
  vector [n_yearEagle]ln_predEagle;
  
  // Predicted Runsize
  real <lower=0> RunSize;
  
  // Cummualative PSS for the year of interest
  real <lower=0> cum_current_PSS;
  
  // log current prediction for PSS
  real ln_curr_predPSS;
  
  // Cummulative current Eagle passage
  real cum_current_Eagle;
  
  // log current eagle prediction
  real ln_curr_predEagle;
  
  // PSS likelihood sd
  real sigma_predPSS;
  
  // Eagle liklihood sd
  real sigma_predEagle;
  
  // First days proportion
  vector <lower=0, upper=1> [n_dayPSS]propCAN;
  
  real <lower=0> phi;
  
  phi = exp(ln_phi);
   
  // Bring propCAN_logit into beta space
  for(t in 1:n_dayPSS){
     
      propCAN[t] = exp(propCAN_logit[t])/(1+exp(propCAN_logit[t]));
      }

  // Loop to get cumulative counts
  for (i in 1:n_yearPSS){
   
   cumHistPSS[i] = 0;
   
    for(j in 1:n_dayPSS){
     
        cumHistPSS[i] += (PSS_mat[j,i] * propCAN[j]);
       } // j loop
    } // t loop
 
 
  // Calculate predPSS from aplha, beta and cumHistPSS for model section
  for (i in 1:n_yearPSS){
    ln_predPSS[i] = alpha + beta * cumHistPSS[i];}
 
 
 // Empiircal SD for updating
 sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS]) - ln_predPSS[loc_pf_years_PSS]);
 
 
 // Bring it back to reality
 RunSize = exp(ln_RunSize);
 
 // Sum curr_PSS after ajusting by propCAN
 
   cum_current_PSS = 0;
 
 for(j in 1:n_dayPSS){
   
  cum_current_PSS += (curr_PSS[j]*propCAN[j]);

    }
 
 // PSS current predictions
 ln_curr_predPSS = alpha + beta * cum_current_PSS;
 
  // Eagle Sonar calculations ////////////////////////////////
 
   // Loop to get cumulative Eagle counts
  for (i in 1:n_yearEagle){
   cumHistEagle[i] = sum(Eagle_mat[,i]);
   }
 
  // Summ current year Eagle passage
  cum_current_Eagle = sum(curr_Eagle);

  for (i in 1:n_yearEagle){
  // Eagle regression w/ slope
   ln_predEagle[i] = alpha_eagle + beta_eagle * cumHistEagle[i];}

  // Empirical SD for update  
  sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS]) - ln_predEagle[loc_pf_years_Eagle]);
 
  //Eagle prediction with slope
  ln_curr_predEagle = alpha_eagle + beta_eagle * cum_current_Eagle;
 
 
 
 
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
  
  // Regression priors
  alpha ~ normal(0,1e6);
  beta ~ normal(0,1);
  sigma ~ normal(0,5);
  
  // Eagle regression priors
  alpha_eagle ~ normal(0,1e6);
  beta_eagle ~ normal(0,10);
  sigma_eagle ~ normal(0,5);
  
  // Random walk GSI priors
  ln_phi ~ normal(-1,1);
  propCAN_logit[1] ~ normal(0,1);// change back to (0,1)
  
  
  // Random walk 
  propCAN_logit[2:n_dayPSS] ~normal(propCAN_logit[1:n_dayPSS-1], phi);

  // Preseason Forcast
  // ln_RunSize ~ normal(Pf - Pf_sigma^2/2, Pf_sigma);
  ln_RunSize ~ normal(Pf - Pf_sigma^2/2, Pf_sigma);
  
  // PSS regression likelihood
  log(totalEOS)~normal(ln_predPSS, sigma);
  
  
  // Eagle regression likelihood
  for(i in 1:n_yearEagle){
    
    log(totalEOS[loc_eagle_years][i])~normal((ln_predEagle[i]), sigma_eagle);
  }

  //Canadian Proportion likelihood
  propCAN ~ beta(paramA, paramB);
     // propCAN ~ beta(2,2);
  // Update with PSS likelihood for the posterior prediction
  // target += normal_lpdf(ln_RunSize | ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);
  target += normal_lpdf(ln_RunSize | ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);
  
  // Update with Eagle likelihood for the posterior prediction
  target += normal_lpdf(ln_RunSize | ln_curr_predEagle - sigma_predEagle^2/2, sigma_predEagle);
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
  
  vector [n_yearEagle]predEagle;
  
  real curr_predEagle;
  
  real ln_post_curr_predEagle;
  
  real post_curr_predEagle;
  
  // ln_prior_pf = normal_rng(Pf - Pf_sigma^2/2, Pf_sigma);
  ln_prior_pf = normal_rng(Pf - Pf_sigma^2/2, Pf_sigma);
  
  prior_pf = exp(ln_prior_pf);
  //
  
  predPSS = exp(ln_predPSS);
  
  //
  curr_predPSS = exp(ln_curr_predPSS);
  
  // ln_post_curr_predPSS = normal_rng(ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);
  ln_post_curr_predPSS = normal_rng(ln_curr_predPSS - sigma_predPSS^2/2, sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
      
  predEagle = exp(ln_predEagle);
  //
  curr_predEagle = exp(ln_curr_predEagle);
  
  ln_post_curr_predEagle = normal_rng(ln_curr_predEagle - sigma_predEagle^2/2, sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);
}









