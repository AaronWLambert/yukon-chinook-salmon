// Yukon King Inseason Forecast version 1.0
//Notes:


//Components:
//  1) Pre-season forecast  (Prior)
//  2) EOS-Can ~ PSS in normal space
//  3) GSI incorporated (Beta likelihood)
//  4) Eagle Sonar regression added to likelihood
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
  
  // Location vector of integers
  int <lower=0> loc_eagle_years[n_yearEagle];
  int <lower=0> loc_pf_years_Eagle[n_yearsPF];
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  
  // Propotion of observed fish at Eagle
  real mean_propEagle;
  
}

transformed data{
  
  // Parameters for GSI beta likelihood
  real paramA[n_dayPSS];
  
  real paramB[n_dayPSS];
  
  // Calculate the A an B parameter for beta dist
  for(x in 1:n_dayPSS){
    
      paramA[x] = ((1-meanpropCAN[x])/(sd_meanpropCAN[x]^2)-(1/meanpropCAN[x]))*meanpropCAN[x]^2;

      paramB[x] = paramA[x] * ((1/meanpropCAN[x])-1);}
  
}

parameters {
   
  // PSS regression parameters
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  // Eagle regression parameters
  real <lower=0> alpha_eagle;
  real <lower=0> beta_eagle;
  real <lower=0> sigma_eagle;
  
  // GSI parameters
  real ln_phi;
  
  // Logit canadian proportion draws
  real propCAN_logit[n_dayPSS];
  
  // Log Canadian run size (Integrated run size estimate)
  real ln_RunSize;
  
}


transformed parameters{
  
   //Empty vector for cum sum of counts for each year to myyear -1
  real cumHistPSS[n_yearPSS];
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearPSS]predPSS;
  
  //Empty vector for cum sum of counts for each year to myyear -1
  vector [n_yearEagle]cumHistEagle;
  
  // Empty vector for predicted PSS counts in transformed parameters
  vector [n_yearEagle]predEagle;
  
  // Variable for cummulative myYear PSS passage
  real cum_current_PSS;
  
  //
  real curr_predPSS;
  
  //
  real cum_current_Eagle;
  
  //
  real curr_predEagle;
  
  real sigma_predPSS;
  
  real sigma_predEagle;
  
  // First days proportion
  vector<lower=0, upper=1> [n_dayPSS]propCAN;
 
  real phi;
 
  phi = exp(ln_phi);
  
  //  Bring propCAN_logit into beta space
   for(t in 1:n_dayPSS){
     
      propCAN[t] = exp(propCAN_logit[t])/(1+exp(propCAN_logit[t]));
      }

  // Loop to get cumulative counts
  for (i in 1:n_yearPSS){
   
   cumHistPSS[i] = 0;
   
    for(j in 1:n_dayPSS){
     
   cumHistPSS[i] += (PSS_mat[j,i] * propCAN[j]);
     }
    }
 
  // Loop to sum Eagle counts up to myDay
  for (i in 1:n_yearEagle){
     cumHistEagle[i] = sum(Eagle_mat[,i]);
    }
 
  // Summ current year Eagle passage
   cum_current_Eagle = sum(curr_Eagle);
   
  //  Calculate predPSS from aplha, beta and cumHistPSS for model section
  for (i in 1:n_yearPSS){
    predPSS[i] = alpha + beta * cumHistPSS[i];}
 
  // Sigma calculation for update
  //  Only years used in PF comparison are used here
  sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS] ./ predPSS[loc_pf_years_PSS]));
 
  // Sum curr_PSS after ajusting by propCAN
  cum_current_PSS = 0;
 
  for(j in 1:n_dayPSS){
   
    cum_current_PSS += (curr_PSS[j]*propCAN[j]);

  }
 
  //
    curr_predPSS = alpha + beta * cum_current_PSS;
 
  // Eagle Component using observed propotions
  if(mean_propEagle > 0){
  curr_predEagle = (cum_current_Eagle+1)/mean_propEagle;
  
  predEagle =  (cumHistEagle+1)/mean_propEagle;
 

  // Empirical SD for update  
  sigma_predEagle = sd(log(totalEOS[loc_pf_years_PSS] ./ predEagle[loc_pf_years_Eagle]));}
 
}

model {
  // PSS regression priors
  alpha ~ normal(0,1e15);
  beta ~ normal(0,1e15);
  sigma ~ normal(0,10);
  
  // Eagle regression priors
  alpha_eagle ~ normal(0,1e15);
  beta_eagle ~ normal(0,1e15);
  sigma_eagle ~ normal(0,10);
  
    // Random walk GSI priors
  ln_phi ~ normal(-1,1);
  propCAN_logit[1] ~ normal(0,1);
  
  // GSI random walk 
  propCAN_logit[2:n_dayPSS] ~normal(propCAN_logit[1:n_dayPSS-1], phi);
  
    // Canadian Proportion likelihood
  propCAN ~ beta(paramA, paramB);
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  // PSS regression likelihood 
  for(i in 1:n_yearPSS){
    
    log(totalEOS[i])~normal(log(predPSS[i]),sigma);
    
  }
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | log(curr_predPSS),sigma_predPSS);
  
  // Update likelihood for the posterior prediction
  if(mean_propEagle > 0){
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
  
  real RunSize;
  
  // Predicted forecast
  ln_prior_pf = normal_rng(Pf,Pf_sigma);
  
  // Expenential of log predicted preseason forecast
  prior_pf = exp(ln_prior_pf);
  
  // Predicted PSS passage
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  if(mean_propEagle > 0){
  ln_post_curr_predEagle = normal_rng(log(curr_predEagle), sigma_predEagle);
  
  post_curr_predEagle = exp(ln_post_curr_predEagle);}
  
  RunSize = exp(ln_RunSize);
}



