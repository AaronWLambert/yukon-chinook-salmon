// Yukon King Inseason Forecast Version 6.0
// Notes: Proportion estimator with run timing with scale estimated.....
//        Normal Likelihood for prop fitting


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
  vector <lower=0> [n_dayPSS]dayPSS;
  int <lower=0> n_dayPSS_all;
  vector <lower=0> [n_dayPSS_all]dayPSS_all;
  
  // Historic PSS years up to myYear -1
  int<lower=0> n_yearPSS;
  int<lower=0> yearPSS[n_yearPSS];
  
  // Matrix of historic counts by days & years
  // int<lower=0> n_hist_counts;
  matrix<lower=0> [n_dayPSS, n_yearPSS]PSS_mat;
  
  // Matrix of historic PSS passage for complete historic years
  matrix <lower=0> [n_dayPSS_all, n_yearPSS]PSS_mat_all;
  // matrix <lower=0> [n_dayPSS_all, n_yearPSS]PSS_mat_prop_all;
      // Matrix of cumulative historic PSS passage by days and years
  matrix <lower=0> [n_dayPSS_all, n_yearPSS]cum_PSS_mat_all;
  
  // Current year counts up to myDay
  int<lower=0> n_curr_PSS;
  real<lower=0> curr_PSS[n_curr_PSS];
  
  // Pointer vector
  int <lower=0> loc_pf_years_PSS[n_yearsPF];
  
  // April air mean temp
  vector [n_yearPSS]may_sst_hist;
  real  may_sst_curr;
  
  // Historical deviations from avg midpoint timing
  vector[n_yearPSS]dev_hist;
  // real avg_midpoint;
  
  real myDay;
  
  // real avg_mid;
  // 
  // real shape;
  
}

transformed data{
  
  // Vector to hold sum of counts for each historic year up to myDay
  vector [n_yearPSS]cumHistPSS;
  
  // Vector for the sum of each historic year to the EOS
  real cumHistPSS_all[n_yearPSS];
  
  // Sum of the current year up to myDay
  real cum_current_PSS;
  
  // matrix [n_dayPSS_all, n_yearPSS]PSS_mat_prop_logit;
  
  // Loop to get cumulative counts
    for (i in 1:n_yearPSS){
        cumHistPSS_all[i] = sum(PSS_mat_all[,i]);
      }
  
  // Loop to get cumulative historic counts
    for (i in 1:n_yearPSS){
        cumHistPSS[i] = sum(PSS_mat[,i]);
      }
  // Sum up the current year PSS passage to myDay
  cum_current_PSS = sum(curr_PSS);
  
  // Logit transform the proportions for likelihood comparison
  // for(y in 1:n_yearPSS){
  //   for(d in 1:n_dayPSS_all){
  //     if(PSS_mat_prop_all[d,y] == 0){
  //     PSS_mat_prop_logit[d,y] = log((PSS_mat_prop_all[d,y]+0.00001) ./ (1-(PSS_mat_prop_all[d,y] +0.00001)));}else{
  //       if(PSS_mat_prop_all[d,y] == 1){
  //         PSS_mat_prop_logit[d,y] = log((PSS_mat_prop_all[d,y]-0.00001) ./ (1-(PSS_mat_prop_all[d,y] - 0.00001)));
  //       }else{
  //         PSS_mat_prop_logit[d,y] = log((PSS_mat_prop_all[d,y]) ./ (1-(PSS_mat_prop_all[d,y])));
          
  //       }
  //     }
  //   }
  // }
  
}

// 
parameters {
   
  // Regression slope param
  real <lower=0> alpha;
  real <lower=0> beta;
  real <lower=0> sigma;
  
  // Midpoint regression param
  real alpha_t;
  real beta_sst;
  real <lower=0> sigma_t;
  
  // logistic curve paramaters for current year
  real <lower=0> scale;
  real <lower=0> mid;
  real <lower=0> shape;
  real <lower=0> sigma_logistic;
  
  real ln_RunSize;
  
}

transformed parameters{
  
  vector [n_yearPSS]predPSS;
  
  real curr_predPSS;
  
  real sigma_predPSS;
  
  vector [n_yearPSS]dev_predMP;
  
  real curr_dev_predMP;
  
  real expectedProp_curr;
  
  vector [n_yearPSS]expectedProp_hist;
  
  real day_adj_curr;
  
  vector [n_yearPSS]day_adj_hist;
  
  vector [n_yearPSS]est_predPSS;
  
  real curr_est_predPSS;
  
  vector [n_yearPSS]predPSS_hist;
  
  vector [n_dayPSS_all]ps_cumm_pred;
  
  // vector  [n_dayPSS_all]ps_cumm_pred_logit;
  
  // vector [n_dayPSS_all]paramA;
  // 
  // vector [n_dayPSS_all]paramB;

  // Mid point regression
  dev_predMP = alpha_t + beta_sst * may_sst_hist;
  
  // This years prediction
  curr_dev_predMP = alpha_t + beta_sst * may_sst_curr;
  
  // Fit logisitic curve to all observed years except myYear at one time
  for(d in 1:n_dayPSS_all){
    ps_cumm_pred[d] = scale * (1/(1 + exp((-(dayPSS_all[d] - mid)) / shape)));}
  
  // Adjust day by estimated deviation for current year and historic years
  day_adj_curr = myDay - curr_dev_predMP;
  
  day_adj_hist = myDay - dev_predMP;
  
  // Regress cumulative EOS PSS against the EOS Can run size
  for (i in 1:n_yearPSS){
    predPSS[i] = alpha + beta * cumHistPSS_all[i];}
  
  // Calculate expected proportions
  expectedProp_curr = 1/(1+exp((-(day_adj_curr - mid))/ shape));
  
  for(y in 1:n_yearPSS){
  expectedProp_hist[y] = 1/(1+exp((-(day_adj_hist[y] - mid))/ shape));}
 
  // Estimated historic years by proportions
  for(y in 1:n_yearPSS){
  est_predPSS[y] = cumHistPSS[y] /expectedProp_hist[y];}
  
  predPSS_hist = alpha + beta * est_predPSS;
  
  // Empirical SD for update
  sigma_predPSS = sd(log(totalEOS[loc_pf_years_PSS]) - predPSS_hist[loc_pf_years_PSS]);
  
  curr_est_predPSS = cum_current_PSS/expectedProp_curr;
 
  //
  curr_predPSS = alpha + beta * curr_est_predPSS;
 
 
  // 
  // // print("ln_predPSS",ln_predPSS);
  // print("ps_cumm_pred",  ps_cumm_pred);
  // print("Ps_cumm_pred_logit", ps_cumm_pred_logit)
  
}

model {
  
  // Midpoint regression priors
  alpha_t ~ normal(0,1e6);
  beta_sst ~ normal(0,100);
  sigma_t ~ normal(0,100);
  
  // Logistic curve fitting priors
  scale ~ normal(150000, 60000);
  mid ~ normal(175,5);
  shape ~ normal(5,2);
  sigma_logistic ~ normal(0,5);
  
  // Regression priors
  alpha ~ normal(0,1e6);
  beta ~ normal(0,10);
  sigma ~ normal(0,5);
  
  
  
  // Preseason Forcast
  ln_RunSize ~ normal(Pf, Pf_sigma);
  
  // Midpoint deviations regression likelihood
  dev_hist ~ normal(dev_predMP, sigma_t);
  
  // Logistic curve fitting likelihood
  for(d in 1:n_dayPSS_all){ 
    for(y in 1:n_yearPSS)
    // if(cum_PSS_mat_all[d,y] > 0){
      
    log(cum_PSS_mat_all[d,y]+1e-5)  ~ normal(log(ps_cumm_pred[d]+1e-5), sigma_logistic);}
    // }
  
  for(i in 1:n_yearPSS){
    log(totalEOS[i])~normal(log(predPSS[i]),sigma);
  }
  
  // Update for the posterior
  target += normal_lpdf(ln_RunSize | log(curr_predPSS), sigma_predPSS);

}

generated quantities{
  
  real ln_prior_pf;
  
  real prior_pf;
  
  real ln_post_curr_predPSS;
  
  real post_curr_predPSS;
  
  real post_curr_predMP;
  
  real RunSize;
  
  ln_prior_pf = normal_rng(Pf, Pf_sigma);
  
  prior_pf = exp(ln_prior_pf);
  //
  // predPSS = exp(ln_predPSS);
  //
  // curr_predPSS = exp(ln_curr_predPSS);
  
  ln_post_curr_predPSS = normal_rng(log(curr_predPSS), sigma_predPSS);
  
  post_curr_predPSS = exp(ln_post_curr_predPSS);
  
  post_curr_predMP = normal_rng(curr_dev_predMP, sigma_t);
  
   // Bring it back to reality
 RunSize = exp(ln_RunSize);
  
  
 }
