# Stan file descriptions
## Stan Scripts
|Script|Description|
|---|---|
|Yukon Inseason Forecast 1.0|This model is the base model using PSS cummulative passage regressed against Candian reconstructed EOS border passage.|
|Yukon Inseason Forecast 2.0|Extension of base model where Candian origin EOS passage is estimated in log space.|
|Yukon Inseason Forecast 2.1|Lognormal correction is incorporated . The correction is *subtracted* in this version.|
|Yukon Inseason Forecast 2.2|Lognormal correction is incorporated. The correction is *added* in this version.|
|Yukon Inseason Forecast 3.0|PSS passage adjusted by mean GSI proportions across days *before* to going into STAN model.|
|Yukon Inseason Forecast 3.1|GIS random walk with Beta distribtution. *(Need heavy restriction on phi_var draw to avoid negative shape parameter)*|
|Yukon Inseason Forecast 3.2|GSI random walk done with logit propcan from normal dist and converted back to proportion for beta likelihood. With reparamterized mean and sd GSI. *Likelihood from mean gsi proportion across that day of all years.*|
|Yukon Inseason Forecast 3.3|GSI random walk done with logit propcan from normal dist. propCAN likelihood from normal dist with meanpropCAN and Sd_meanPropCAN. *Likelihood from mean gsi proportion across that day of all years.*|
|Yukon Inseason Forecast 3.4|GSI random walk done with logit propcan from normal dist and converted back to proportion for beta likelihood. With reparamterized mean and sd GSI. *Likelihood from mean GSI across mean strata dates.*|
|Yukon Inseason Forecast 3.5|GSI random walk done with logit propcan from normal dist. propCAN likelihood from normal dist with meanpropCAN and Sd_meanPropCAN.*Likelihood from mean GSI across mean strata date.*|

## Note
This project is in progress and all results are preliminary.
