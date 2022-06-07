# yukon-chinook-salmon
## Bayesian Inseason Projection For Canadian-origin Yukon Chinook Salmon 
## Purpose
A management tool created to generate an inseason projection for each day during the summer Chinook salmon run on the Yukon River to estimate the end-of-season escapement of Canadian-origin Chinook salmon into Canda. This method utilizes a Bayesian updating approach using Rstan in R. Users can generate projections retrospectively for model testing and selecting purposes or to generate inseason projections for the current year.
## Folders
|File|Description|
|---|---|
|Data|Data files that are used in the projection model. These data are generally preprocessed in an R script and fed into the Stan model.|
|R|Scripts used to preprocess data and to call the Stan model using Rstan. The most general script is called Yukon Inseason Forecast, which can easily be used to generate predictions for a single day in any year of interest.|
|Figs|Plots and figures generated for analysis and visualization.|
|Output|Model generated results for retrospective testing. These files are not included in the repository as they are too large to upload.|
|Stan| Stan scripts that are called in the Rstan model. These include different versions of the model utilizing methods that are currently being assessed for their predictive performance.|

## Note
This project is in progress and all results are preliminary.
