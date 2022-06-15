
# R Script Subfolder Descriptions
## Bayesian Inseason Projection For Canadian-origin Yukon Chinook Salmon 
## Folders
|File|Description|
|---|---|
|Can vs PSS Trial.R| Script to explore the relationship between PSS Chinook passage and reconstructed end-of-season Canadian origin Chinook salmon. R^2 values were compared for different GSI proportion adjustments (i.e. naive vs. precise) and unadjusted PSS counts.|
|Dataframe preprocess| Script to import PSS passage information downloaded from ADFG excel sheet into R and pivot to long format for model use. When using the model inseason, this is the first scritp that will be called. Future iterations may have this as a function.|
|Retro Function| Script with function to calculate retrospective testing stats such as MAPE, PE, and RMSE.|
|Yukon Inseason Forecast| Script that all other scripts are built from. This is the base model which contains all the code found in functions and loops for running different versions of the model up to version 3.| 
|Retrospective Testing With Plots| Script that runs the model for historical years for the purpose of retorospective testing. Model outputs are saved/loaded and then the retorspective function is called and plots are generated from the retorpective function outputs.|
|Streamlined Function with Plot Script| Simplified script to run a single years projection and generate density plots. This script will be updated over time to include other figures and metrics.|
|Works in Progress| Subfolder containing scratch work and experimental scripts that may or may not be functional. For my refernce only|


## Note
This project is in progress and all results are preliminary.
