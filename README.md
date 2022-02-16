# CapturingThePoolDilution_Effect_In_Group_Testing_Regression_A_Bayesian_Approach

This repository provides code which implements the Bayesian group testing regression methodology described in 'Capturing the Pool Dilution Effect in Group
Testing Regression: A Bayesian Approach' by Stella Self, Christopher McMahan and Stefani Mokalled.

Simulations_Individual_Testing implements the methodoldogy for simulated individual testing data under gamma biomarker distributons.

Simulations_Dorfman_Testing implements the methodology for simulated group testing data under the Dorfman retesting approach described in the manuscript with two rounds of model fitting under gamma biomarker distributions.

Simulations_Array_Testing implements the methodology for simulated group testing data under the square array retesting approach described in the manuscript with two rounds of model fitting under gamma biomarker distributions.

Data_Application is a general implementation for the method for use on real data under gamma biomarker distributions.
It requires two input files, one containing the individual level covariates and one containing the pool level observations, pool sizes, and the individual who contributed to each pool. 
It outputs estimates of the model parameters and a list of positive pools.
After these positive pools have been retested in any manner desired, the pool data file can be updated and the code applied again.
This process can be repeated as many times as desired.

Individual_Data provides an example of the individual data file necessary to run the Data_Application file.

Pool_Data provides an example of the pool data file necessary to run the Data_Application file.
