# AAS_Plant-Growth-Stages-and-Weather-Index-Insurance-Design
Data and codes for replication shared by Jing Zou, Martin Odening, and Ostap Okhrin.

Hello my dear readers:

There are four folders in this dataset.

#########################################################################################################
If you would like to replicate the core result table of the paper, please open this folder

### PSANOVA&MGCV

Inside this folder, the data table is "WeatherIndices&Yields.xlsx", 
which contains all the weather indices and yield-related data of each county in our sample.
"wc-PSANOVA.r" and "wc-MGCV.r" are R scripts for running whole-cycle models using two methods.
"ps-PSANOVA.r" and "ps-MGCV.r" are R codes for running phase-division models using two methods.
You could shift between quadratic utility and exponential utility by changing the response variable.
You could also play with "alpha" to compare the models with three levels of risk aversion.
Cautions:
When you calculate the scale-dependent Root Mean Square Error (RMSE) of models under exponential utility, 
please be careful because converting fitted values back needs two steps! 

Yes! Do not forget the visualization of whole-cycle and phase-division weather-yield loss relation!


##########################################################################################################
To tell you the truth, tuning the parameter, the number of knots, does not bring much performance gain!
However, to prove that we really did the cross validation procedure, we provide this folder

### CrossValidation
In this folder, training (16 years), validation (4 years) and test (3 years) datasets are separated according to the year order 
in three excel tables because weather indices and yield losses are time series.


##########################################################################################################
After replicating the result table, we move on to hedging effectiveness part in this folder

### HedgingEffectiveness

Here you could see clearly how three evaluation criteria Expected Utility (EU) ratio, mean root square loss (MRSL) 
and variance of revenue (VAR) are calculated. 

Yes! Do not forget the EU ratio plots!


##########################################################################################################
For readers who are more curious about how weather indices and yield losses are calculated from the raw data,
we encourage you to investigate this folder

### CalculateWeatherIndex&YieldLoss
Also, you could obtain the time series plots of yield, detrended yield, yield loss, GDD and CR with the R codes.


Thank you very much for reading our paper and even replicating the study!

