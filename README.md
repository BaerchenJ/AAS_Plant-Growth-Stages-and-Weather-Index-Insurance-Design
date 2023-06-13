# AAS_Plant-Growth-Stages-and-Weather-Index-Insurance-Design
The dataset shared by Jing Zou, Martin Odening, and Ostap Okhrin comprises four folders:

1. PSANOVA&MGCV: This folder contains the data table "WeatherIndices&Yields.xlsx," which encompasses weather indices and yield-related data for each county in the sample. It also includes R scripts, namely "wc-PSANOVA.r" and "wc-MGCV.r," used for running whole-cycle models employing two methods. Additionally, the folder contains R codes, "ps-PSANOVA.r" and "ps-MGCV.r," utilized for running phase-division models using the same methods. The response variable can be modified to toggle between quadratic utility and exponential utility. By adjusting the "alpha" parameter, models with three levels of risk aversion can be compared. Caution is advised when calculating the scale-dependent Root Mean Square Error (RMSE) of models under exponential utility, as it requires a two-step process. The folder also includes visualizations that portray the relationship between weather and yield loss for both whole-cycle and phase-division models.

2. CrossValidation: This folder contains separate excel tables for the training (16 years), validation (4 years), and test (3 years) datasets. The data is organized chronologically, considering that weather indices and yield losses exhibit time series characteristics.

3. HedgingEffectiveness: Within this folder, you will find the computation of three evaluation criteria: Expected Utility (EU) ratio, mean root square loss (MRSL), and variance of revenue (VAR). Additionally, it provides plots illustrating the EU ratio.

4. CalculateWeatherIndex&YieldLoss: For readers seeking to comprehend the calculation process of weather indices and yield losses from the raw data, this folder presents relevant information. It includes R codes for generating time series plots of yield, detrended yield, yield loss, growing degree days (GDD), and cumulative rainfall (CR).

Raw data sources:
County-level soybean yield data and phenology dates of 96 counties in Illinois from 1998 to 2020 are from the National Agricultural Statistics Services (NASS).
Meteorological data are from
- Thornton, M. M., Shrestha, R., Wei, Y., Thornton, P. E., Kao, S., & Wilson, B. E. (2020). Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 4. ORNL DAAC, Oak Ridge, Tennessee, USA.

The authors sincerely appreciate your interest in their study and value the time you have dedicated to reviewing and potentially replicating their research.

Sincerely,
Jing Zou, Martin Odening, and Ostap Okhrin
