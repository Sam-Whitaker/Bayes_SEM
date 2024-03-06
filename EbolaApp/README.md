# Ebola application
Here you will find the code for replicating the Ebola application.

ebolaData.csv - Weekly incidence counts of Ebola in Seirra Leone spanning  

ebolaInf.r - Functions for inference for Ebola application via particle filter

ebolaForecast.r - Function to perform forecasting for Ebola application using particle filter output

ebolaRun.r - Code to run particle filter and produce forecasts

Summaries (mean and 95% credible interval) for model parameters through time:

![plot](~/summariesEbola.pdf)

Five one-step-ahead forecasts for last five (non-zero) observations (using particle filter output)

![plot2](~/forecastsEbola.pdf)