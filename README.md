# Bayes_SEM
R code for implementing the methods relating to: INSERT_LINK

The paper considers sequential Bayesian inference for S(E)IR model parameters and latent cumulative incidence process given noisy observations. A sequential inference scheme is developed, making use of a computationally cheap approximation of the most natural Markov process model with a novel bridge construct for propagating state trajectories conditional on the next observation. The developed scheme also allows for stochastic infection and reporting rates to be implemented.

The code in this repository can be used to reproduce the analysis for the applications in Section 4 of the paper.

## Organisation

Code for each application (4.1 Simulation Study, 4.2 Ebola application and 4.3 Covid application) can be found in the corresponding directories. To run the inference scheme start with "run.r".

Each directory also includes the following:
- hazard.r - Containing functions to calculate the conditioned hazard;
- bridge.r - Containing functions to implement the bridge in parallel;
- particle_filter.r - Containing a function to implement Algorithm 2;
- forecast.r - Containing a function to generate forecasts from the output of run.r;
- forecast_run.r - Code to generate and store the forecasts.