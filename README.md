# Sequential_Bayes_Epi
R code for implementing the methods relating to: https://arxiv.org/abs/2405.13537

Cite as: Whitaker, S. A., Golightly, A., Gillespie, C. S., and Kypraios, T. (2024). Sequential Bayesian inference for stochastic epidemic models of cumulative incidence. arXiv preprint  	arXiv:2405.13537.

The code in this repository can be used to reproduce the analysis for the applications in Section 4 of the paper.

## Abstract

Epidemics are inherently stochastic, and stochastic models provide an appropriate way to describe and analyse such phenomena. Given temporal incidence data consisting of, for example, the number of new infections or removals in a given time window, a continuous-time discrete-valued Markov process provides a natural description of the dynamics of each model component, typically taken to be the number of susceptible, exposed, infected or removed individuals. Fitting the SEIR model to time-course data is a challenging problem due incomplete observations and, consequently, the intractability of the observed data likelihood. Whilst sampling based inference schemes such as Markov chain Monte Carlo are routinely applied, their computational cost typically restricts analysis to data sets of no more than a few thousand infective cases. Instead, we develop a sequential inference scheme that makes use of a computationally cheap approximation of the most natural Markov process model. Crucially, the resulting model allows a tractable conditional parameter posterior which can be summarised in terms of a set of low dimensional statistics. This is used to rejuvenate parameter samples in conjunction with a novel bridge construct for propagating state trajectories conditional on the next observation of cumulative incidence. The resulting inference framework also allows for stochastic infection and reporting rates. We illustrate our approach using synthetic and real data applications.

## Organisation

Code for each application (4.1 Simulation Study, 4.2 Ebola and 4.3 Covid-19 in New York) can be found in the corresponding directories. To run the inference scheme start with "run.R".

Each directory includes the following files:
- run.R - Containing the code to run the inference scheme;
- data.csv - A file containing the data for the analysis;
- hazard.R - Containing functions to calculate the conditioned hazard;
- bridge.R - Containing functions to implement the bridge in parallel;
- particle_filter.R - Containing a function to implement Algorithm 2;
- forecast.R - Containing a function to generate forecasts from the output of run.r;
- forecast_run.R - Code to generate and store the forecasts.

## Timings

The following timings are provided as a rough guide for the user. Note that these timings are based on running the schemes on a PC with a 2.6 GHz clock speed, parallelised across 20 cores.

- 4.1 Simulation Study: 16 minutes
- 4.2 Ebola: 5 hours
- 4.3 Covid-19 in New York: 50 minutes