# Regularised B-splines projected Gaussian Process priors to estimate time-trends of age-specific COVID-19 deaths related to vaccine roll-out

Mélodie Monod, Alexandra Blenkinsop, Andrea Brizzi, Yu Chen, Vidoushee Jogarah, Yuanrong Wang, Samir Bhatt and Oliver Ratmann (2021). Regularised B-splines projected Gaussian Process priors to estimate time-trends of age-specific COVID-19 deaths related to vaccine roll-out, http://arxiv.org/abs/2106.12360.

## Abstract
The COVID-19 pandemic has caused severe public health consequences in the United States.
In this study, we use a hierarchical Bayesian model to estimate the age-specific COVID-19 attributable deaths over time in the United States. The model is specified by a novel non-parametric spatial approach, a low-rank Gaussian Process (GP) projected by regularised B-splines. We show that this projection defines a new GP with attractive smoothness and computational efficiency properties, derive its kernel function, and discuss the penalty terms induced by the projected GP. Simulation analyses and benchmark results show that the spatial approach performs better than standard B-splines and Bayesian P-splines and equivalently well as a standard GP, for considerably lower runtimes. The B-splines projected GP priors that we develop are likely an appealing addition to the arsenal of Bayesian regularising priors.
We apply the model to weekly, age-stratified COVID-19 attributable deaths reported by the US Centers for Disease Control, which are subject to censoring and reporting biases. Using the B-splines projected GP, we can estimate longitudinal trends in COVID-19 associated deaths across the US by 1-year age bands. These estimates are instrumental to calculate age-specific mortality rates, describe variation in age-specific deaths across the US, and for fitting epidemic models. Here, we couple the model with age-specific vaccination rates to show that lower vaccination rates in younger adults aged 18-64 are associated with significantly stronger resurgences in COVID-19 deaths, especially in Florida and Texas. These results underscore the critical importance of medically able individuals of all ages to be vaccinated against COVID-19 in order to limit fatal outcomes. 

## About this repository
| Folder    | Description |
|-----------|------------------------------------------------------|
| misc   | Clean the incomplete age and state-specific mortality data reported by the CDC |
| inst | Reproduce the result of the "Real data analysis" section  |
| simulations      | Reproduce the result of the "Simulations" section |
| predictions | Reproduce the results of the "Real world benchmarks" |
| inst/results | Median and 95\% credible interval of the age-specific COVID-19 attributable deaths predicted by our approach |



## Installation 
A ```yml``` file is provided and can be used to build a conda virtual environment containing all dependencies. Create the environment using:
```bash
$ conda env create -f BSplinesProjectedGPs.yml
```
Then activate the environment for use:
```bash
$ source activate BSplinesProjectedGPs
```

This release has been checked on Ubuntu version 18.04.2.

## Warranty

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.


## Cite

Please cite 
* M Monod et al. Regularised B-splines projected Gaussian Process priors to estimate time-trends of age-specific COVID-19 deaths related to vaccine roll-out. (2021)
