# Regularised B-splines projected Gaussian Process priors to estimate time-trends in age-specific COVID-19 deaths

Mélodie Monod, Alexandra Blenkinsop, Andrea Brizzi, Yu Chen, Vidoushee Jogarah, Yuanrong Wang, Samir Bhatt and Oliver Ratmann (2022). Regularised B-splines projected Gaussian Process priors to estimate time-trends in age-specific COVID-19 deaths, **Bayesian analysis**, https://dx.doi.org/10.1214/22-BA1334.

## Abstract
The COVID-19 pandemic has caused severe public health consequences in the United States.
In this study, we use a hierarchical Bayesian model to estimate the age-specific COVID-19 attributable deaths over time in the United States. The model is specified by a novel non-parametric spatial approach over time and age, a low-rank Gaussian Process (GP) projected by regularised B-splines. We show that this projection defines a new GP with attractive smoothness and computational efficiency properties, derive its kernel function, and discuss the penalty terms induced by the projected GP. Simulation analyses and benchmark results show that the B-splines projected GP may perform better than standard B-splines and Bayesian P-splines, and equivalently well as a standard GP at considerably lower runtimes. 
We apply the model to weekly, age-stratified COVID-19 attributable deaths reported by the US Centers for Disease Control, which are subject to censoring and reporting biases. Using the B-splines projected GP, we can estimate longitudinal trends in COVID-19 associated deaths across the US by 1-year age bands. These estimates are instrumental to calculate age-specific mortality rates, describe variation in age-specific deaths across the US, and for fitting epidemic models. Here, we couple the model with age-specific vaccination rates to show that vaccination rates were significantly associated with the magnitude of resurgences in COVID-19 deaths during the summer 2021. With counterfactual analyses, we quantify the avoided COVID-19 deaths under lower vaccination rates and avoidable COVID-19 deaths under higher vaccination rates. The B-splines projected GP priors that we develop are likely an appealing addition to the arsenal of Bayesian regularising priors. 

## About this repository
| Folder    | Description |
|-----------|------------------------------------------------------|
| misc   | Obtain data to run the COVID-19 case study (age-specific weekly COVID-19 deaths reported by the CDC, all-ages weekly COVID-19 deaths reported by the JHU, vaccination data) |
| inst | Reproduce the result of the COVID-19 case study |
| simulations      | Reproduce the simulation analysis |
| predictions | Reproduce the benchmark results |

## Median and 95\% credible interval of the predicted age-specific COVID-19 attributable deaths adjusted for under-reporting and reporting-delays
In ```inst/results/predictions/predicted_weekly_deaths.rds```


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
* M Monod et al. (2022) Regularised B-splines projected Gaussian Process priors to estimate time-trends in age-specific COVID-19 deaths. Bayesian Analysis. https://doi.org/10.1214/22-ba1334
