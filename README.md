# Regularised B-splines projected Gaussian Process priors to estimate the age profile of COVID-19 deaths before and after vaccine roll-out

Mélodie Monod, Alexandra Blenkinsop, Andrea Brizzi, Yu Chen, Vidoushee Jogarah, Yuanrong Wang, Samir Bhatt and Oliver Ratmann (2021). Regularised B-splines projected Gaussian Process priors to estimate the age profile of COVID-19 deaths before and after vaccine roll-out, http://arxiv.org/abs/2106.12360.

## Abstract
The COVID-19 pandemic has caused severe public health consequences in the United States. The United States began a vaccination campaign at the end of 2020 targeting primarily elderly residents before extending access to younger individuals. With both COVID-19 infection fatality ratios and vaccine uptake being heterogeneous across ages, an important consideration is whether the age contribution to deaths shifted over time towards younger age groups. In this study, we use a Bayesian non-parametric spatial approach to estimate the age-specific contribution to COVID-19 attributable deaths over time. The proposed spatial approach is a low-rank Gaussian Process projected by regularised B-splines. Simulation analyses and benchmark results show that the spatial approach performs better than a standard B-splines approach and equivalently well as a standard Gaussian Process, for considerably lower runtimes. We find that COVID-19 has been especially deadly in the United States. The mortality rates among individuals aged 85+ ranged from 1\% to 5\% across the US states. Since the beginning of the vaccination campaign, the number of weekly deaths reduced in every US state with a faster decrease among individuals aged 75+ than individuals aged 0-74. Simultaneously to this reduction, the contribution of individuals age 75+ to deaths decreased, with important disparities in the timing and rapidity of this decrease across the country.

## About this repository
| Folder    | Description |
|-----------|------------------------------------------------------|
| misc   | Clean the incomplete age and state-specific mortality data reported by the CDC |
| inst | Reproduce the result of the "Real data analysis" section  |
| simulations      | Reproduce the result of the "Simulations" section |
| predictions | Reproduce the results of the "Real world benchmarks" |
| documents | paper |


## System Requirements
- [R](https://www.r-project.org/) version >= 3.6.1

## Installation 
A ```yml``` file is provided and can be used to build a conda virtual environment containing all R dependencies. Create the environment using:
```bash
$ conda env create -f BSplinesProjectedGPs.yml
```
Then activate the environment for use:
```bash
$ source activate BSplinesProjectedGPs
```

Conda is notorious for not working with RStan since it causes issues with the compiler linking. It is safer to install rstan and StanHeaders through R directly by following the steps [here](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-Source#mac) after loading the virtual environment.


## Warranty

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.


## Cite

Please cite 
* M Monod et al. Regularised B-splines to quantify population-level trends in the age composition of COVID-19 attributable deaths before and after the 2021 vaccine drive in the United States.
