# Regularised B-splines to quantify population-level trends in the age composition of COVID-19 attributable deaths before and after the 2021 vaccine drive in the United States

Mélodie Monod, Alexandra Blenkinsop, Andrea Brizzi, Yu Chen, Vidoushee Jogarah, Yuanrong Wang, Samir Bhatt and Oliver Ratmann (2021). Regularised B-splines to quantify population-level trends in the age composition of COVID-19 attributable deaths before and after the 2021 vaccine drive in the United States, LINK

## Abstract
XXX

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
$ cd CDC-covid19-agespecific-mortality-data
$ conda env create -f covid19Vaccination.yml
```
Then activate the environment for use:
```bash
$ source activate covid19Vaccination
```


## Warranty

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.


## Cite

Please cite 
* M Monod et al. Regularised B-splines to quantify population-level trends in the age composition of COVID-19 attributable deaths before and after the 2021 vaccine drive in the United States.
