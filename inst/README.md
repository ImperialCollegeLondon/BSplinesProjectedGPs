## This folder contains...
| Folder    | Description |
|-----------|------------------------------------------------------|
| scripts   | R scripts to reproduce the results of the paper |
| functions | Functions needed to run the scripts |
| data      | Datasets, including cleaned CDC data |
| stan-models | Stan models of the weekly COVID-19 attributable deaths by age model |


## Additional dependencies
Please install the following dependencies:
```bash
source activate BSplinesProjectedGPs
R -e 'options(unzip = "internal");devtools::install_github("zeehio/facetscales");install.packages("usmap")'
```

## Data
Datasets in ```data/``` generated in ```misc/``` include:
* All-ages daily deaths reported by JHU, ```JHU_data_$DATE.rds```
* Age-specific weekly deaths reported by the CDC, ```CDC-data_$DATE.rds```
* Vaccination data by age reported by the CDC, ```vaccination-prop-$DATE.rds```
* Vaccination data regardless of age reported by the CDC, ```vaccination-prop-pop-$DATE.rds```

Datasets extracted from other source include:
* Age-specific daily deaths reported by the DoH and extracted by the [Imperial College COVID-19 Response Team](https://github.com/ImperialCollegeLondon/US-covid19-agespecific-mortality-data), ```DeathsByAge_US_$DATE.csv```
* Population counts in the United States from the 2018 Census, ```us_population.csv```
* Share of COVID-19 deaths linked to nursing homes and other long term care facilities published by the [New York Time](https://www.nytimes.com/interactive/2020/us/coronavirus-nursing-homes.html), ```NYT_sharedeaths_carehomes.csv```

## Instructions 
The entry point to run the model is ```run-model.sh```.

### Header
You need to specify 
* the repository directory, ```INDIR```
* the output directory (to store the results), ```CWD```  
* the stan model ID under ```STAN_MODEL```
* the US state(s) you which to consider under ```STATES```
with,
```bash
INDIR="repositorydirectory"
CWD="resultsdirectory"
STAN_MODEL="stanmodelid"
STATES='CA'
```

Note the correspondence between the stan model ID and the model, 
| stan model id    | Model |
|-----------|------------------------------------------------------|
| 220209a | Regularised B-splines projected Gaussian Process  |
| 220209b   | Standard Gaussian Process |
| 220209c     | Standard B-splines |
| 220209d     | Bayesian P-splines |

For example, if you wish to use a standard Gaussian Process, specify
```bash
STAN_MODEL="220209b"
```

### Usage
From the repository directory, on the terminal console execute, 
```bash
$ source activate BSplinesProjectedGPs
$ cd inst/
$ ./run-model.sh
```
This attaches to your run a unique JOBID and generates in the output directory two bash scripts. One for running the model and one for processing the results. Execute these bash scripts one after the other,
```bash
$ cd $CWD/$STAN_MODEL-JOBID
$ ./$STAN_MODEL-JOBID.sh 
```
when finished, 
```bash
$ ./$STAN_MODEL-JOBID-postprocessing.sh 
```

