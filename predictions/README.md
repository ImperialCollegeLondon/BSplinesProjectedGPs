## This folder contains...
| Folder    | Description |
|-----------|------------------------------------------------------|
| stan-models   | Stan models |
| data      | Training and test datasets  |


## System Requirements
- [R](https://www.r-project.org/) version >= 3.6.1


## Instruction 

### Header
In ```run_predictions.R```, specify the repository directory, the directory to store the results and the model under
```bash
indir = "repositorydirectory" 
outdir = "resultsdirectory"
model = "modelid"
```
Please note the correspondence between the model ID and the the model
| Model ID    | Mode |
|-----------|------------------------------------------------------|
| BS-GP-I   | Standard B-splines |
| BS-GP-SE      | Gaussian Process projected by regularized B-splines  |
For example, if you want to use standard B-splines, please specify, 
```bash
model = "BS-GP-I"
```

### Usage
To run the predictions, from the command line use
```bash
$ Rscript run_predictions.R
``` 
