## This folder contains...
| Folder    | Description |
|-----------|------------------------------------------------------|
| stan-models   | Stan models |
| data      | Training and test datasets  |


## Instruction 

### Header
In the Rscript ```run_predictions.R```, specify 
* the repository directory ```indir```, 
* the directory to store the results ```outdir``` and,
* the model ID ```model``` 
with
```R
indir = "repositorydirectory" 
outdir = "resultsdirectory"
model = "modelid"
```
Please note the correspondence between the model ID and the method used
| Model ID    | Method |
|-----------|------------------------------------------------------|
| B-SPLINES_2D   | Standard B-splines |
| P-SPLINES_2D  | Bayesian P-splines |
| GP-B-SPLINES_2D     | Regularised B-splines projected Gaussian Process  |

For example, if you wish to use standard B-splines specify, 
```R
model = "B-SPLINES_2D"
```

### Usage
From the repository directory, on the terminal console execute, 
```bash
$ source activate BSplinesProjectedGPs
$ cd predictions/
$ Rscript run_predictions.R
``` 
