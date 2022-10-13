## This folder contains...
| Folder    | Description |
|-----------|------------------------------------------------------|
| stan-models   | Stan models |
| functions      | Functions needed to run the script   |
| data      | Simulated data   |


## Instruction 

### Header
In the Rscript ```run_simulations_2D.R```, specify the 
* repository directory ```indir``` and,
* the directory to store the results ```outdir``` 
with,
```R
indir = "repositorydirectory" 
outdir = "resultsdirectory"
```

### Usage
From the repository directory, on the terminal console execute, 
```bash
$ source activate BSplinesProjectedGPs
$ cd simulations
$ Rscript run_simulations_count_2D.R
``` 
