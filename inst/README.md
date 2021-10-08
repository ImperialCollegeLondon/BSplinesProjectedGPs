## This folder contains...
| Folder    | Description |
|-----------|------------------------------------------------------|
| scripts   | R scripts to reproduce the results of the paper |
| functions | Functions needed to run the scripts |
| data      | Datasets, including cleaned CDC data |
| stan-models | Stan models of the weekly COVID-19 attributable deaths by age model |


## System Requirements
- macOS or UNIX, the code was developed on macOS Mojave 10.14
- [R](https://www.r-project.org/) version >= 3.6.1

## instructions 
The entry point to run the model on one US state is ```run-model-one-location-bash.sh``` and for all 50 US states it is ```run-model-bash.R```. 

### Header
In both files, you will need to specify the repository directory, the directory to store the results and the stan models under 
```bash
INDIR="repositorydirectory"
CWD="resultsdirectory"


### Usage: one US state 
To execute the model for one US state from the terminal console run, 
```bash
$ chmod +x run-model-one-location-bash.sh
$ ./run-model-one-location-bash.sh
```

### Usage: 50 US states on a laptop
To execute the model for 50 US states from the terminal console run, 
```bash
$ Rscript run-model-bash.R
```
This will print a JOBID and generate in the specified output directory one bash script for each state. Execute these bash script one by one, or alternatively run them in the background,
```bash
$ cd CWD/STAN_MODEL-JOBID
$ ./startme-STAN_MODEL-JOBID-1.sh 
```

### Usage: 50 US states in high-throughput computing
Run the Rscript
```bash
$ run-model-bash.R
```
that will generate a general submission script encapsulating the submission scripts of each state in individual PBS arrays.

You will need to modify the PBS header function at the start of this Rscript:
```R
#	function to make PBS header
make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal", hpc.q="pqcovid19c", hpc.array=1, hpc.log = NULL )
{	
  pbshead <- "#!/bin/sh"
  tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
  pbshead <- paste(pbshead, tmp, sep = "\n")
  tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")	
  pbshead <- paste(pbshead, tmp, sep = "\n")
  pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(hpc.array>1)
  {
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  }				
  if(!is.na(hpc.q))
  {
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  }		
  if(!is.null(hpc.log)){
    pbshead <- paste(pbshead, paste("#PBS -o", hpc.log), sep = "\n")
  }
  
  pbshead	<- paste(pbshead, hpc.load, sep = "\n")
  pbshead
}
```
