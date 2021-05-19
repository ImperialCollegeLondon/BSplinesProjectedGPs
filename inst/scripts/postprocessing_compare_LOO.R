
library(rstan)
library(data.table)
library(dplyr)
library(loo)

indir = "~/git/CDC-covid19-agespecific-mortality-data/inst" # path to the repo
outdir = '/rds/general/user/mm3218/home/git/CDC-covid19-agespecific-mortality-data/inst/results'
location.index = 1
stan_model = "210429b2"
stan_model.2 = "210505c1"
JOBID = 25741
JOBID.2 = 25996

# set directories
run_tag = paste0(stan_model, "-", JOBID)
run_tag.2 = paste0(stan_model.2, "-", JOBID.2)
outdir.table = file.path(outdir, run_tag, "table", run_tag)
outdir.table.2 = file.path(outdir, run_tag.2, "table", run_tag.2)


outdir.fit = file.path(outdir, run_tag, "fits")
locations = readRDS( file.path(outdir.fit, paste0("location_", run_tag,".rds")) )
Code = locations[location.index,]$code

LOO.1 = readRDS( paste0(outdir.table, "-LOO_", Code, ".rds") )
LOO.2 = readRDS( paste0(outdir.table.2, "-LOO_", Code, ".rds") )

comp <- loo_compare(LOO.1, LOO.2)
print(comp, simplify = FALSE, digits = 3)
