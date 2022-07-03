#### PACKAGES ####
library("extraDistr")
library("data.table")
library("LaplacesDemon")
library("plyr")
library("ggplot2")
library("compiler")
library("R.utils")
library("purrr")
library("doParallel")
library("prodlim")
library("tidyverse")
library("igraph")
library("ramcmc")
library("lognorm")
library('metRology')
library(abind)
library(compiler)

indir <- "~/git/BSplinesProjectedGPs"
stan_model <- '220209a'
JOBID <- '1081'

if(1){
  outdir <- "~/git/BSplinesProjectedGPs/src/results"
  datadir <- "~/git/BSplinesProjectedGPs/inst/results"
}

if(0){
  outdir <- "/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/src/results"
  datadir <- "/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results"
}


## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-datadir')
  stopifnot(args_line[[7]]=='-stan_model')
  stopifnot(args_line[[9]]=='-JOBID')
  
  indir <- args_line[[2]]  
  outdir <- args_line[[4]]
  datadir <- args_line[[6]]
  stan_model <- args_line[[8]]
  JOBID <- as.numeric(args_line[[10]])

} 

# load functions
source(file.path(indir, 'inst', 'src', "mcmc.R"))
source(file.path(indir, 'inst', 'src', "mcmc-utils.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
outdir.table = file.path(outdir, run_tag, "table", run_tag)

# code
locations = readRDS( file.path(outdir.fit.post, paste0("location_", run_tag,".rds")) )
loc_name = locations[code %in% states,]$loc_label
Code = locations[code %in% states, ]$code
df_state = data.table(loc_label = loc_name, code = Code, state_index = 1:length(Code))

# load image 
load(file.path(outdir.data, paste0("stanin_",run_tag,".RData")))

# resurgence deats
r_pdeaths  = vector(mode = 'list', length = length(locs))
for(i in seq_along(locs)){
  r_pdeaths[[i]] = readRDS(paste0(outdir.table, '-r_pdeathsPosteriorSamples_', locs[i], '.rds'))
}
r_pdeaths = do.call('rbind', r_pdeaths)

#format r_pdeaths
r_pdeaths <- r_pdeaths[order(iterations, state_index, age_index, week_index)]
r_pdeaths_list <- vector(mode = 'list', length = max(r_pdeaths$iterations))
for( i in sort(unique(r_pdeaths$iterations))){
  tmp <- array(0, dim = c(max(r_pdeaths$state_index), max(r_pdeaths$age_index), max(r_pdeaths$week_index)))
  for(j in sort(unique(r_pdeaths$week_index))){
    tmp[,,j] <- as.matrix(dcast(r_pdeaths[iterations == i & week_index == j], state_index~age_index, value.var = 'value')[, -1])
  }
  r_pdeaths_list[[i]] <- copy(tmp)
}
r_pdeaths <- copy(r_pdeaths_list)

# iterations and warmup
iter = max(r_pdeaths$iterations)
warmup = 5

# adapt 
adapt = T

# sampling
MetrHastrw_outputs = MetrHastrw(iter, warmup, r_pdeaths_list, stan_data, adapt)

# save
file = paste0(outdir.table, '-vaccination_PosteriorSamples.rds')
cat('\n Save file ', file, ' \n')
saveRDS(MetrHastrw_outputs, file = file)

