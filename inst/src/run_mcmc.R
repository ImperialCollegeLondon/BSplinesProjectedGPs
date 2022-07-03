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
# library("ramcmc")
# library("lognorm")
# library('metRology')

indir <- "~/git/BSplinesProjectedGPs/inst"
stan_model <- '220209a'
JOBID <- '1081'
states = strsplit('FL,NY,TX,PA,IL,OH,GA,NC,MI',',')[[1]]

if(0){
  outdir <- "~/git/BSplinesProjectedGPs/inst/results"
  datadir <- "~/git/BSplinesProjectedGPs/inst/results"
}

if(1){
  outdir <- "/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results"
  datadir <- "/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results"
}


## command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-states')
  stopifnot(args_line[[7]]=='-stan_model')
  stopifnot(args_line[[9]]=='-JOBID')
  
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  states <-  strsplit((args_line[[6]]),',')[[1]]
  stan_model <- args_line[[8]]
  JOBID <- as.numeric(args_line[[10]])
} 

# load functions
source(file.path(indir, 'src', "mcmc.R"))
source(file.path(indir, 'src', "mcmc-utils.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
datadir.table = file.path(outdir, run_tag, "table", run_tag)
datadir.fit.post = file.path(outdir, run_tag, "fits")
datadir.data = file.path(outdir, run_tag, "data")
outdir.fit.post = file.path(outdir, run_tag, "fits")
locs <- states

# load image 
load(file.path(datadir.data, paste0("stanin_",run_tag,".RData")))

# code
locations = readRDS( file.path(outdir.fit.post, paste0("location_", run_tag,".rds")) )
loc_name = locations[code %in% locs,]$loc_label
Code = locations[code %in% locs, ]$code
df_state = data.table(loc_label = loc_name, code = Code, state_index = 1:length(Code))

# resurgence deats
r_pdeaths  = vector(mode = 'list', length = length(Code))
for(i in seq_along(Code)){
  r_pdeaths[[i]] = readRDS(paste0(datadir.table, '-r_pdeathsPosteriorSamples_', Code[i], '.rds'))
}
r_pdeaths = do.call('rbind', r_pdeaths)
set(r_pdeaths, NULL, 'state_index', NULL)
r_pdeaths <- merge(r_pdeaths, df_state, by = c('code', 'loc_label'))

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

# select chain
# posterior[, chain := findInterval(iterations, seq(1, max_iterations,length.out = 8)), by = 'iterations']

# iterations and warmup
iter = length(r_pdeaths)
warmup = 5

# adapt 
adapt = T

# add data for all countries
stan_data = prepare_stan_data(deathByAge, loc_name, ref_date); data <- tmp
stan_data = add_JHU_data(stan_data, df_week, Code)
stan_data = add_vaccine_age_strata(stan_data, df_age_vaccination)
resurgence_dates <- find_resurgence_dates(JHUData, deathByAge, locations$code)[code %in% Code]
stan_data = add_resurgence_period(stan_data, df_week, resurgence_dates)
stan_data = add_vaccine_prop(stan_data, df_week, Code, vaccine_data, resurgence_dates)

# sampling
MetrHastrw_outputs = MetrHastrw(iter, warmup, r_pdeaths, stan_data, adapt)

# save
file = paste0(datadir.table, '-vaccination_PosteriorSamples.rds')
cat('\n Save file ', file, ' \n')
saveRDS(MetrHastrw_outputs, file = file)

