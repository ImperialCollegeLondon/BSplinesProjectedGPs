library(rstan)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

indir ="~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path('~/Downloads/', "results")
location.index = 1
stan_model = "210429b2"
JOBID = round(runif(1,1,1000))

if(0)
{
  outdir = '/rds/general/user/mm3218/home/git/CDC-covid19-agespecific-mortality-data/results'
  JOBID = 18389
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-location.index')
  stopifnot(args_line[[7]]=='-stan_model')
  stopifnot(args_line[[9]]=='-JOBID')
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  location.index <- as.numeric(args_line[[6]])
  stan_model <- args_line[[8]]
  JOBID <- as.numeric(args_line[[10]])
}

# stan model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
path.to.stan.model = file.path(indir, "stan-models", paste0("CDC-covid-tracker_", stan_model, ".stan"))

# path to CDC and JHU data
path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2021-04-25.rds"))
path.to.JHU.data = file.path(indir, "data", paste0("jhu_death_data_padded_2021-04-11.rds"))

# load functions
source(file.path(indir, "functions", "summary_functions.R"))
source(file.path(indir, "functions", "plotting_functions.R"))
source(file.path(indir, "functions", "stan_utility_functions.R"))

# tag and directories
run_tag = paste0(stan_model, "-", JOBID)

outdir.fit = file.path(outdir, run_tag, "fits")
outdir.data = file.path(outdir, run_tag, "data")
outdir.fig = file.path(outdir, run_tag, "figure", run_tag)
if(!dir.exists(outdir.fit)) dir.create( outdir.fit, recursive = T)
if(!dir.exists( dirname(outdir.fig)) ) dir.create( dirname(outdir.fig), recursive = T)

cat("\n outfile.dir is ", file.path(outdir, run_tag), '\n')

# max age considered
age_max = 105
  
# load CDC data
deathByAge = readRDS(path.to.CDC.data) # cdc data 

# load JHU data
JHUData = readRDS(path.to.JHU.data)

# Create age maps
create_map_age(age_max)

# find locations 
locations = unique(select(deathByAge, loc_label, code)) 
saveRDS(locations, file = file.path(outdir.fit, paste0("location_", run_tag,".rds")))
loc_name = locations[location.index,]$loc_label
Code = locations[location.index,]$code
cat("Location ", as.character(loc_name), "\n")

# plot data 
if(1){
  plot_data(deathByAge = deathByAge, Code = Code, outdir = outdir.fig)
  compare_CDC_JHU_error_plot(CDC_data = deathByAge,
                             JHU_data = JHUData, 
                             var.daily.deaths.CDC = 'daily.deaths', 
                             outdir = outdir.fig)
}

# reference date
ref_date = as.Date('2020-08-29')
cat("The reference date is", as.character(ref_date), "\n")

# Prepare stan data
cat("\n Prepare stan data \n")
stan_data = prepare_stan_data(deathByAge, loc_name, ref_date); data = tmp

if(grepl('210319d2|210319d3', stan_model)){
  cat("\n Using a GP \n")
  stan_data$age = matrix(stan_data$age, nrow = 106, ncol = 1)
}
if(grepl('210426|210429', stan_model)){
  cat("\n Using splines \n")
  stan_data = add_splines_stan_data(stan_data, spline_degree = 3, n_knots = 8)
}
if(grepl('210426a|210426b|210426f|210426g|210429b|210429a', stan_model)){
  cat("\n Adding adjacency matrix on splines parameters \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$W, m = stan_data$num_basis)
}
if(grepl('210408', stan_model)){
  cat("\n Adding adjacency matrix on week and age \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$W, m = stan_data$A)
}
if(grepl('210416|210422a|210422e|210426a|210426g', stan_model)){
  cat("\n Adding nodes index \n")
  stan_data = add_nodes_stan_data(stan_data)
}
if(grepl('210416|210422d|210422e|210426f|210426g', stan_model)){
  cat("\n With RW2 prior on splines parameters \n")
  stan_data = add_diff_matrix(stan_data, n = stan_data$W, m = stan_data$num_basis)
}
if(grepl('210429a1|210429b1|210429d1', stan_model)){
  cat("\n With prior for lambda \n")
  stan_data = add_prior_parameters_lambda(stan_data, distribution = 'gamma')
}
if(grepl('210429a2|210429b2|210429d2', stan_model)){
  cat("\n With prior for lambda \n")
  stan_data = add_prior_parameters_lambda(stan_data, distribution = 'log_normal')
}

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.data, paste0("stanin_", Code, "_",run_tag,".RData")) )


# fit 
cat("\n Start sampling \n")
model = rstan::stan_model(path.to.stan.model)

if(0){
  
  fit_cum <- rstan::sampling(model,data=stan_data,iter=2000,warmup=200,chains=1,
                             seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99))
}


fit_cum <- rstan::sampling(model,data=stan_data,iter=2000,warmup=200,chains=8,
                           seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99))

# save
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_", Code, "_",run_tag,".rds"))
cat('\n Save file', file, '\n')
while(!file.exists(file)){
  tryCatch(saveRDS(fit_cum, file=file), error=function(e){cat("ERROR :",conditionMessage(e), ", let's try again \n")})
}



