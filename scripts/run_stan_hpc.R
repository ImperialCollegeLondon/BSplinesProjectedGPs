library(rstan)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

indir ="~/git/CDC-covid19-agespecific-mortality-data" # path to the repo
outdir = file.path('~/Downloads/', "results")
location.index = 4
stan_model = "210416a"
JOBID = round(runif(1,1,1000))

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
#options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
path.to.stan.model = file.path(indir, "stan-models", paste0("CDC-covid-tracker_", stan_model, ".stan"))

# path to CDC and JHU data
path.to.CDC.data.1 = file.path(indir, "data", paste0("CDC-data-1_2021-04-13.rds"))
path.to.CDC.data.2 = file.path(indir, "data", paste0("CDC-data-2_2021-04-13.rds"))
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
deathByAge_1 = readRDS(path.to.CDC.data.1) # cdc data before 2020-09-05 with first age specification 
deathByAge_2 = readRDS(path.to.CDC.data.2) # cdc data after 2020-09-05 with second age specification 
plot_data(deathByAge_1, deathByAge_2, outdir = outdir.fig)

# compare to JHU data
JHUData = readRDS(path.to.JHU.data)
compare_CDC_JHU_error_plot(CDC_data_1 = deathByAge_1, CDC_data_2 = deathByAge_2,
                           JHU_data = JHUData, 
                           var.daily.deaths.CDC = 'daily.deaths', 
                           outdir = outdir.fig)

# Create age maps
create_map_age(age_max)

# find locations 
locations = unique(deathByAge_1$loc_label) 
loc_name = locations[location.index]
cat("Location ", as.character(loc_name), "\n")

# reference date
ref_date = min(deathByAge_2$date)
cat("The reference date is", as.character(ref_date), "\n")

# Prepare stan data
cat("\n Prepare stan data \n")
stan_data_1 = prepare_stan_data(deathByAge_1, loc_name, ref_date); df_week1 = df_week; data1=tmp
stan_data_2 = prepare_stan_data(deathByAge_2, loc_name, ref_date, last_date_previous_spec); df_week2 = df_week; data2=tmp
stan_data = merge_stan_data(stan_data_1, stan_data_2)

if(grepl('210319d2|210319d3', stan_model)){
  cat("\n Using a GP \n")
  stan_data$age = matrix(stan_data$age, nrow = 106, ncol = 1)
}
if(grepl('210326|210329|210330|210406|210409|210412|210415|210416', stan_model)){
  cat("\n Using splines \n")
  stan_data = add_splines_stan_data(stan_data)
}
if(grepl('210406|210409|210412b|210415b|210416', stan_model)){
  cat("\n Adding adjacency matrix on splines parameters \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$W, m = stan_data$num_basis)
}
if(grepl('210408', stan_model)){
  cat("\n Adding adjacency matrix on week and age \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$W, m = stan_data$A)
}
if(grepl('210408b|210409|210412a|210412b|210415b|210416', stan_model)){
  cat("\n Adding nodes index \n")
  stan_data = add_nodes_stan_data(stan_data)
}
if(grepl('210416', stan_model)){
  cat("\n With RW2 prior on splines parameters \n")
  stan_data = add_diff_matrix(stan_data, n = stan_data$W, m = stan_data$num_basis)
}

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.data, paste0("stanin_", Code, "_",run_tag,".RData")) )


# fit 
cat("\n Start sampling \n")
model = rstan::stan_model(path.to.stan.model)
fit_cum <- rstan::sampling(model,data=stan_data,iter=2000,warmup=200,chains=3, seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.95))

# save
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_", Code, "_",run_tag,".rds"))
cat('\n Save file', file, '\n')
while(!file.exists(file)){
  tryCatch(saveRDS(fit_cum, file=file), error=function(e){cat("ERROR :",conditionMessage(e), ", let's try again \n")})
}

