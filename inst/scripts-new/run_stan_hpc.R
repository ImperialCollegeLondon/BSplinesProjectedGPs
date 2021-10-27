library(rstan)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

indir ="~/git/covid19Vaccination/inst" # path to the repo
outdir = file.path('~/Downloads/', "results")
states = strsplit('CA,FL,NY,TX,WA',',')[[1]]
stan_model = "211019b6a"
JOBID = 3541

if(0)
{
  outdir = '/rds/general/user/mm3218/home/git/covid19Vaccination/results'
  JOBID = 18389
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
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


# stan model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
path.to.stan.model = file.path(indir, "stan-models", paste0("CDC-covid-tracker_", stan_model, ".stan"))

# path to data
path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2021-09-29.rds"))
path.to.JHU.data = file.path(indir, "data", paste0("jhu_data_2021-09-30.rds"))
path_to_scraped_data = file.path(indir, "data", paste0("DeathsByAge_US_2021-03-21.csv"))
path_to_vaccine_data = file.path(indir, "data", paste0("vaccination-prop-2021-10-14.rds"))
path.to.pop.data = file.path(indir, "data", paste0("us_population_withnyc.rds"))

# load functions
source(file.path(indir, "functions-new", "summary_functions.R"))
source(file.path(indir, "functions-new", "plotting_functions.R"))
source(file.path(indir, "functions-new", "stan_utility_functions.R"))

# tag and directories
run_tag = paste0(stan_model, "-", JOBID)

outdir.fit = file.path(outdir, run_tag, "fits")
outdir.data = file.path(outdir, run_tag, "data")
outdir.fig = file.path(outdir, run_tag, "figure", run_tag)
if(!dir.exists(outdir.fit)) dir.create( outdir.fit, recursive = T)
if(!dir.exists( dirname(outdir.fig)) ) dir.create( dirname(outdir.fig), recursive = T)
print(outdir.fit)
print(outdir.data)
print(outdir.fig)
cat("\n outfile.dir is ", file.path(outdir, run_tag), '\n')

# load CDC data
deathByAge = readRDS(path.to.CDC.data) # cdc data 

# load JHU and scraped data data
JHUData = readRDS(path.to.JHU.data)
scrapedData = read.csv(path_to_scraped_data)

# load vaccine data
vaccine_data = readRDS(path_to_vaccine_data)

# load population count 
pop_data = as.data.table( reshape2::melt( readRDS(path.to.pop.data), id.vars = c('Region', 'code', 'Total')) )
setnames(pop_data, c('Region', 'variable', 'value'), c('loc_label', 'age', 'pop'))

# define resurgence period
start_resurgence <- as.Date('2021-07-17')
pick_resurgence <- as.Date('2021-09-04')

# Create age maps
age_max = 105
create_map_age(age_max)

# find locations 
locations = unique(select(deathByAge, loc_label, code)) 
saveRDS(locations, file = file.path(outdir.fit, paste0("location_", run_tag,".rds")))
loc_name = locations[code %in% states,]$loc_label
Code = locations[code %in% states, ]$code
df_state = data.table(loc_label = loc_name, code = Code, state_index = 1:length(Code))
cat("Location ", as.character(loc_name), "\n")

# plot data 
if(1){
  plot_data(deathByAge = deathByAge, Code = Code, outdir = outdir.fig)
  plot_vaccine_data(deathByAge = deathByAge, vaccine_data = vaccine_data, pop_data = pop_data, outdir = outdir.fig)
  compare_CDC_JHU_DoH_error_plot(CDC_data = deathByAge,
                                 JHUData = JHUData, 
                                 scrapedData = scrapedData,
                                 var.weekly.deaths.CDC = 'weekly.deaths', 
                                 outdir = outdir.fig,
                                 Code = Code)
}

# reference date
ref_date = as.Date('2020-12-05')
cat("The reference date is", as.character(ref_date), "\n")


# Prepare stan data
cat("\n Prepare stan data \n")
stan_data = prepare_stan_data(deathByAge, loc_name, ref_date); data <- tmp

if(grepl('211014|211019|211020|211025|211026|211027', stan_model)){
  cat("\n Using 2D splines \n")
  stan_data = add_2D_splines_stan_data(stan_data, spline_degree = 3, n_knots_rows = 12, n_knots_columns = 10)
}
if(grepl('211015', stan_model)){
  cat("\n Adding adjacency matrix on 2D splines parameters \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$num_basis_row, m = stan_data$num_basis_column)
}
if(grepl('211014b|211019|211020|211025|211026|211027', stan_model)){
  cat("\n With vaccine effects \n")
  stan_data = add_resurgence_period(stan_data, df_week, start_resurgence, pick_resurgence)
  stan_data = add_vaccine_prop(stan_data, df_week, Code, vaccine_data, start_resurgence, pick_resurgence)
  stan_data = add_JHU_data(stan_data, df_week, Code)
}
if(1){
  cat("\n With Gamma prior for lambda \n")
  stan_data = add_prior_parameters_lambda(stan_data, distribution = 'gamma')
}
if(grepl('211019b6a', stan_model)){
  cutoff_1864 = round(mean(range(stan_data$prop_vac_start[[1]])), 2)
  cutoff_65p =  round(mean(range(stan_data$prop_vac_start[[2]])), 2)
  stan_data = add_vaccine_prop_indicator(stan_data, cutoff_1864, cutoff_65p)
}
if(grepl('211025', stan_model)){
  cat("\n Add sequence of vaccinated \n")
  stan_data$prop_vac_sequence = seq(0, 1, 0.05)
  stan_data$P = length(stan_data$prop_vac_sequence)
}

print("A = 12, W = 10")

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.data, paste0("stanin_",run_tag,".RData")) )

# initial values
stan_init <- list()
stan_init$rho_gp1 <- rep(1.25, stan_data$M)
stan_init$rho_gp2 <- rep(1.25, stan_data$M)

# fit 
cat("\n Start sampling \n")
model = rstan::stan_model(path.to.stan.model)

if(0){
  
  fit_cum <- rstan::sampling(model,data=stan_data,iter=100,warmup=10,chains=1,
                             seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99))
}


fit_cum <- rstan::sampling(model,data=stan_data,iter=2500,warmup=500,chains=8,
                           seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99), 
                           init = rep(list(stan_init), 8))

# save
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_",run_tag,".rds"))
cat('\n Save file', file, '\n')
while(!file.exists(file)){
  tryCatch(saveRDS(fit_cum, file=file), error=function(e){cat("ERROR :",conditionMessage(e), ", let's try again \n")})
}



