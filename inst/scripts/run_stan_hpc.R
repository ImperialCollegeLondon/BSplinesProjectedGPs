library(rstan)
library(data.table)
library(dplyr)
library(jcolors)
library(gridExtra)
library(ggpubr)

indir ="~/git/BSplinesProjectedGPs/inst" # path to the repo
outdir = file.path('~/Downloads/', "results")
states = 'CA'

stan_model = "220209a"
JOBID = 3541

if(0)
{
  outdir = '/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/results'
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
with_cmdstan <- F
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
path.to.stan.model = file.path(indir, "stan-models", paste0("CDC-covid-tracker_", stan_model, ".stan"))

# path to data
path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2022-02-06.rds"))
path.to.JHU.data = file.path(indir, "data", paste0("jhu_data_2022-02-06.rds"))
path_to_scraped_data = file.path(indir, "data", paste0("DeathsByAge_US_2021-03-21.csv"))
path_to_vaccine_data = file.path(indir, "data", paste0("vaccination-prop-2022-02-06.rds"))
path_to_vaccine_data_pop = file.path(indir, "data", paste0("vaccination-prop-pop-2022-02-06.rds"))
path.to.pop.data = file.path(indir, "data", paste0("us_population.csv"))
path.to.nyt.data = file.path(indir, "data", paste0("NYT_sharedeaths_carehomes.csv"))

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
if(!dir.exists(outdir.data)) dir.create( outdir.data, recursive = T)
if(!dir.exists( dirname(outdir.fig)) ) dir.create( dirname(outdir.fig), recursive = T)
print(outdir.fit)
print(outdir.data)
print(outdir.fig)
cat("\n outfile.dir is ", file.path(outdir, run_tag), '\n')

if(grepl('cmdstan', stan_model)){ 
  with_cmdstan = T
  tmp <- Sys.getenv("PBS_JOBID")
  tmp <- ifelse(tmp!='', gsub("(.+)\\[(.*)", "\\1", tmp), '' )
  job_id <- ifelse(tmp!='', tmp, as.character(abs(round(rnorm(1) * 1e6))) )
  tmp <- Sys.getenv("PBS_JOBID_INDEX")
  job_id <- ifelse(tmp!='', paste0(job_id, '[', tmp, ']'), job_id)
  job_dir <- file.path(outdir, run_tag, paste0(run_tag,'-',job_id)) 
  dir.create(job_dir)
  cat("\n job_dir is ", job_dir, '\n')
}

# load CDC data
deathByAge = readRDS(path.to.CDC.data) # cdc data 

# load JHU and scraped data data
JHUData = readRDS(path.to.JHU.data)
scrapedData = read.csv(path_to_scraped_data)

# load vaccine data
vaccine_data = readRDS(path_to_vaccine_data)
vaccine_data_pop = readRDS(path_to_vaccine_data_pop)

# load population count 
pop_data = as.data.table( read.csv(path.to.pop.data) )
setnames(pop_data, 'location', 'loc_label')

# load NYT data
nyt_data <- process.nyt_data(as.data.table(read.csv(path.to.nyt.data)))

# Create age maps
age_max = 105
create_map_age(age_max)

# find locations 
locations = unique(select(deathByAge, loc_label, code)) 
locations = locations[order(code)]
saveRDS(locations, file = file.path(outdir.fit, paste0("location_", run_tag,".rds")))
loc_name = locations[code %in% states,]$loc_label
Code = locations[code %in% states, ]$code
df_state = data.table(loc_label = loc_name, code = Code, state_index = 1:length(Code))
cat("Location ", as.character(loc_name), "\n")

# plot data 
if(1){
  p <- plot_data(deathByAge = deathByAge, Code = Code, outdir = outdir.fig)
  plot_vaccine_data(deathByAge = deathByAge, vaccine_data = vaccine_data[!date %in% max(date):(max(date) - 6)], pop_data = pop_data, Code, outdir = outdir.fig)
  compare_cumulative_CDC_JHU_DoH_error_plot(deathByAge = copy(deathByAge), JHUData = copy(JHUData), scrapedData = copy(scrapedData), outdir = outdir.fig)
  compare_weekly_CDC_JHU_DoH_error_plot(deathByAge = copy(deathByAge), JHUData = copy(JHUData), outdir = outdir.fig)
}


# reference date
ref_date = as.Date('2020-12-05')
cat("The reference date is", as.character(ref_date), "\n")

# Prepare stan data
cat("\n Prepare stan data \n")
stan_data = prepare_stan_data(deathByAge, loc_name, ref_date); data <- tmp
stan_data = add_JHU_data(stan_data, df_week, Code)
stan_data = add_vaccine_age_strata(stan_data, df_age_vaccination)
resurgence_dates <- find_resurgence_dates(JHUData, deathByAge, locations$code)[code %in% Code]
stan_data = add_resurgence_period(stan_data, df_week, resurgence_dates)
if(grepl('220209a|220209c|220209d|220607a|220208a|220131a|220615|220616|220617|220624|220701|220702|220712|220713', stan_model)){
  cat("\n Using 2D splines \n")
  knots_rows = c(df_age_reporting$age_from, max(df_age_continuous$age_to))
  stan_data = add_2D_splines_stan_data(stan_data, spline_degree = 3, n_knots_columns = 16, knots_rows = knots_rows)
}
if(grepl('220209d', stan_model)){
  cat("\n Adding adjacency matrix on 2D splines parameters \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$num_basis_row, m = stan_data$num_basis_column)
  stan_data = add_nodes_stan_data(stan_data)
}
if(1){
  cat("\n With Gamma prior for lambda \n")
  stan_data = add_prior_parameters_lambda(stan_data, distribution = 'exponential')
}
if(grepl('220208|220131|220615|220616|220617|220624|220701|220702|220712|220713', stan_model)){
  cat("\n With vaccine effects \n")
  stan_data = add_vaccine_prop(stan_data, df_week, Code, vaccine_data, resurgence_dates)
}
stan_init = add_init(stan_data)



## save image before running Stan
tmp <- names(.GlobalEnv)
if(with_cmdstan){
  file <- file.path(job_dir, paste0(basename(job_dir), '_stanin.RData'))
}else{
  # stan model
  model = rstan::stan_model(path.to.stan.model)
  
  file <- file.path(outdir.data, paste0("stanin_",run_tag,".RData"))
}
cat('Writing ', file, '\n')
save(list=tmp, file= file)

# fit 
cat("\n Start sampling \n")
if(0){
  nchain = 1
  fit_cum <- rstan::sampling(model,data=stan_data,iter=15,warmup=5,chains=nchain,
                             seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99), 
                             init = rep(list(stan_init), nchain))
}

if(with_cmdstan){
  file = file.path(job_dir, paste0(basename(job_dir),"_cmdstanin.R"))
  cat('Writing ', file, '\n')
  rstan::stan_rdump( names(stan_data), file=file, envir=list2env(stan_data))  
  
  file = file.path(job_dir, paste0(basename(job_dir),"_cmdstaninit.R"))
  cat('Writing ', file, '\n')
  rstan::stan_rdump( names(stan_init), file=file, envir=list2env(stan_init))
}else{
  library(foreach)
  library(doParallel)
  
  fit_cum <- rstan::sampling(model,data=stan_data,iter=1500,warmup=500,chains=8,
                             seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99), 
                             init = rep(list(stan_init), 8))
  
  # save
  file = file.path(outdir.fit, paste0("fit_cumulative_deaths_",run_tag,".rds"))
  
  cat('\n Save file', file, '\n')
  while(!file.exists(file)){
    tryCatch(saveRDS(fit_cum, file=file), error=function(e){cat("ERROR :",conditionMessage(e), ", let's try again \n")})
  }
}




