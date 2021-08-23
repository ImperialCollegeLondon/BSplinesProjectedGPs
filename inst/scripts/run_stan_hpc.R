library(rstan)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

indir ="~/git/covid19Vaccination/inst" # path to the repo
outdir = file.path('~/Downloads/', "results")
location.index = 1
stan_model = "210823a"
JOBID = round(runif(1,1,1000))

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
path.to.CDC.data = file.path(indir, "data", paste0("CDC-data_2021-08-03.rds"))
path.to.JHU.data = file.path(indir, "data", paste0("jhu_data_2021-08-03.rds"))
path_to_scraped_data = file.path(indir, "data", paste0("DeathsByAge_US_2021-03-21.csv"))
path_to_vaccine_data = file.path(indir, "data", paste0("vaccination-prop-2021-08-21.rds"))

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

# load JHU and scraped data data
JHUData = readRDS(path.to.JHU.data)
scrapedData = read.csv(path_to_scraped_data)

# load vaccine data
vaccine_data = readRDS(path_to_vaccine_data)

# Create age maps
create_map_age(age_max)

# find locations 
locations = unique(select(deathByAge, loc_label, code)) 
saveRDS(locations, file = file.path(outdir.fit, paste0("location_", run_tag,".rds")))
loc_name = locations[location.index,]$loc_label
Code = locations[location.index,]$code
cat("Location ", as.character(loc_name), "\n")

# plot data 
if(0){
  plot_data(deathByAge = deathByAge, Code = Code, outdir = outdir.fig)
  plot_vaccine_data(deathByAge = deathByAge, vaccine_data = vaccine_data, outdir = outdir.fig)
  compare_CDC_JHU_DoH_error_plot(CDC_data = deathByAge,
                                    JHU_data = JHUData, 
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
stan_data = prepare_stan_data(deathByAge, loc_name, ref_date); data = tmp

if(grepl('210429a1|210429b1|210505b|210513a', stan_model)){
  cat("\n Using 1D splines \n")
  stan_data = add_1D_splines_stan_data(stan_data, spline_degree = 3, n_knots = 8)
}
if(grepl('210529b|210529c|210808|210823', stan_model)){
  cat("\n Using 2D splines \n")
  stan_data = add_2D_splines_stan_data(stan_data, spline_degree = 3, n_knots_rows = 12, n_knots_columns = 4)
}
if(grepl('210429a1|210429b1', stan_model)){
  cat("\n Adding adjacency matrix on 1D splines parameters \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$num_basis, m = stan_data$W)
}
if(grepl('210429a2|210429b2', stan_model)){
  cat("\n Adding adjacency matrix on 2D splines parameters \n")
  stan_data = add_adjacency_matrix_stan_data(stan_data, n = stan_data$num_basis_row, m = stan_data$num_basis_column)
}
if(grepl('210429a', stan_model)){
  cat("\n Adding nodes index \n")
  stan_data = add_nodes_stan_data(stan_data)
}
if(grepl('210429f|210429g', stan_model)){
  cat("\n With RW2 prior on splines parameters \n")
  stan_data = add_diff_matrix(stan_data, n = stan_data$num_basis, m = stan_data$W)
}
if(grepl('210808|210823', stan_model)){
  cat("\n With vaccine effects \n")
  stan_data = add_vaccine_prop(stan_data, df_week, Code, vaccine_data = vaccine_data)
}
if(1){
  cat("\n With Gamma prior for lambda \n")
  stan_data = add_prior_parameters_lambda(stan_data, distribution = 'gamma')
}

print("A = 12, W = 4")

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.data, paste0("stanin_", Code, "_",run_tag,".RData")) )


# fit 
cat("\n Start sampling \n")
model = rstan::stan_model(path.to.stan.model)

if(0){
  
  fit_cum <- rstan::sampling(model,data=stan_data,iter=100,warmup=10,chains=1,
                             seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99))
}


fit_cum <- rstan::sampling(model,data=stan_data,iter=2500,warmup=500,chains=8,
                           seed=JOBID,verbose=TRUE, control = list(max_treedepth = 15, adapt_delta = 0.99))

# save
file = file.path(outdir.fit, paste0("fit_cumulative_deaths_", Code, "_",run_tag,".rds"))
cat('\n Save file', file, '\n')
while(!file.exists(file)){
  tryCatch(saveRDS(fit_cum, file=file), error=function(e){cat("ERROR :",conditionMessage(e), ", let's try again \n")})
}




