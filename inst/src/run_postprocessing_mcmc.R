cat(" \n -------------------------------- \n \n Running postprocessing-figures.r \n \n -------------------------------- \n")

suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(bayesplot, quietly = TRUE))
suppressMessages(library(ggplot2, quietly = TRUE))
suppressMessages(library(ggnet, quietly = TRUE))
suppressMessages(library(coda, quietly = TRUE))
suppressMessages(library(gridExtra, quietly = TRUE))
library(scales)

indir <- "~/git/BSplinesProjectedGPs/inst"
stan_model <- '220209a'
JOBID <- '1081'
states = strsplit('FL,NY,TX,PA,IL,OH,GA,NC,MI',',')[[1]]

if(1){
  outdir <- "~/git/BSplinesProjectedGPs/inst/results"
}

if(0){
  outdir <- "/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results"
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
source(file.path(indir, 'src', 'functions', "postprocessing-mcmc.R"))

# set directories
run_tag = paste0(stan_model, "-", JOBID)
datadir.table = file.path(outdir, run_tag, "table", run_tag)
datadir.fit.post = file.path(outdir, run_tag, "fits")
datadir.data = file.path(outdir, run_tag, "data")
outdir.fit.post = file.path(outdir, run_tag, "fits")
outdir.fig.post = file.path(outdir, run_tag, "figure")
locs <- states

# load image 
load(file.path(datadir.data, paste0("stanin_",run_tag,".RData")))

# load posteror
file = paste0(datadir.table, '-vaccination_PosteriorSamples.rds')
MetrHastrw_outputs = readRDS(file)
posterior <- format_posterior(MetrHastrw_outputs)

# convergence and mixing analysis
posterior[, list(ESS = coda::effectiveSize(value) ), by = 'variable']
tryCatch(
  {
save_convergence_diagnostics(posterior, datadir.table)
  },
error=function(cond) {
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
}
)   

# confidence intervals
confidence_intervals <- save_confidence_intervals(posterior)

# forest plots
min_age_index_vac = 3
df_age_vaccination2 = df_age_vaccination[age_index >= 3]
df_age_vaccination2[, age_index := age_index - min_age_index_vac + 1]

names = c('slope_resurgence0', 'slope_resurgence_re', 'intercept_resurgence0', 'intercept_resurgence_re', 
          'vaccine_effect_intercept_cross', 'vaccine_effect_intercept_diagonal', 'vaccine_effect_intercept0',
          'vaccine_effect_slope_cross', 'vaccine_effect_slope_diagonal', 'vaccine_effect_slope0'
)
math_name = c('psi^"base"*""', 'psi^"state"*""', 'chi^"base"*""', 'chi^"state"*""', 
              'chi^"vacc-cross"*""', 'chi^"vacc"*""',  'chi^"vacc0"*""',
              'psi^"vacc-cross"*""', 'psi^"vacc"*""', 'psi^"vacc0"*""')
groups = c('slope', 'slope', 'baseline', 'baseline', 'baseline', 'baseline', 'baseline','slope', 'slope', 'slope')
groups_levels = c('baseline', 'slope')

tmp <- make_forest_plot_table(confidence_intervals, df_age_vaccination2, df_state, names, math_name, groups, groups_levels)
forest_plot <- plot_forest_plot(tmp, outdir.fig.post)


cat(" \n -------------------------------- \n \n Ending postprocessing-figures.r \n \n -------------------------------- \n")
