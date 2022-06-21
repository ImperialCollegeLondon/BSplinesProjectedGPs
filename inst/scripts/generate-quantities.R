require(data.table)
require(rstan)
require(EnvStats)

## command line parsing if any
args = list()
args[['location.index']] = 1
args[['indir.results']] = '/rds/general/project/ratmann_covid19/live/age_renewal_usa/base_age_fsq_mobility_201015f8_cmdstanv-40states_tau10_Oct29_Levin/base_age_fsq_mobility_201015f8_cmdstanv-40states_tau10_Oct29_Levin-2558246[1].pbs'


args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{  
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-indir.results')
  stopifnot(args_line[[4]]=='-location.index')
  args <- list()    
  args[['indir']] <- args_line[[2]]  
  args[['indir.results']] <- args_line[[4]]  
  args[['location.index']] <- args_line[[6]]  
}

#	print args
str(args)

#	read args
indir.results <- args$indir.results
location.index <- as.integer(args$location.index)

#	read Stan input data and add location.index to parallelise computations
cat('\nReading Stan input data...')
infile.stanin <- list.files(indir.results, pattern='.*_stanin.RData$', recursive=TRUE)
stopifnot(length(infile.stanin)>0)
stopifnot(length(infile.stanin)<=1)
tmp <- load(file.path(indir.results, infile.stanin))
stopifnot(c('args','stan_data')%in%tmp)
stan_data$LOCATION_PROCESSING_IDX <- location.index

#	reset args
args[['work_dir']] <- getwd()

pkg.dir <- args$indir

#	read stanfits
cat('\nReading Stanfit ...')
infile.stanfits <- list.files(indir.results, pattern='.*_stanout.RData$', recursive=TRUE)
stopifnot(length(infile.stanfits)>0)
stopifnot(length(infile.stanfits)<=1)
tmp <- load(file.path(indir.results, infile.stanfits[1]))

#	find stan model and stan model for generating quantities
file_stanModel <- gsub('.*inst/(.*)','\\1',args$file_stanModel)
file_stanModel <- file.path(indir, file_stanModel)
file_stanModel_gqs <- gsub('cmdstan','gqs',file_stanModel)
stopifnot( file.exists(file_stanModel_gqs) )

cat('\nCompiling gqs model file ...')
m2 <- rstan::stan_model(file_stanModel_gqs)

cat('\nGenerating quantities ...')
draws <- as.matrix(fit)
fit2 <- rstan::gqs(m2, data=stan_data, draws=draws)
fit.gqs <- rstan::extract(fit2)

file <- file.path(indir.results, paste0(basename(args$job_dir), '_location',location.index,'_stangqs.RDS'))
saveRDS(fit.gqs, file = file.path( args[['work_dir']], dirname(file), basename(file)))

cat('\nFinished base-ages-generate-quantities.R ...')