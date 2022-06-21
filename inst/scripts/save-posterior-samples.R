cat(" \n -------------------------------- \n \n Running save-posterior-samples.r \n \n -------------------------------- \n")

suppressMessages(library(rstan, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(abind, quietly = TRUE))

`%notin%` <- Negate(`%in%`)

# for dev purpose: Melodie
if(0){
  args_dir <- list()
  args_dir[['stanModelFile']] <- 'base_age_fsq_mobility_201015i4_cmdstanv'
  args_dir[['out_dir']] <- '/rds/general/project/ratmann_covid19/live/age_renewal_usa/base_age_fsq_mobility_201015i4_cmdstanv-40states_Oct29_Levin7_schoolbound6_v2'
  args_dir[['JOBID']] <- '40states_Oct29_Levin7_schoolbound6_v2'
  args_dir[['numb_chains']] <- 8
}

# save args for report before loading those from running session 
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-stanModelFile')	
  stopifnot(args_line[[3]]=='-out_dir')
  stopifnot(args_line[[5]]=='-JOBID')
  stopifnot(args_line[[7]]=='-numb_chains')	
  args_dir <- list()
  args_dir[['stanModelFile']] <- args_line[[2]]
  args_dir[['out_dir']] <- args_line[[4]]
  args_dir[['JOBID']] <- args_line[[6]]
  args_dir[['numb_chains']] <- args_line[[8]]
} 

## start script
args_dir[['work_dir']] <- getwd()

cat(" \n --------------------------------  with arguments -------------------------------- \n")
str(args_dir)

cat(" \n -------------------------------- check that HMC chains have run -------------------------------- \n")

do <- data.table(F=list.files(args_dir$out_dir, pattern='_stanout.RData', recursive=TRUE, full.name=TRUE))
do[, STANMF:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\1',basename(F))]
do[, JOBID:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\2',basename(F))]
do[, JOB_ID:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\3',basename(F))]
do <- subset(do, grepl(args_dir$stanModelFile,STANMF) & grepl(args_dir$JOBID, JOBID))
cat(paste("\n", nrow(do),"/",args_dir$numb_chains, "chains are finished \n"))

#if(nrow(do) != args_dir$numb_chains) stop()

outfile.base <- unique( do[, file.path(dirname(dirname(F)), paste0(STANMF,'-',JOBID))] )
outfile.base2 <- unique( do[, file.path(dirname(dirname(F)))] )

args_dir[['JOBID']] <- do[1, gsub('^([^-]+)-([^-]+)-([^-]+)_([^-]+)$','\\3',basename(F))]

stopifnot(length(outfile.base)==1 )

cat(" \n -------------------------------- load jobs outputs -------------------------------- \n")

#	load all input variables for this analysis run
z <- load( gsub('_stanout.RData','_stanin.RData',do[1,F]) )

str(args)

#	reading job output, merge separate stanfits into one consolidated stanfit object
rf <- vector('list', nrow(do))
median_lp_ <- vector('numeric', nrow(do))
for(i in seq_len(nrow(do)))
{
  cat('Loading output in ',do[i,F],'\n')
  z <- load(do[i,F])
  stopifnot('fit' %in% z)
  median_lp_[i] = median(rstan:::extract(fit)$lp__)
  if(all(rstan::summary(fit)$summary[,1] == 0) & all(is.na(rstan::summary(fit)$summary[,2]))) next
  rf[[i]] <- fit
}
fit <- rstan:::sflist2stanfit(rf[lapply(rf,length)>0])
re <- rstan::extract(fit)


cat(" \n -------------------------------- save: fit -------------------------------- \n")
cat("\n save file:", paste0(outfile.base,'-stanout-fit.RDS'))
saveRDS(fit, file.path(args_dir[['work_dir']], args_dir[['out_dir']], basename(paste0(outfile.base,'-stanout-fit.RDS')))
fit <- NULL
gc()


cat(" \n -------------------------------- load: generated quantities -------------------------------- \n")
do <- data.table(F=list.files(args_dir$out_dir, pattern='_stangqs.RDS', recursive=TRUE, full.name=TRUE))
do[, STANMF:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\1',basename(F))]
do[, JOBID:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\2',basename(F))]
do[, JOB_ID:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\3',basename(F))]
do <- subset(do, grepl(args_dir$stanModelFile,STANMF) & grepl(args_dir$JOBID, JOBID))

do[, JOB_ID2:= gsub('.*\\[(.+)\\].*','\\1',basename(F))]

if(nrow(do)!=0){
  rf.gqs <- list()
  re.gqs = list()
  
  chains = as.numeric(unique(gsub(".*\\[(.+)\\].*", "\\1", do[,JOB_ID])))
  locations = as.numeric(unique(gsub(".*_location(.+)_.*", "\\1", do[,JOB_ID])))
  
  stan_data <- gqs_add_stan_data_for_flows(stan_data,dates)
  
  for(i in chains)
  {
    
    rf.gqs[[i]] <- list()
    
    for(Location in locations){
      index = which( grepl(paste0("_location", Location, "_"), do[,JOB_ID]) & grepl(paste0("\\[", i, "\\]"), do[,JOB_ID]) )
      cat('Loading output in ', do[index,F],'\n')
      rf.gqs[[i]][[Location]] <- readRDS(do[index,F])
      
      re.gqs[[i]] = list()
      for (var in names(rf.gqs[[i]][[Location]])){
        re.gqs[[i]][[var]] =  array(unlist(lapply(rf.gqs[[i]], "[[", var)), dim = c(dim(rf.gqs[[i]][[1]][[var]]), length(locations)))
      }
    }
  }
  
  vars = names(re.gqs[[chains[1]]])[names(re.gqs[[chains[1]]]) %notin% names(re)]
  
  for (var in vars){
    listvar = list( do.call("abind",list(lapply(re.gqs, "[[", var), along = 1)) )
    names(listvar) = var
    re <- c(re, listvar)
  }
  
    gc()
}


cat(" \n -------------------------------- processing job outputs -------------------------------- \n")



#
#	processing  quantities
cat(" \n -------------------------------- processing basic quantities: start -------------------------------- \n")

file <- paste0(outfile.base,'-posterior-samples.RDS')
saveRDS(re, file = file)

gc()
cat(" \n -------------------------------- processing basic quantities: end -------------------------------- \n")
