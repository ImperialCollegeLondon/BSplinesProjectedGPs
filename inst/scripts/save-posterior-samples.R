cat(" \n -------------------------------- \n \n Running save-posterior-samples.r \n \n -------------------------------- \n")

suppressMessages(library(rstan, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(abind, quietly = TRUE))

`%notin%` <- Negate(`%in%`)

# for dev purpose: Melodie
if(0){
  args_dir <- list()
  args_dir[['stanModelFile']] <- 'cmdstan_220616a'
  args_dir[['out_dir']] <- '~/Downloads/results/'
  args_dir[['JOBID']] <- '9228'
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
run_tag = paste0(args_dir$stanModelFile, "-", args_dir$JOBID)
args_dir$file.dir <- file.path(args_dir$out_dir, run_tag)

cat(" \n --------------------------------  with arguments -------------------------------- \n")
str(args_dir)

cat(" \n -------------------------------- check that HMC chains have run -------------------------------- \n")

do <- data.table(F=list.files(args_dir$file.dir, pattern='_stanout.RData', recursive=TRUE, full.name=TRUE))
do[, STANMF:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\1',basename(F))]
do[, JOBID:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\2',basename(F))]
do[, JOB_ID:= gsub('^([^-]+)-([^-]+)-([^-]+)$','\\3',basename(F))]
do <- subset(do, grepl(args_dir$stanModelFile,STANMF) & grepl(args_dir$JOBID, JOBID))
cat(paste("\n", nrow(do),"/",args_dir$numb_chains, "chains are finished \n"))

#if(nrow(do) != args_dir$numb_chains) stop()

outfile.base <- unique( do[, file.path(dirname(dirname(F)), paste0(STANMF,'-',JOBID))] )
outfile.base2 <- unique( do[, file.path(dirname(dirname(F)))] )

# args_dir[['JOBID']] <- do[1, gsub('^([^-]+)-([^-]+)-([^-]+)_([^-]+)$','\\3',basename(F))]

stopifnot(length(outfile.base)==1 )

cat(" \n -------------------------------- load jobs outputs -------------------------------- \n")

#	load all input variables for this analysis run
data_file <- gsub('_stanout.RData','_stanin.RData',do[1,F])
z <- load( data_file )

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

cat(" \n -------------------------------- load: generated quantities -------------------------------- \n")
do <- data.table(F=list.files(args_dir$file.dir, pattern='_stangqs.RDS', recursive=TRUE, full.name=TRUE))
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
  
  # permute the dimension to match previous dim
  if("lambda" %in% vars) re[["lambda"]] <- aperm(re[["lambda"]], c(1, 3, 2))
  if("phi" %in% vars) re[["phi"]] <- aperm(re[["phi"]], c(1, 4, 2, 3))
  if("alpha" %in% vars) re[["alpha"]] <- aperm(re[["alpha"]], c(1, 4, 2, 3))
  if("alpha_reduced" %in% vars) re[["alpha_reduced"]] <- aperm(re[["alpha_reduced"]], c(1, 4, 2, 3))
  if("phi_reduced" %in% vars) re[["phi_reduced"]] <- aperm(re[["phi_reduced"]], c(1, 4, 2, 3))
  if("phi_reduced_vac" %in% vars) re[["phi_reduced_vac"]] <- aperm(re[["phi_reduced_vac"]], c(1, 4, 2, 3))
  if("phi_predict_reduced_vac" %in% vars) re[["phi_predict_reduced_vac"]] <- aperm(re[["phi_predict_reduced_vac"]], c(1, 4, 2, 3))
  if("alpha_reduced_vac" %in% vars) re[["alpha_reduced_vac"]] <- aperm(re[["alpha_reduced_vac"]], c(1, 4, 2, 3))
  if("f" %in% vars) re[["f"]] <- aperm(re[["f"]], c(1, 4, 2, 3))
  if("E_pdeaths" %in% vars) re[["E_pdeaths"]] <- aperm(re[["E_pdeaths"]], c(1, 4, 2, 3))
  if("r_pdeaths" %in% vars) re[["r_pdeaths"]] <- aperm(re[["r_pdeaths"]], c(1, 4, 2, 3))
  if("log_lik" %in% vars) re[["log_lik"]] <- apply(re[["log_lik"]], 1:2, sum)
  if("deaths_predict" %in% vars) re[["deaths_predict"]] <- aperm(re[["deaths_predict"]], c(1, 4, 2, 3))
  if("deaths_predict_state_age_strata" %in% vars) re[["deaths_predict_state_age_strata"]] <- aperm(re[["deaths_predict_state_age_strata"]], c(1, 4, 2, 3))
  if("deaths_predict_vac_age_strata" %in% vars) re[["deaths_predict_vac_age_strata"]] <- aperm(re[["deaths_predict_vac_age_strata"]], c(1, 4, 2, 3))
  if("log_r_pdeaths_predict" %in% vars) re[["log_r_pdeaths_predict"]] <- aperm(re[["log_r_pdeaths_predict"]], c(1, 4, 2, 3))
  if("r_pdeaths_predict" %in% vars) re[["r_pdeaths_predict"]] <- aperm(re[["r_pdeaths_predict"]], c(1, 4, 2, 3))
  if("E_pdeaths_predict" %in% vars) re[["E_pdeaths_predict"]] <- aperm(re[["E_pdeaths_predict"]], c(1, 4, 2, 3))
  if("E_pdeaths_predict_resurgence_cumulative" %in% vars) re[["E_pdeaths_predict_resurgence_cumulative"]] <- aperm(re[["E_pdeaths_predict_resurgence_cumulative"]], c(1, 4, 2, 3))
  if("E_pdeaths_counterfactual" %in% vars) re[["E_pdeaths_counterfactual"]] <- aperm(re[["E_pdeaths_counterfactual"]], c(1, 2, 5, 3, 4))
  if("E_pdeaths_counterfactual_resurgence_cumulative" %in% vars) re[["E_pdeaths_counterfactual_resurgence_cumulative"]] <- aperm(re[["E_pdeaths_counterfactual_resurgence_cumulative"]], c(1, 2, 5, 3, 4))
  if("diff_E_pdeaths_counterfactual" %in% vars) re[["diff_E_pdeaths_counterfactual"]] <- aperm(re[["diff_E_pdeaths_counterfactual"]], c(1, 2, 5, 3, 4))
  if("perc_E_pdeaths_counterfactual" %in% vars) re[["perc_E_pdeaths_counterfactual"]] <- aperm(re[["perc_E_pdeaths_counterfactual"]], c(1, 2, 5, 3, 4))

  # name first column as iterations 
  names(dimnames(re[["lambda"]])) = c('iterations', rep('', length(dim(re[["lambda"]])) - 1))
  names(dimnames(re[["phi"]])) = c('iterations', rep('', length(dim(re[["phi"]])) - 1))
  names(dimnames(re[["alpha"]])) = c('iterations', rep('', length(dim(re[["alpha"]])) - 1))
  names(dimnames(re[["alpha_reduced"]])) = c('iterations', rep('', length(dim(re[["alpha_reduced"]])) - 1))
  names(dimnames(re[["alpha_reduced_vac"]])) = c('iterations', rep('', length(dim(re[["alpha_reduced_vac"]])) - 1))
  names(dimnames(re[["phi_reduced"]])) = c('iterations', rep('', length(dim(re[["phi_reduced"]])) - 1))
  names(dimnames(re[["phi_reduced_vac"]])) = c('iterations', rep('', length(dim(re[["phi_reduced_vac"]])) - 1))
  names(dimnames(re[["phi_predict_reduced_vac"]])) = c('iterations', rep('', length(dim(re[["phi_predict_reduced_vac"]])) - 1))
  names(dimnames(re[["f"]])) = c('iterations', rep('', length(dim(re[["f"]])) - 1))
  names(dimnames(re[["E_pdeaths"]])) = c('iterations', rep('', length(dim(re[["E_pdeaths"]])) - 1))
  names(dimnames(re[["r_pdeaths"]])) = c('iterations', rep('', length(dim(re[["r_pdeaths"]])) - 1))
  names(dimnames(re[["log_lik"]])) = c('iterations', rep('', length(dim(re[["log_lik"]])) - 1))
  names(dimnames(re[["deaths_predict"]])) = c('iterations', rep('', length(dim(re[["deaths_predict"]])) - 1))
  names(dimnames(re[["deaths_predict_state_age_strata"]])) = c('iterations', rep('', length(dim(re[["deaths_predict_state_age_strata"]])) - 1))
  names(dimnames(re[["deaths_predict_vac_age_strata"]])) = c('iterations', rep('', length(dim(re[["deaths_predict_vac_age_strata"]])) - 1))
  names(dimnames(re[["log_r_pdeaths_predict"]])) = c('iterations', rep('', length(dim(re[["log_r_pdeaths_predict"]])) - 1))
  names(dimnames(re[["r_pdeaths_predict"]])) = c('iterations', rep('', length(dim(re[["r_pdeaths_predict"]])) - 1))
  names(dimnames(re[["E_pdeaths_predict"]])) = c('iterations', rep('', length(dim(re[["E_pdeaths_predict"]])) - 1))
  names(dimnames(re[["E_pdeaths_predict_resurgence_cumulative"]])) = c('iterations', rep('', length(dim(re[["E_pdeaths_predict_resurgence_cumulative"]])) - 1))
  names(dimnames(re[["E_pdeaths_counterfactual"]])) = c('iterations', rep('', length(dim(re[["E_pdeaths_counterfactual"]])) - 1))
  names(dimnames(re[["E_pdeaths_counterfactual_resurgence_cumulative"]])) = c('iterations', rep('', length(dim(re[["E_pdeaths_counterfactual_resurgence_cumulative"]])) - 1))
  names(dimnames(re[["diff_E_pdeaths_counterfactual"]])) = c('iterations', rep('', length(dim(re[["diff_E_pdeaths_counterfactual"]])) - 1))
  names(dimnames(re[["perc_E_pdeaths_counterfactual"]])) = c('iterations', rep('', length(dim(re[["perc_E_pdeaths_counterfactual"]])) - 1))
  
  gc()
}


cat(" \n -------------------------------- processing job outputs -------------------------------- \n")



#
#	processing  quantities
cat(" \n -------------------------------- processing basic quantities: start -------------------------------- \n")

outdir.fit = file.path(args_dir$out_dir, run_tag,"fits")
outdir.data = file.path(args_dir$out_dir, run_tag,"data")

file <- file.path(outdir.fit, paste0("fit_cumulative_deaths_",run_tag,".rds"))
cat('\n Save file ', file)
saveRDS(fit, file = file)

file <- file.path(outdir.fit, paste0("posterior_samples_",run_tag,".rds"))
cat('\n Save file ', file)
saveRDS(re, file = file)

file <- file.path(outdir.data, paste0("stanin_",run_tag,".RData"))
cat('\n Save file ', file)
file.copy(data_file, file)

gc()
cat(" \n -------------------------------- processing basic quantities: end -------------------------------- \n")
