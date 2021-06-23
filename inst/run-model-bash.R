args = list(STAN_MODEL="210429h1",
            CWD="/rds/general/user/mm3218/home/git/covid19Vaccination/inst/results/",
            INDIR="/rds/general/user/mm3218/home/git/covid19Vaccination/inst/",
            locations = 1:50, 
            nchains = 8,
            JOBID = round(runif(1,1,10000)))

dir = paste0(args$CWD, '/', args$STAN_MODEL,'-', args$JOBID)
cat('directory is ', dir)
system(paste0('mkdir ', dir))

# PBS array
cmds = vector(mode = 'list', length = length(args$locations))
for(i in args$locations){
  
  # arguments
  
  cmds[[i]] = paste0('PWD=$(pwd)/temporary\n')
  cmds[[i]] = paste0(cmds[[i]], 'mkdir -p temporary\n')
  cmds[[i]] = paste0(cmds[[i]], 'JOBID=', args$JOBID, '\n')
  cmds[[i]] = paste0(cmds[[i]], 'STAN_MODEL=', args$STAN_MODEL, '\n')
  cmds[[i]] = paste0(cmds[[i]], 'CWD=', args$CWD, '\n')
  cmds[[i]] = paste0(cmds[[i]], 'INDIR=', args$INDIR, '\n')
  
  # main directories
  cmds[[i]] = paste0(cmds[[i]], 'mkdir $PWD/$STAN_MODEL-$JOBID\n')
  
  # directories for fits, figures and tables
  cmds[[i]] = paste0(cmds[[i]], 'mkdir $PWD/$STAN_MODEL-$JOBID/fits\n')
  cmds[[i]] = paste0(cmds[[i]], 'mkdir $PWD/$STAN_MODEL-$JOBID/data\n')
  cmds[[i]] = paste0(cmds[[i]], 'mkdir $PWD/$STAN_MODEL-$JOBID/figure\n')
  cmds[[i]] = paste0(cmds[[i]], 'mkdir $PWD/$STAN_MODEL-$JOBID/table\n')
  
  cmds[[i]] = paste0(cmds[[i]], 'echo "-- running model --"\n')
  tmp = paste0('Rscript ', args$INDIR,'/scripts/run_stan_hpc.R -indir ',args$INDIR, ' -outdir $PWD -location.index ', i, ' -stan_model ',args$STAN_MODEL, ' -JOBID $JOBID')
  cmds[[i]] = paste0(cmds[[i]], tmp, '\n')
  
  cmds[[i]] = paste0(cmds[[i]], 'echo "-- postprocessing --"\n')
  tmp = paste0('Rscript ', args$INDIR,'/scripts/postprocessing_assess_mixing.R -indir ',args$INDIR, ' -outdir $PWD -location.index ', i, ' -stan_model ',args$STAN_MODEL, ' -JOBID $JOBID')
  cmds[[i]] = paste0(cmds[[i]], tmp, '\n')
  tmp = paste0('Rscript ', args$INDIR,'/scripts/postprocessing_figures.R -indir ',args$INDIR, ' -outdir $PWD -location.index ', i, ' -stan_model ',args$STAN_MODEL, ' -JOBID $JOBID')
  cmds[[i]] = paste0(cmds[[i]], tmp, '\n')
  
  cmds[[i]] = paste0(cmds[[i]], 'cp -R --no-preserve=mode,ownership "$PWD"/* $CWD\n')
  
  tmp = paste0('Rscript ', args$INDIR,'/scripts/knit_report.R -indir ',args$INDIR, ' -outdir $CWD -stan_model ',args$STAN_MODEL, ' -JOBID $JOBID')
  cmds[[i]] = paste0(cmds[[i]], tmp, '\n')
  

  # write file	
  file <- paste0(args$CWD, '/', args$STAN_MODEL,'-', args$JOBID, '/', 'startme-', args$STAN_MODEL,'-', args$JOBID, '-',i,'.sh')
  cat(cmds[[i]], file=file)
  # set permissions
  Sys.chmod(file, mode='775')	
  
}






