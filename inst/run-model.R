args = list(STAN_MODEL="210429a1",
            CWD="/rds/general/user/mm3218/home/git/covid19Vaccination/inst/results/",
            INDIR="/rds/general/user/mm3218/home/git/covid19Vaccination/inst/",
            locations = 1:50, 
            nchains = 8,
            JOBID = round(runif(1,1,10000)))


make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal", hpc.q="pqcovid19c", hpc.array=1, hpc.log = NULL )
{	
  pbshead <- "#!/bin/sh"
  tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
  pbshead <- paste(pbshead, tmp, sep = "\n")
  tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")	
  pbshead <- paste(pbshead, tmp, sep = "\n")
  pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if(hpc.array>1)
  {
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  }				
  if(!is.na(hpc.q))
  {
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  }		
  if(!is.null(hpc.log)){
    pbshead <- paste(pbshead, paste("#PBS -o", hpc.log), sep = "\n")
  }
  
  pbshead	<- paste(pbshead, hpc.load, sep = "\n")
  pbshead
}


pbshead = make.PBS.header(hpc.walltime = 30, hpc.nproc = args$nchains, hpc.mem = '100gb', 
                          hpc.array = length(args$locations), hpc.log = file.path(args$CWD, paste0(args$STAN_MODEL, '-', args$JOBID, '/')) )


pbshead = paste0(pbshead, '\n mkdir ', args$CWD, '/', args$STAN_MODEL,'-', args$JOBID, '\n')

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

  tmp = paste0('Rscript ', args$INDIR,'/scripts/postprocessing_union.R -indir ',args$INDIR, ' -outdir $CWD -stan_model ',args$STAN_MODEL, ' -JOBID $JOBID')
  cmds[[i]] = paste0(cmds[[i]], tmp, '\n')
  tmp = paste0('Rscript ', args$INDIR,'/scripts/knit_report.R -indir ',args$INDIR, ' -outdir $CWD -stan_model ',args$STAN_MODEL, ' -JOBID $JOBID')
  cmds[[i]] = paste0(cmds[[i]], tmp, '\n')
}


#	make array job
for(i in args$locations)
{
  cmds[[i]] <- paste0(i,')\n',cmds[[i]],';;\n')
}
cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')			
cmd		<- paste(pbshead,cmd,sep='\n')


#	submit job
outfile		<- file.path(args$CWD, paste0('bash_', args$STAN_MODEL,'-',args$JOBID,'.pbs'))
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))	


