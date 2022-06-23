require(data.table)

cmdstan_dir <- '/rds/general/user/mm3218/home/git/cmdstan'
out_dir <- '/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results'
hmc_chains_n <- 8
source_dir <- '/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/'

if(0){
  out_dir <- '~/Downloads/results'
  cmdstan_dir <- '~/git/cmdstan'
  source_dir <- '~/git/BSplinesProjectedGPs/inst/'
}

#	function to make PBS header
make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal", hpc.array=1 )
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
  pbshead	<- paste(pbshead, hpc.load, sep = "\n")
  pbshead
}

if(1)
{
  countries <- paste0(c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME",
                        "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",
                        "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY"), collapse = ',')
  #countries <- 'CA,FL,NY,TX'
  
  n_countries <- length(unlist(strsplit(countries,',')))
  hpc.nproc.cmdstan <- n_countries	
  args <- data.table(
    source_dir= source_dir,
    cmdstan_dir = cmdstan_dir,
    out_dir= out_dir,
    script_file= 'scripts/run_stan_hpc.R',
    script_converting_file = "scripts/stan-convert-csv-to-rda.R",
    script_generate_quantities_file = "scripts/generate-quantities.R",
    stanModelFile= 'cmdstan_220616a',
    hmc_stepsize= 0.02,
    hmc_num_samples= 1500,
    hmc_num_warmup= 1000,			
    seed= 42,
    chain= 1,
    JOBID= round(runif(1,1,10000)),
    countries = countries,			
    cmdstan = 1L
  )	
}

if(exists('hpc.nproc.cmdstan'))
{
  stopifnot( hpc.nproc.cmdstan<=length(unlist(strsplit(args$countries, split=','))) )	
}

if(1)
{
  tmp <- data.table(chain=1:hmc_chains_n)		
  tmp[, seed:= round(runif(seq_len(nrow(tmp)))*1e6)]		
  set(args, NULL, colnames(tmp), NULL)
  tmp[, dummy:= 1L]
  args[, dummy:= 1L]
  args <- merge(args, tmp, by='dummy')
  set(args, NULL, 'dummy', NULL)	
}

# make commands
cmds <- vector('list', nrow(args))
for(i in seq_len(nrow(args)))
{
  cmd				<- ''			
  #	general housekeeping
  cmd				<- paste0(cmd,"CWD=$(pwd)\n")
  cmd				<- paste0(cmd,"echo $CWD\n")	
  tmpdir.prefix	<- paste0('cvd_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir			<- paste0("$CWD/",tmpdir.prefix)
  cmd				<- paste0(cmd,"mkdir -p ",tmpdir,'\n')	
  #	generate data set and run if not using cmdstan
  cmd 			<- paste0( cmd, 'echo "----------- Generating input data: ------------"\n')
  tmp 			<- paste0('Rscript ', file.path(args$source_dir[i],args$script_file[i]), 
                   ' -indir ', args$source_dir[i],'',
                   ' -outdir ', tmpdir,'',
                   ' -states "', args$countries[i],'"',
                   ' -stan_model "', args$stanModelFile[i],'"',
                   ' -JOBID ', args$JOBID[i]
  )
  cmd				<- paste0(cmd, tmp, '\n')
  #	if using cmdstan  
  if(args$cmdstan[i]==1) 
  {
    cmd <- paste0(cmd, 'echo "----------- Building Stan model file: ------------"\n')
    #	clean up any existing model code
    cmd <- paste0(cmd, 'rm ', file.path('$CWD',paste0("CDC-covid-tracker_", args$stanModelFile[i],'.*')), ' \n')
    #	copy stan model 
    cmd	<- paste0(cmd, 'cp -R ',file.path(args$source_dir[i], 'stan-models', paste0("CDC-covid-tracker_",args$stanModelFile[i],'.stan')), ' .\n')
    #	build model		
    cmd <- paste0(cmd, 'cd ', args$cmdstan_dir[i], '\n')
    cmd <- paste0(cmd, 'make ', file.path('$CWD',paste0("CDC-covid-tracker_",args$stanModelFile[i])), ' \n')
    cmd <- paste0(cmd, 'cd $CWD\n')
    #	set up env variables
    cmd <- paste0( cmd, 'STAN_DATA_FILE=$(find ', tmpdir, ' -name "*cmdstanin.R")\n')
    cmd <- paste0( cmd, 'JOB_DIR=${STAN_DATA_FILE%/*}\n')
    cmd <- paste0( cmd, 'JOB_DIR_NAME=${JOB_DIR##*/}\n')
    cmd <- paste0( cmd, 'STAN_OUT_FILE=', file.path('$JOB_DIR', '${JOB_DIR##*/}_stanout.csv'),' \n')
    #	run model
    cmd <- paste0( cmd, 'echo "----------- env variables are: ------------"\n')
    cmd <- paste0( cmd, 'echo $JOB_DIR\n')
    cmd <- paste0( cmd, 'echo $JOB_DIR_NAME\n')
    cmd <- paste0( cmd, 'echo $STAN_DATA_FILE\n')
    cmd <- paste0( cmd, 'echo $STAN_OUT_FILE\n')
    cmd <- paste0( cmd, 'echo "----------- Starting Stan sampling: ------------"\n')
    #	
    tmp <- paste0( './',paste0("CDC-covid-tracker_",args$stanModelFile[i]),' ',
                   'sample num_samples=',args$hmc_num_samples[i],' num_warmup=',args$hmc_num_warmup[i],' save_warmup=0 thin=1 ',
                   'adapt delta=0.99 ',
                   'algorithm=hmc engine=nuts max_depth=15 stepsize=',args$hmc_stepsize[i],' ',
                   'data file=$STAN_DATA_FILE ',
                   'random seed=',args$seed[i],' ',
                   'output file=$STAN_OUT_FILE' )
    cmd <- paste0(cmd, tmp, '\n')
    # convert csv to rdata
    cmd		<- paste0( cmd, 'echo "----------- Converting Stan output to RDA file: ------------"\n')
    tmp		<- paste0('Rscript ', file.path(args$source_dir[i],args$script_converting_file[i]), 
                   ' -csv_file "', "$STAN_OUT_FILE",'"',
                   ' -rda_file "', file.path('$JOB_DIR','${JOB_DIR##*/}_stanout.RData'),'"'
    )
    cmd		<- paste0(cmd, tmp, '\n')		
  }			
  
  #	general housekeeping
  cmd 	<- paste0( cmd, 'echo "----------- Copy files to out directory: ------------"\n')
  tmpdir2	<- file.path(args$out_dir[i], paste0(args$stanModelFile[i],'-',args$JOBID[i]))
  if(i==1)
  {
    dir.create(tmpdir2)		  		
  }
  cmd		<- paste0(cmd,"mkdir -p ",tmpdir2,'\n')
  cmd		<- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', file.path(tmpdir, paste0(args$stanModelFile[i],'-',args$JOBID[i])),'"/* ', tmpdir2,'\n')
  cmd		<- paste0(cmd, 'chmod -R g+rw ', tmpdir2,'\n')
  
  #	generate quantities 
  cmd <- paste0( cmd, 'echo "----------- Generating quantities: ------------"\n')
  cmd <- paste0( cmd, 'JOB_DIR2="',tmpdir2,'"/"$JOB_DIR_NAME" \n')
  cmd <- paste0( cmd, 'echo $JOB_DIR2\n')
  tmp <- length(unlist(strsplit(args$countries[i], split=',')))
  cmd <- paste0( cmd, paste0("echo {1..",tmp,"} | tr ' ' '\\n' | ") )
  tmp <- ifelse(args$cmdstan[i]==1, hpc.nproc.cmdstan, 1)
  stopifnot(is.numeric(tmp))
  cmd <- paste0( cmd, paste0('xargs -P ',tmp,' -n 1 -I {} ') )
  tmp <- paste0('Rscript ', file.path(args$source_dir[i],args$script_generate_quantities_file[i]),
                ' -indir "', args$source_dir[i], '"',
                ' -indir.results "$JOB_DIR2"',
                ' -location.index {}')		
  cmd <- paste0(cmd, tmp,'\n')
  
  # create post-processing shell script for central analyses
  if(i==1)
  {
    cmd2 <- make.PBS.header(	hpc.walltime=23, 
                             hpc.select=1, 
                             hpc.nproc=9, 
                             hpc.mem= "550gb", 
                             hpc.load= "module load anaconda3/personal", 
                             hpc.array= 1)
    cmd2 <- paste0(cmd2,'\n')
    # set up env variables
    cmd2 <- paste0(cmd2,'SCRIPT_DIR=',args$source_dir[i],'\n',			
                   'OUT_DIR=',dirname(tmpdir2),'\n',
                   'JOBID=',args$JOBID[i],'\n',
                   'STAN_MODEL_FILE=',args$stanModelFile[i],'\n',
                   'NUMB_CHAINS=', max(args$chain),'\n',
                   'OVERWRITE=0\n'
    )
    # save posterior samples
    tmp <- paste0('Rscript ', file.path('$SCRIPT_DIR','scripts','save-posterior-samples.R'),
                  ' -stanModelFile $STAN_MODEL_FILE -out_dir $OUT_DIR -JOBID $JOBID -numb_chains $NUMB_CHAINS')
    cmd2 <- paste0(cmd2,tmp,'\n')
    tmp = paste0('Rscript ', file.path('$SCRIPT_DIR','/scripts/postprocessing_assess_mixing.R'), 
                 ' -indir $SCRIPT_DIR -outdir $OUT_DIR -states ',  args$countries[i], ' -stan_model $STAN_MODEL_FILE -JOBID $JOBID')
    cmd2 = paste0(cmd2, tmp, '\n')
    tmp = paste0('Rscript ', file.path('$SCRIPT_DIR','/scripts/postprocessing_figures.R'), 
                 ' -indir $SCRIPT_DIR -outdir $OUT_DIR -states ',  args$countries[i], ' -stan_model $STAN_MODEL_FILE -JOBID $JOBID')
    cmd2 = paste0(cmd2, tmp, '\n')
    tmp = paste0('Rscript ', file.path('$SCRIPT_DIR','/scripts/postprocessing_union.R'), 
                 ' -indir $SCRIPT_DIR -outdir $OUT_DIR -states ',  args$countries[i], ' -stan_model $STAN_MODEL_FILE -JOBID $JOBID')
    cmd2 = paste0(cmd2, tmp, '\n')
    
    # write submission file	
    post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
    cat(cmd2, file=post.processing.file)
    # set permissions
    Sys.chmod(post.processing.file, mode='644')		
  }
  
  #	schedule post-processing	
  cmd		<- paste0( cmd, 'echo "----------- Post-processing: ------------"\n')
  tmp		<- paste("if [ $(find ",tmpdir2," -name '*_stanout.RData' | wc -l) -ge ",max( args$chain )," ]; then\n",sep='')
  cmd		<- paste(cmd,tmp,sep='')	
  post.processing.file <- file.path(tmpdir2, 'post_processing.sh')
  cmd 	<- paste0(cmd, '\tcd ', dirname(post.processing.file),'\n')
  cmd 	<- paste0(cmd,'\tqsub ', basename(post.processing.file),'\n')
  cmd		<- paste0(cmd,"fi\n")
  cmd		<- paste(cmd, "rm -rf $CWD/", basename(args$source_dir[i]),'\n',sep='')
  cmds[[i]]	<- cmd	
}	


pbshead <- make.PBS.header(	hpc.walltime=63, 
                            hpc.select=1, 
                            hpc.nproc=hpc.nproc.cmdstan, 
                            hpc.mem= paste0(hpc.nproc.cmdstan*9,'gb'), 
                            hpc.load= paste0("module anaconda3/personal\nexport STAN_NUM_THREADS=",hpc.nproc.cmdstan,"\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"),
                            hpc.array= length(cmds) )

#	make array job
for(i in seq_len(nrow(args)))
{
  cmds[[i]] <- paste0('echo PBS_JOBID_INDEX=', i,' > .Renviron\n', cmds[[i]])
  cmds[[i]] <- paste0(i,')\n',cmds[[i]],';;\n')
}
cmd		<- paste0('case $PBS_ARRAY_INDEX in\n',paste0(cmds, collapse=''),'esac')			
cmd		<- paste(pbshead,cmd,sep='\n')

#	submit job
outfile		<- gsub(':','',paste("cvd",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
outfile		<- file.path(args$out_dir[1], outfile)
cat(cmd, file=outfile)
cmd 		<- paste("qsub", outfile)
cat(cmd)
cat(system(cmd, intern= TRUE))	
cat(cmd)
