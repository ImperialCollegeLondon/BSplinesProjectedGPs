#!/bin/sh

JOBID=$$
STAN_MODEL="211014"
CWD="/rds/general/user/mm3218/home/git/covid19Vaccination/inst/results/"
INDIR="/rds/general/user/mm3218/home/git/covid19Vaccination/inst/"
STATES='CA,TX'
  
cat > $CWD/bash_$STAN_MODEL-$JOBID.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=40:59:00
#PBS -l select=1:ncpus=8:ompthreads=1:mem=96gb
#PBS -j oe
module load anaconda3/personal
  
PWD=\$(pwd)
CWD=$CWD
INDIR=$INDIR
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID
  
# main directory
mkdir \$PWD/\$STAN_MODEL-\$JOBID
  
# directories for fits, figures and tables
mkdir \$PWD/\$STAN_MODEL-\$JOBID/fits
mkdir \$PWD/\$STAN_MODEL-\$JOBID/data
mkdir \$PWD/\$STAN_MODEL-\$JOBID/figure
mkdir \$PWD/\$STAN_MODEL-\$JOBID/table
  
Rscript \$INDIR/scripts/run_stan_hpc.R -indir \$INDIR -outdir \$PWD -states STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
  
Rscript \$INDIR/scripts/postprocessing_assess_mixing.R -indir \$INDIR -outdir \$PWD -states STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
Rscript \$INDIR/scripts/postprocessing_figures.R -indir \$INDIR -outdir \$PWD -states STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
  
cp -R --no-preserve=mode,ownership \$PWD/* \$CWD
  
EOF
  
cd $CWD
qsub bash_$STAN_MODEL-$JOBID.pbs
