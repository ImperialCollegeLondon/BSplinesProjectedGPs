#!/bin/sh

CWD="/Users/melodiemonod/Downloads/results/"
INDIR="/Users/melodiemonod/temp/covid19Vaccination/inst/"

JOBID=$$
STAN_MODEL="210529b"
LOCATION_INDEX=9
PWD=$(pwd)

# main directory
mkdir $PWD/$STAN_MODEL-$JOBID
  
# directories for fits, figures and tables
mkdir $PWD/$STAN_MODEL-$JOBID/fits
mkdir $PWD/$STAN_MODEL-$JOBID/data
mkdir $PWD/$STAN_MODEL-$JOBID/figure
mkdir $PWD/$STAN_MODEL-$JOBID/table
  
Rscript $INDIR/scripts/run_stan_hpc.R -indir $INDIR -outdir $PWD -location.index $LOCATION_INDEX -stan_model $STAN_MODEL -JOBID $JOBID
  
Rscript $INDIR/scripts/postprocessing_assess_mixing.R -indir $INDIR -outdir $PWD -location.index $LOCATION_INDEX -stan_model $STAN_MODEL -JOBID $JOBID
Rscript $INDIR/scripts/postprocessing_figures.R -indir $INDIR -outdir $PWD -location.index $LOCATION_INDEX -stan_model $STAN_MODEL -JOBID $JOBID
  
cp -R --no-preserve=mode,ownership $PWD/* $CWD

