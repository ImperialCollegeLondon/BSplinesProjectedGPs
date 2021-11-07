#!/bin/sh

JOBID=$$
STAN_MODEL="211031b1"
CWD="/Users/melodiemonod/git/covid19Vaccination/inst/results-try/"
INDIR="/Users/melodiemonod/git/covid19Vaccination/inst/"
STATES='CA,FL,NY,TX'
#STATES='CA'
  
cat > $CWD/bash_$STAN_MODEL-$JOBID.pbs <<EOF
  
#!/bin/sh
  
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
  
Rscript \$INDIR/scripts/run_stan_hpc.R -indir \$INDIR -outdir \$PWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
  
cp -R --no-preserve=mode,ownership \$PWD/* \$CWD
  
cd \$CWD

EOF

cat > $CWD/bash_$STAN_MODEL-$JOBID-postprocessing.pbs <<EOF
  
#!/bin/sh
  
CWD=$CWD
PWD=$CWD
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
  
Rscript \$INDIR/scripts/postprocessing_assess_mixing.R -indir \$INDIR -outdir \$PWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
Rscript \$INDIR/scripts/postprocessing_figures.R -indir \$INDIR -outdir \$PWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
Rscript \$INDIR/scripts/postprocessing_union.R -indir \$INDIR -outdir \$PWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
  
EOF

chmod +x $CWD/bash_$STAN_MODEL-$JOBID.pbs
chmod +x $CWD/bash_$STAN_MODEL-$JOBID-postprocessing.pbs



