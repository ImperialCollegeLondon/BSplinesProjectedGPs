#!/bin/sh

JOBID=$$
STAN_MODEL="220209a"
CWD="/Users/melodiemonod/git/BSplinesProjectedGPs/inst/results-try/"
INDIR="/Users/melodiemonod/git/BSplinesProjectedGPs/"
STATES='CA'

mkdir $CWD
  
cat > $CWD/bash_$STAN_MODEL-$JOBID.pbs <<EOF
  
#!/bin/sh
  
CWD=$CWD
INDIR=$INDIR/inst
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID
  
source activate BSplinesProjectedGPs

# main directory
mkdir \$CWD/\$STAN_MODEL-\$JOBID
  
# directories for fits, figures and tables
mkdir \$CWD/\$STAN_MODEL-\$JOBID/fits
mkdir \$CWD/\$STAN_MODEL-\$JOBID/data
mkdir \$CWD/\$STAN_MODEL-\$JOBID/figure
mkdir \$CWD/\$STAN_MODEL-\$JOBID/table
  
Rscript \$INDIR/scripts/run_stan_hpc.R -indir \$INDIR -outdir \$CWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
  
cd \$CWD

EOF

cat > $CWD/bash_$STAN_MODEL-$JOBID-postprocessing.pbs <<EOF
  
#!/bin/sh
  
CWD=$CWD
INDIR=$INDIR/inst
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID
  
source activate BSplinesProjectedGPs

# main directory
mkdir \$CWD/\$STAN_MODEL-\$JOBID
  
# directories for fits, figures and tables
mkdir \$CWD/\$STAN_MODEL-\$JOBID/fits
mkdir \$CWD/\$STAN_MODEL-\$JOBID/data
mkdir \$CWD/\$STAN_MODEL-\$JOBID/figure
mkdir \$CWD/\$STAN_MODEL-\$JOBID/table
  
Rscript \$INDIR/scripts/postprocessing_assess_mixing.R -indir \$INDIR -outdir \$CWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
Rscript \$INDIR/scripts/postprocessing_figures.R -indir \$INDIR -outdir \$CWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
Rscript \$INDIR/scripts/postprocessing_union.R -indir \$INDIR -outdir \$CWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID
  
EOF

chmod +x $CWD/bash_$STAN_MODEL-$JOBID.pbs
chmod +x $CWD/bash_$STAN_MODEL-$JOBID-postprocessing.pbs

echo "run model by executing: \$ $CWD/bash_$STAN_MODEL-$JOBID.pbs"
echo "then produce results by executing: \$ $CWD/bash_$STAN_MODEL-$JOBID-postprocessing.pbs"
