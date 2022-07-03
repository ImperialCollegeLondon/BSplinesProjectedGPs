#!/bin/sh

JOBID="3068"
STAN_MODEL="220209a"
CWD="/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/inst/results/"
INDIR="/rds/general/user/mm3218/home/git/BSplinesProjectedGPs/"
#STATES='CA,FL,NY,TX'
#STATES='CA,FL,NY,TX,PA,IL,OH,GA,NC,MI,NJ,VA,WA,AZ,MA,TN,IN,MD,MO,WI,CO,MN,SC,AL,LA,KY,OR,UT,IA,NV,AR,MS,KS,NM,NE,ID,WV,HI,NH,ME'
STATES='CA,FL,NY,TX,PA,IL,OH,GA,NC,MI'

mkdir $INDIR

cat > $INDIR/inst/src/bash_$STAN_MODEL-$JOBID-mcmc.pbs <<EOF

#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal

CWD=$CWD
PWD=$CWD
INDIR=$INDIR/inst
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID

Rscript \$INDIR/src/run_mcmc.R -indir \$INDIR -outdir \$PWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID

cd \$INDIR/src
qsub bash_\$STAN_MODEL-\$JOBID-mcmc-postprocessing.pbs

EOF

cat > $INDIR/inst/src/bash_$STAN_MODEL-$JOBID-mcmc-postprocessing.pbs <<EOF

#!/bin/sh
#PBS -l walltime=40:59:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal

CWD=$CWD
PWD=$CWD
INDIR=$INDIR/inst
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID

Rscript \$INDIR/src/run_postprocessing_mcmc.R -indir \$INDIR -outdir \$PWD -states $STATES -stan_model \$STAN_MODEL -JOBID \$JOBID

EOF


cd $INDIR/inst/src
qsub bash_$STAN_MODEL-$JOBID-mcmc.pbs

