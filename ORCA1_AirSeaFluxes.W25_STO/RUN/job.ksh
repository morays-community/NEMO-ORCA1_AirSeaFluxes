#!/usr/bin/env bash

###################################
## Definitions for batch system
#SBATCH -A rkm@cpu                    # Accounting information
#SBATCH --job-name=ORCA1_W25_STO      # Job name
#SBATCH --partition=cpu_p1            # Partition Name
#SBATCH --ntasks=300                  # Total number of MPI processes
#SBATCH --hint=nomultithread          # 1  MPI process per node  (no hyperthreading)
#SBATCH --time=12:00:00                # Maximum execution time (HH:MM:SS)
#SBATCH --output=ORCA1_W25_STO.out    # Name of output listing file
#SBATCH --error=ORCA1_W25_STO.err     # Name of error listing file (the same)
###################################

# Process distribution
NPROC_NEMO=280
NPROC_PYTHON=20

## -------------------------------------------------------
##   End of user-defined section - modify with knowledge
## -------------------------------------------------------
set -x
ulimit -s unlimited

## Load environment
source ${HOME}/.bash_profile

# job information 
cat << EOF
------------------------------------------------------------------
Job submit on $SLURM_SUBMIT_HOST by $SLURM_JOB_USER
JobID=$SLURM_JOBID Running_Node=$SLURM_NODELIST 
Node=$SLURM_JOB_NUM_NODES Task=$SLURM_NTASKS
------------------------------------------------------------------
EOF

## Move to config directory
CONFIG_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}
cd ${CONFIG_DIR}

## Create execution directory and move there
XXD=`date +%F%H%M%S`
echo " XXD " $XXD
mkdir -p $CONFIG_DIR/OUT/$XXD
cd $CONFIG_DIR/OUT/$XXD
echo "RUN directory " `pwd`

## Get input files for NEMO
cp $CONFIG_DIR/grid_cfsites.nc . || exit 2
DATA1DIR=$SCRATCH/FORCING/ORCA1_CLIM
for file in $DATA1DIR/*
do
ln -s $file . || exit 2
done
for file in $DATA1DIR/BULK/*
do
ln -s $file . || exit 2
done
for file in $DATA1DIR/GRID/*
do
ln -s $file . || exit 2
done

## Get input namelist and xml files
for file in $CONFIG_DIR/*namelist*_ref $CONFIG_DIR/*namelist*_cfg $CONFIG_DIR/*.xml
do
    cp $file . || exit 3
done

## Get Executables
cp $CONFIG_DIR/nemo nemo.exe  || exit 5
cp $CONFIG_DIR/*.py . || exit 5

set -e
ls -l

## Run Eophis preproduction
python3 ./main.py --exec preprod
mv eophis.out eophis_preprod.out
mv eophis.err eophis_preprod.err

# write multi-prog file
touch run_file
echo 0-$((NPROC_NEMO - 1)) ./nemo.exe >> run_file
echo ${NPROC_NEMO}-$((NPROC_NEMO + NPROC_PYTHON - 1)) python3 ./main.py >> run_file

# run coupled NEMO-Python
time srun --multi-prog ./run_file

# Post-process
python3 ./compute_yearly_mean.py ORCA1_W25_STO_1y_20000101_20391231_ocebudget.nc ORCA1_W25_STO_40y_ocebudget.nc
python3 ./compute_yearly_mean.py ORCA1_W25_STO_1y_20000101_20391231_ocebudget.nc ORCA1_W25_STO_10y_ocebudget.nc 10
python3 ./compute_monthly_mean.py ORCA1_W25_STO_1m_20000101_20391231_grid_T.nc ORCA1_W25_STO_1m40y_averaged_grid_T.nc
python3 ./compute_monthly_mean.py ORCA1_W25_STO_1m_20000101_20391231_grid_T.nc ORCA1_W25_STO_1m10y_averaged_grid_T.nc 10

# Plots
python3 ./plots_res.py
python3 ./plots_diff.py
