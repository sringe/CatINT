#!/bin/bash
#SBATCH -p owners,iric
#SBATCH --exclusive
#SBATCH --job-name="$PWD"
#SBATCH --output=opt_work.log
#SBATCH --error=err_work.log
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=0:10:00
#SBATCH --qos=normal
#SBATCH --mem-per-cpu=4000

NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`
NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`
NCPU=`echo " $NTASKS * $NNODES " | bc`
module load openmpi/2.0.0
mpiexec -n $NCPU python test_tp_new.py
