#!/bin/bash
#SBATCH -p owners,iric
#SBATCH --exclusive
#SBATCH --job-name="$PWD"
#SBATCH --output=opt_work.log
#SBATCH --error=err_work.log
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --qos=normal
#SBATCH --mem-per-cpu=4000

NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`
NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`
NCPU=`echo " $NTASKS * $NNODES " | bc`
#module load openmpi/2.0.0
#module load python/2.7.5
module load anaconda/anaconda.4.4.0.python2.7
module load mpich/3.1.4/intel
#module load openmpi/1.8.
#mpiexec -n 16 
python test_tp_new.py
