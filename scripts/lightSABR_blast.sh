#!/bin/bash

#SBATCH --nodes=1   # Number of nodes to use
#SBATCH --ntasks-per-node=1   # Use 1 processor core (single task, since BLAST is not parallelized)
#SBATCH --time=2-0:0:0   # Walltime limit (DD-HH:MM:SS)
#SBATCH --mem=32G   # Maximum memory per node
#SBATCH --job-name="lightSABR_BLAST_corn"   # Job name to display in squeue
#SBATCH --mail-user=bolivar@iastate.edu   # Email address
#SBATCH --mail-type=ALL   # Send an email when the job starts
#SBATCH --chdir="/work/adina/bolivar/lightSABR" # Set the working directory of the batch script to directory before it is executed. Absolute or relative paths
#SBATCH --output="/work/adina/bolivar/lightSABR/hpc/slurm-%j-corn_signif_asv_blast.out"   # Job standard output file (%j will be replaced by the slurm job id)
#SBATCH --error="/work/adina/bolivar/lightSABR/hpc/slurm-%j-corn_signif_asv_blast.out"   # Job standard error file (%j will be replaced by the slurm job id)


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK # Set OMP_NUM_THREADS to the number of CPUs per task we asked for.

##Modules/Singularity
module purge
#module load micromamba/1.4.2-lcemqbe # Latest in Nova
module load blast-plus/2.13.0-py310-irl7uxo

# Basic session info
echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS
echo sbatch invoked from: $SLURM_SUBMIT_DIR

# Run BLAST job
# Paths
PROJ_DIR=/work/adina/bolivar/lightSABR
IN=$PROJ_DIR/data/output/processed/sequences/sabr_2023_corn_signif_asv_fcm.fa
OUT=$PROJ_DIR/data/output/processed/sequences/sabr_2023_corn_signif_asv_blast.tsv

# Execution
./scripts/blast_nt_remote.sh "$IN" "$OUT"

# End job
module purge

echo End Job