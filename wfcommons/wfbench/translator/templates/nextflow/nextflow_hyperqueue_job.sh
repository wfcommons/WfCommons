#!/bin/bash
#SBATCH --job-name=nextflowjob
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --account=<ACCOUNT>
#SBATCH --nodes=<NUM_NODES>
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=<MAX_NUM_CPU_PER_TASK>
#SBATCH --time=01:00:00

# Set variables
export PATH=$PATH:/$PWD
export NXF_DISABLE_CHECK_LATEST=true

# Set the directory which hyperqueue will use 
export HQ_SERVER_DIR=${PWD}/.hq-server
mkdir -p ${HQ_SERVER_DIR}

# Start the server in the background (&) and wait until it has started
hq server start &
until hq job list &>/dev/null ; do sleep 1 ; done

# Start the workers in the background and wait for them to start
srun --overlap --cpu-bind=none --mpi=none hq worker start --cpus=${SLURM_CPUS_PER_TASK} &
hq worker wait "${SLURM_NTASKS}"

# Ensure Nextflow uses the right executor and knows how many jobs it can submit
# The `queueSize` can be limited as needed. 
echo "executor {
  queueSize = $(( SLURM_CPUS_PER_TASK*SLURM_NNODES ))
  name = 'hq'
  cpus = $(( SLURM_CPUS_PER_TASK*SLURM_NNODES ))
}" > ${PWD}/nextflow.config

# run the Nextflow pipeline 
nextflow run workflow.nf -c nextflow.config --pwd $PWD &> output.log

# Wait for all jobs to finish, then shut down the workers and server
hq job wait all
hq worker stop all
hq server stop
