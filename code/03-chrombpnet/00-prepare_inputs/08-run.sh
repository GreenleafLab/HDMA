#!/bin/bash

# Purpose: runner script which executes all commands to get GC-matched negative
# or background training regions, using
# an array job to limit the number of jobs in the queue.
# 1) 08-write_cmds.sh first loops over datasets and generates all commands
# 2) 08-run.sh is then submitted to the scheduler, which creates an array job over all cmds in 08-commands.sh
# 3) 08-get_negatives.sh is called by each command.

#SBATCH --job-name=08-get_negatives.sh
#SBATCH --output=../../logs/03-chrombpnet/00/08/%x-%j.out
#SBATCH --partition=akundaje
#SBATCH --mem=20G
#SBATCH --array=1-1015:5
#SBATCH -n 5
#SBATCH --time=02:00:00

cmdfile="08-commands.sh"

for i in {0..4}; do
	id=$((SLURM_ARRAY_TASK_ID+i))
	cmd=$(sed -n "${id}p" $cmdfile)
	bash ${cmd} &
done

# # important to make sure the job doesn't exit before the background tasks are done
# # https://www.sherlock.stanford.edu/docs/advanced-topics/job-management/#minimizing-the-number-of-jobs-in-queue
wait



