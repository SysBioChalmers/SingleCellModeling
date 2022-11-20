#!/bin/bash
#SBATCH -A C3SE2022-1-16 -p vera
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=gustajo@chalmers.se
#SBATCH --mail-type=end  # send mail when job ends


# instead of using the "SBATCH -o run_init-X.log" line here,
# include the flag "-o run_init-%A_%a.log" in the sbatch submission command:
# sbatch -o logs/run_evaless_newalg2-%A-%a.log --array=1-40 run_evaless_newalg2.sh

module load MATLAB/2019a
module load GCCcore/10.3.0
module load GCCcore/11.2.0
module load Gurobi/9.5.0

#use the upper line the first time it is run. If it generates some models and then for example runs out of time, the second line can be used to continue on the previous run
matlab -nodesktop -nodisplay -r "getTaskEssentialGenes('DepMap/ftINIT2/depmap_models_newalg-', ${SLURM_ARRAY_TASK_ID}); exit" < /dev/null &

wait

