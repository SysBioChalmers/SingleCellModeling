#!/bin/bash
#SBATCH -A C3SE2022-1-16 -p vera
#SBATCH -n 20
#SBATCH -t 5-00:00:00
#SBATCH --mail-user=gustajo@chalmers.se
#SBATCH --mail-type=end  # send mail when job ends


# instead of using the "SBATCH -o run_init-X.log" line here,
# include the flag "-o run_init-%A_%a.log" in the sbatch submission command:
# sbatch -o run_init_depmap-%A-%a.log --array=1-10 run_gen_depmap_models.sh

module load MATLAB/2019a
module load Gurobi/8.0.0

#use the upper line the first time it is run. If it generates some models and then for example runs out of time, the second line can be used to continue on the previous run
#matlab -nodesktop -nodisplay -nojvm -r "generate_DepMap_models('test'); exit" < /dev/null &
matlab -nodesktop -nodisplay -nojvm -r "generate_DepMap_models(${SLURM_ARRAY_TASK_ID}); exit" < /dev/null &
#matlab -nodesktop -nodisplay -nojvm -r "generate_DepMap_models(${SLURM_ARRAY_TASK_ID}, true); exit" < /dev/null &

wait

