#!/bin/bash
#SBATCH -A C3SE2022-1-16 -p vera
#SBATCH -n 20
#SBATCH -t 5-00:00:00
#SBATCH --mail-user=gustajo@chalmers.se
#SBATCH --mail-type=end  # send mail when job ends


# submit batch job using the line below:
# sbatch -o gen_gtex_ind-%A-%a.log --array=1-10 gen_gtex_ind.sh


module load MATLAB/2019a
module load GCCcore/10.3.0
module load GCCcore/11.2.0
module load Gurobi/9.5.0

matlab -nodesktop -nodisplay -r "generate_gtexind_models(${SLURM_ARRAY_TASK_ID}); exit" &

wait

