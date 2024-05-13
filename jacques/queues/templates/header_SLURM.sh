#SBATCH -J {name}
#SBATCH -e {msgf}/{name}.msg
#SBATCH -o {msgf}/{name}.msg
#SBATCH -p {queue}
#SBATCH -A {account}
#SBATCH -N {nodes}
#SBATCH -c {cores}
#SBATCH --mem={memory}
#SBATCH --array={array_first}-{array_last}

ID=$SLURM_ARRAY_TASK_ID
