#SBATCH -J {name}
#SBATCH -e {msgf}/{name}.msg
#SBATCH -o {msgf}/{name}.msg
#SBATCH -p {queue}
#SBATCH -N {nodes}
#SBATCH --ntasks-per-node={cores}
#SBATCH --mem={memory}
