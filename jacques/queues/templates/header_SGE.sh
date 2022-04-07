#$ -N {name}
#$ -e {msgf}/{name}.msg
#$ -o {msgf}/{name}.msg
#$ -q {queue}
#$ -pe mp{cores} {cores}
#$ -R yes
#$ -cwd
