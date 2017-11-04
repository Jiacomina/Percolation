#PBS -l nodes=8:ppn=1
#PBS -m abe
#PBS -M 21719676@student.uwa.edu.au
source /etc/bash.bashrc
mpirun sim
