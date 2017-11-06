#PBS -l nodes=1:ppn=12
#PBS -m abe
#PBS -M 21719676@student.uwa.edu.au
source /etc/bash.bashrc
mpirun sim
