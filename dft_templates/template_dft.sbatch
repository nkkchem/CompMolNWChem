#!/bin/csh -f
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -t 00:30:00
#SBATCH -A emslc51198
#SBATCH -o ./INCHIKEY/dft/INCHIKEY.out
#SBATCH -e ./INCHIKEY/dft/INCHIKEY.err.%j
#SBATCH -J INCHIKEY
#SBATCH --export ALL

############################################################################
# Print out some information for refund purposes
############################################################################

echo "refund: UserID = thom510"
echo "refund: SLURM Job ID = ${SLURM_JOBID}"
echo "refund: Number of nodes          = 1"
echo "refund: Number of cores per node = 1"
echo "refund: Number of cores          = 1"
echo "refund: Amount of time requested = 0:30"
echo "refund: Directory = ${PWD}"
echo " "
echo Processor list
echo " "
echo "${SLURM_JOB_NODELIST}"
echo " "

############################################################################
# Actually run the job
############################################################################

source /etc/profile.d/modules.csh
module purge
module load nwchem/7.0.0

cd /scratch

setenv ARMCI_DEFAULT_SHMMAX 131072
#this disables xeon phi offload
setenv NWC_RANKS_PER_DEVICE 0
setenv OMP_NUM_THREADS 1
setenv MKL_NUM_THREADS 1
setenv NWC_RANKS_PER_DEVICE 0
setenv ARMCI_OPENIB_DEVICE mlx4_0
setenv OFFLOAD_INIT on_offload

setenv MPIRETURN 999
srun --mpi=pmi2 -n $SLURM_NPROCS -K1    /dtemp/scicons/bin/nwchem7.0.0 PATH/INCHIKEY.nw
setenv MPIRETURN $?

############################################################################
# End of the job script
############################################################################

exit $MPIRETURN

