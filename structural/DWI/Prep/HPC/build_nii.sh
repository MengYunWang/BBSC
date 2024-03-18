#!/bin/bash
#
#SBATCH --account=p33_norment
#SBATCH --time=1:00:00
#SBATCH --job-name=dcm2nii
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#


source /cluster/bin/jobsetup
module purge

set -o errexit
echo "We run SBATCH  ..."

echo "Start..."
echo "dcm 2 nifti ..."


cd ${1}
/cluster/projects/p33/users/ivanma/bin/dcm2niix -i y -z y ./
#rm -f *.dcm




