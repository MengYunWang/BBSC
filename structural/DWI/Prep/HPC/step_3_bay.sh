#!/bin/bash
# Post-process UKB diffusion data using our pipeline
#
# Author: Ivan I. Maximov
# version 1.0
# 30.01.2019

#SBATCH --account=p33_norment
#SBATCH --time=1:00:00
#SBATCH --job-name=BAYES
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

module load FSL/6.0.7.2
source ${FSLDIR}/etc/fslconf/fsl.sh

module load MATLAB/2020b

echo "We run SBATCH with $1"

# define variables here
RESN="${1}/RESULT"

gzip -df ${RESN}/data.nii.gz
gzip -df ${RESN}/brain_mask.nii.gz

cp -f ${1}/data.bval ${RESN}/
cp -f ${1}/data.bvec ${RESN}/

matlab -nodesktop -nosplash -singleCompThread -r "baymodelfit('${RESN}'); exit"

gzip -f ${RESN}/*.nii

fcor=`ls ${RESN}/data_Bayes_*.nii.gz`

for i in ${fcor}; do
    echo ${i}
    fslmaths ${i} -mul ${RESN}/brain_mask.nii.gz ${i}
done

echo "... done"

