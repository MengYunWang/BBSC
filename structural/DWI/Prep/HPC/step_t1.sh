#!/bin/bash
# Run denoise.m and unring.m as a step 1 in the pipeline
# Ivan I. Maximov, Skollenborg
# 22.11.2018, ver 1.0
# 11.12.2018, ver. 1.1
# 11.08.2019, ver. 2.0

#SBATCH --account=p33_norment
#SBATCH --time=72:00:00
#SBATCH --job-name=UKB_t1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

source /cluster/bin/jobsetup
module purge
module load FSL/6.0.7.2
source ${FSLDIR}/etc/fslconf/fsl.sh

set -o errexit

echo "Step 1: BEDPOSTX ..."
echo "FOLDER ${1}"

##########################################
# WRITE you actual path
##########################################
rr="/cluster/projects/p33/users/ivanma/Max/"

ln -sf ${rr}/${1}/RESULT/data.bval ${rr}/${1}/RESULT/bvals
ln -sf ${rr}/${1}/RESULT/data.bvec ${rr}/${1}/RESULT/bvecs
ln -sf ${rr}/${1}/RESULT/brain_mask.nii.gz ${rr}/${1}/RESULT/nodif_brain_mask.nii.gz

echo "Transformation to and back of MNI using FA image ..."
fsl_reg ${1}/RESULT/data_DKI_fa ${FSLDIR}/data/standard/FMRIB58_FA_1mm ${1}/RESULT/nat2std -e -FA
invwarp -r ${1}/RESULT/data_DKI_fa -w ${1}/RESULT/nat2std_warp -o ${1}/RESULT/std2nat_warp

echo "Create a reference volume ..."
flirt -in ${1}/RESULT/data_DKI_fa -ref ${1}/RESULT/data_DKI_fa -applyisoxfm 1 -out ${1}/RESULT/refVol
fslmaths ${1}/RESULT/refVol -mul 0 ${1}/RESULT/refVol

bedpostx ${1}/RESULT --model=3

echo "Done"





