#!/bin/bash
#SBATCH --account=p33
#SBATCH --time=48:00:00
#SBATCH --job-name=ABCD_2
#SBATCH --mem-per-cpu=8G
#
# TOPUP for DWI/IVIM data
# Ivan I. Maximov, Ullevål Oslo
# 24.08.2018, version 1.0
# 23.11.2018, ver 1.1

# Define variables here. All nii files should be noise/gibbs ringing corrected

source /cluster/bin/jobsetup
module purge

set -o errexit

module load FSL/6.0.7.2

source ${FSLDIR}/etc/fslconf/fsl.sh

TOPUP_PARAM="topup_params.txt"

if [ -z ${1} ]; then
	echo "No input folder ${1} for the correction. Please, check it."
	exit
else
    if [ -d ${1} ]; then
	IN="${1}"
    else
	echo "Cannot find folder ${1}"
	exit
    fi
fi

echo "Folder ${1}"

# TOPUP txt file
if [ ! -e ${TOPUP_PARAM} ]; then
	echo "No ${TOPUP_PARAM} topup file."
	exit
fi

echo "File check passed..."


mkdir -p ${IN}/AP
rm -f ${IN}/AP/*

mkdir -p ${IN}/PA
rm -f ${IN}/PA/*

mkdir -p ${IN}/TMP

TMP="${IN}/TMP"
AP="${IN}/AP"
PA="${IN}/PA"


if [ ! -e "${TMP}/b0_cor_fieldcoef.nii.gz" ]; then
    
    rm -f ${IN}/TMP/*
    
    fslsplit ${IN}/AP_c ${AP}/
    fslsplit ${IN}/PA_c ${PA}/

    flist=`ls ${AP}/*.nii.gz`
    for i in ${flist}; do
        fslslice ${i} ${TMP}/
        rm -f ${TMP}/_slice_0000.nii.gz
        fname=`basename -s .nii.gz ${i}`
        echo ${fname}
        fslmerge -a ${AP}/${fname}_n ${TMP}/*.nii.gz
        rm -f ${TMP}/*.nii.gz
    done

    fslmerge -a ${IN}/AP_n ${AP}/*_n.nii.gz
    rm -f ${AP}/*.nii.gz

    flist=`ls ${PA}/*.nii.gz`
    for i in ${flist}; do
        fslslice ${i} ${TMP}/
        rm -f ${TMP}/_slice_0000.nii.gz
        fname=`basename -s .nii.gz ${i}`
        echo ${fname}
        fslmerge -a ${PA}/${fname}_n ${TMP}/*.nii.gz
        rm -f ${TMP}/*.nii.gz
    done

    fslmerge -a ${IN}/PA_n ${PA}/*_n.nii.gz
    rm -f ${PA}/*.nii.gz

    fslroi ${IN}/AP_n ${TMP}/ap_2 0 2
    fslroi ${IN}/PA_n ${TMP}/pa_2 0 2
    fslmerge -a ${TMP}/b0all ${TMP}/ap_2 ${TMP}/pa_2
    rm -f ${TMP}/ap_2.nii.gz ${TMP}/pa_2.nii.gz
fi

# TOPUP
if [ ! -e "${TMP}/b0_cor_fieldcoef.nii.gz" ]; then
    ${FSLDIR}/bin/topup --imain=${TMP}/b0all --datain=${TOPUP_PARAM} --config=b02b0.cnf --out=${TMP}/b0_cor --iout=${TMP}/b0_new -v
fi
echo "... is done"

${FSLDIR}/bin/fslmaths ${TMP}/b0_new -Tmean ${TMP}/b0m
${FSLDIR}/bin/bet2 ${TMP}/b0m ${TMP}/brain -m -f 0.3

echo "Clean up temp files ..."
rm -f ${TMP}/brain.nii.gz

############################################
echo "eddy run"
BVAL=`ls ${IN}/*_b1000_*.bval`
BVEC=`ls ${IN}/*_b1000_*.bvec`
echo "BVAL ${BVAL}"
echo "BVEC ${BVEC}"

INDEX="index.txt"
PARAM="eddy_params.txt"

echo "eddy begins ..."
eddy_cpu --imain=${IN}/AP_n --mask=${TMP}/brain_mask --index=${INDEX} --acqp=${PARAM} --bvecs=${BVEC} --bvals=${BVAL} --out=${TMP}/data_eddy --topup=${TMP}/b0_cor --repol --data_is_shelled -v 

#############################################
echo "eddy run"
BVAL=`ls ${IN}/*_b2800_*.bval`
BVEC=`ls ${IN}/*_b2800_*.bvec`
echo "BVAL ${BVAL}"
echo "BVEC ${BVEC}"

INDEX="index.txt"
PARAM="eddy_params2.txt"

echo "eddy begins ..."
eddy_cpu --imain=${IN}/PA_n --mask=${TMP}/brain_mask --index=${INDEX} --acqp=${PARAM} --bvecs=${BVEC} --bvals=${BVAL} --out=${TMP}/data2_eddy --topup=${TMP}/b0_cor --repol --data_is_shelled -v 


rm -f ${IN}/AP_n.nii.gz ${IN}/PA_n.nii.gz

echo "Folder ${1} is finished"








