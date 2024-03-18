#!/bin/bash
#SBATCH --account=p33_norment
#SBATCH --time=24:00:00
#SBATCH --job-name=ABCD_E
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#
# Metrics SMT/DKI/RSI/NODDI
# Ivan I. Maximov, Ullevål Oslo
# 24.08.2018, version 1.0
# 23.11.2018, ver 1.1

source /cluster/bin/jobsetup
module load MATLAB/2020b
module load FSL/6.0.7.2
source ${FSLDIR}/etc/fslconf/fsl.sh

# set -o errexit

RES="${1}/RESULT"
TMP="${1}/TMP"
AP="${1}/AP"
PA="${1}/PA"
rm -f ${AP}/* ${PA}/*

if [ ! -d "${1}/RESULT" ]; then
    mkdir -p ${1}/RESULT
    rm -f ${1}/RESULT/*
fi

# eddy corrected files to merge
if [ ! -e "${RES}/data.nii.gz" ]; then
    bl1=`ls ${1}/*_b1000_*.bval`
    bl2=`ls ${1}/*_b2800_*.bval`
    echo "BVAL 1000 ${bl1}"
    echo "BVAL 2800 ${bl2}"
    matlab -nodesktop -nosplash -nojvm -singleCompThread -r "try merge_bvecs('${bl1}','${bl2}','${TMP}/data_eddy.eddy_rotated_bvecs','${TMP}/data2_eddy.eddy_rotated_bvecs'); catch disp('Error in merge_bvecs.'); end; exit;"
    fslroi ${TMP}/data_eddy ${AP}/b1 0 1
    fslroi ${TMP}/data2_eddy ${AP}/b2 0 1
    flirt -in ${AP}/b2 -ref ${AP}/b1 -omat ${AP}/b1_2_b2.mat -out ${AP}/b2_n 
    fslsplit ${TMP}/data2_eddy ${PA}/
    flist=`ls ${PA}/*.nii.gz`
    for i in ${flist}; do
	k=`basename ${i}`
	flirt -in ${i} -ref ${AP}/b1 -applyxfm -init ${AP}/b1_2_b2.mat -out ${PA}/${k}
    done
    
    fslmerge -a ${RES}/data ${TMP}/data_eddy ${PA}/00??.nii.gz
    cp ${TMP}/brain_mask.nii.gz ${RES}/
fi

rm -f ${AP}/*
rm -f ${PA}/*

BVEC="${1}/data.bvec"
BVAL="${1}/data.bval"

export PATH=${PATH}:/cluster/projects/p33/users/ivanma/software/smt/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cluster/projects/p33/users/ivanma/software/smt/bin

if [ ! -e ${RES}/data_smt_b0.nii.gz ]; then
    fitmicrodt --bvecs ${BVEC} --bvals ${BVAL} --mask ${RES}/brain_mask.nii.gz ${RES}/data.nii.gz ${RES}/data_smt_{}.nii.gz
fi

if [ ! -e ${RES}/data_smt_mc_b0.nii.gz ]; then
    fitmcmicro --bvecs ${BVEC} --bvals ${BVAL} --mask ${RES}/brain_mask.nii.gz ${RES}/data.nii.gz ${RES}/data_smt_mc_{}.nii.gz
fi



if [ ! -e "${RES}/data_DKI_fa.nii.gz" ]; then
    echo "DKI"
    gzip -df ${RES}/data.nii.gz
    gzip -df ${RES}/brain_mask.nii.gz
    matlab -nodesktop -nosplash -singleCompThread -r "try DKI_estimation('${BVAL}','${BVEC}','${RES}/data.nii','${RES}/brain_mask.nii'); catch disp('DKI estimation failed'); end; exit"
    gzip -f ${RES}/*.nii
else
    echo "DKI is done"
fi

if [ ! -e "${RES}/data_RSI_CI.nii.gz" ]; then
    echo "RSI"
    gzip -df ${RES}/data.nii.gz
    gzip -df ${RES}/brain_mask.nii.gz
    matlab -nodesktop -nosplash -singleCompThread -r "try RSI_estimation('${RES}/data.nii','${BVEC}','${BVAL}','${RES}/brain_mask.nii'); catch disp('RSI estimation failed'); end; exit"
    gzip -f ${RES}/*.nii
else
    echo "RSI is done"
fi

if [ ! -e "${RES}/data_NODDI_OD.nii.gz" ]; then
    echo "NODDI"
    gzip -df ${RES}/data.nii.gz
    gzip -df ${RES}/brain_mask.nii.gz
    matlab -nodesktop -nosplash -singleCompThread -r "try NODDI_AMICO('${RES}/data.nii','${RES}/brain_mask.nii','${BVEC}','${BVAL}'); catch disp('NODDI estimation failed'); end; exit"
    mv -f ${RES}/AMICO/NODDI/FIT_ICVF.nii ${RES}/data_NODDI_IC.nii
    mv -f ${RES}/AMICO/NODDI/FIT_ISOVF.nii ${RES}/data_NODDI_ISO.nii
    mv -f ${RES}/AMICO/NODDI/FIT_OD.nii ${RES}/data_NODDI_OD.nii
    rm -fdr ${RES}/AMICO
    gzip -f ${RES}/*.nii
else
    echo "NODDI is done"
fi

echo "Folder ${1} is finished"








