#!/bin/bash
# Run denoise.m and unring.m as a step 1 in the pipeline
# Ivan I. Maximov, Ullevål
# 22.11.2018, ver 1.0
# 11.12.2018, ver. 1.1

#SBATCH --account=p33_norment
#SBATCH --time=12:00:00
#SBATCH --job-name=UKstep1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

#source /cluster/software/fsl/6.0.1/setup.source
module load FSL/6.0.7.2
module load MATLAB/2017b

source ${FSLDIR}/etc/fslconf/fsl.sh

echo "Step 1: denoising and Gibbs ringing corrections ..."
# check input
if [ -z "${1}" ]; then
    echo "No input folder name. Stop."
    exit
else
    if [ -d "$1" ]; then
	echo "Processing $1 ..."
	IN="${1}"
    else
	echo "No such folder. Stop."
	exit
    fi
fi

echo "Ungzipping ..."
gzip -df ${IN}/*DTI_b*.nii.gz

echo "Probe AP"
AP=`ls ${IN}/*_b1000_*.bval`
if [ ! -e "${AP}" ]; then
    echo "There is no AP ${AP}"
    exit
else
    temp=${AP%%.*}
    AP="${temp}.nii"
    echo ${AP}
fi

echo "Probe PA"
PA=`ls ${IN}/*_b2800_*.bval`
if [ ! -e "${PA}" ]; then
    echo "There is no PA ${PA}"
    exit
else
    temp=${PA%%.*}
    PA="${temp}.nii"
    echo ${PA}
fi

echo "Denoise ${AP}"
matlab -nodesktop -nosplash -singleCompThread -r "try denoise('${AP}'); catch disp('DENOISE: Problem with ${AP}'); end; exit;"

echo "Denoise ${PA}"
matlab -nodesktop -nosplash -singleCompThread -r "try denoise('${PA}'); catch disp('DENOISE: Problem with ${PA}'); end; exit;"

AP=`ls ${IN}/*_b1000_*_n.nii`
echo "Unring ${AP}"
export LD_LIBRARY_PATH=/cluster/projects/p33/groups/imaging/ukbio/diffusion3/matlab/unring:${LD_LIBRARY_PATH}

PA=`ls ${IN}/*_b2800_*_n.nii`
echo "Unring ${PA}"
matlab -nodesktop -nosplash -singleCompThread -r "try unring('${PA}'); catch disp('UNRING: Problem with ${PA}'); end; exit;"

matlab -nodesktop -nosplash -singleCompThread -r "try unring('${AP}'); catch disp('UNRING: Problem with ${AP}'); end; exit;"

AP=`ls ${IN}/*_b1000_*_r.nii`
fslswapdim ${AP} -x y z ${IN}/AP_c

PA=`ls ${IN}/*_b2800_*_r.nii`
fslswapdim ${PA} -x y z ${IN}/PA_c

rm -f ${IN}/*_n.nii ${IN}/*_r.nii
gzip ${IN}/*.nii



