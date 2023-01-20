#!/bin/bash

# dir=$1
# cd ./$1
ls *.zip > filelist.txt 
ls *.zip -l | sed -e 's/\.zip$//' > subjectlist.txt

for i in $(cat filelist.txt)  #unpack zipfiles in folder 
   do
   echo "Unpacking $i"
   tar -xvpf $i
done

for j in $(cat subjectlist.txt)  
do 
   ./dcm2niix $j/003_Sag_FSPGR_3D/
   ./dcm2niix $j/005_fMRI_default_restate/

done




#mkdir zipfiles # make folder to store HCP zip files in after they are unzipped
#mv *.zip* ./zipfiles   # move zip files to folder zipfiles.
# rm filelist.txt # remove filelist.txt
#ls | grep -E '^[0-9]+$' > subjectlist.txt # Make subjectlist.txt to use in HCP processing scripts
