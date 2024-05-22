
# initialized parameters
set proj_mesh = fsaverage6
set down_mesh = fsaverage5
set sm = 6


# volume to surface
set regfile = ${subject}_bld${runfolder}${reg_stem}.dat
set output = $surffolder/$hemi.${BOLD}_$proj_short.nii.gz
mri_vol2surf --mov ${BOLD}.nii.gz --reg $regfile --hemi  $hemi \
--projfrac 0.5 --trgsubject $proj_mesh --o $output --reshape --interp trilinear

# smooth the data
set input = $surffolder/$hemi.${BOLD}_${proj_short}.nii.gz
set output = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}.nii.gz
mri_surf2surf --hemi $hemi --s $proj_mesh --sval $input --cortex \
--fwhm-trg 6 --tval $output --reshape

# adding back the medialwall 
set input1 = $surffolder/$hemi.${BOLD}_${proj_short}.nii.gz
set input2 = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}.nii.gz
set tmp_output = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}_fillin_medialwall.nii.gz
CBIG_preproc_fsaverage_medialwall_fillin '$hemi' 'fsaverage6' '$input1' '$input2' '$tmp_output'

# downsample to fsaverage5
set scale = $proj_res
while($scale > $down_res) 
				@ new_scale = $scale - 1
				if($scale == 7) then
					set srcsubject = fsaverage
				else
					set srcsubject = fsaverage$scale
				endif
				
				set trgsubject = fsaverage$new_scale
set input = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}.nii.gz
set curr_input = $input
set output = $surffolder/$hemi.${BOLD}_${proj_short}_sm${sm}_${down_short}.nii.gz
mri_surf2surf --hemi $hemi --srcsubject $srcsubject \
--sval $curr_input --nsmooth-in 1 --trgsubject $trgsubject --tval $output --reshape