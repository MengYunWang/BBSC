# -*- coding: utf-8 -*-
"""

"""

import os

os.chdir('/Users/wang/Desktop/Research_projects/BBSC/Functional/Pre_processing/sub-1/ses_groupICA/ICA_15')
path = os.getcwd()

from nilearn import datasets, image, surface, plotting

fsaverage = datasets.fetch_surf_fsaverage()
ICA2plot = image.load_img('melodic_IC.nii.gz')

figs_path = os.path.join (path, 'figs')
isExist = os.path.exists(figs_path)
if not isExist:
  
  # Create a new directory because it does not exist 
  os.makedirs(figs_path)


for i, ica_data in enumerate(image.iter_img(ICA2plot)):

    ICA_surf_right = surface.vol_to_surf(ica_data, fsaverage.pial_right)
    ICA_surf_left = surface.vol_to_surf(ica_data, fsaverage.pial_left)


    plotting.plot_surf_stat_map(fsaverage.infl_right, ICA_surf_right,
                                 hemi='right', threshold=5., view='medial', colorbar=False,
                                 bg_map=fsaverage.sulc_right, output_file=os.path.join(figs_path, 'ICA%d_right_medial.png' % i)
                                 )

    plotting.plot_surf_stat_map(fsaverage.infl_right, ICA_surf_right, 
                                 hemi='right', threshold=5., view='lateral', colorbar=False,
                                 bg_map=fsaverage.sulc_right, output_file=os.path.join(figs_path, 'ICA%d_right_lateral.png' % i)
                                 )

    plotting.plot_surf_stat_map(fsaverage.infl_left, ICA_surf_left, 
                                 hemi='left', threshold=5., view='medial', colorbar=False,
                                 bg_map=fsaverage.sulc_left, output_file=os.path.join(figs_path, 'ICA%d_left_medial.png' % i)
                                 )

    plotting.plot_surf_stat_map(fsaverage.infl_left, ICA_surf_left, 
                                 hemi='left', threshold=5., view='lateral', colorbar=False,
                                 bg_map=fsaverage.sulc_left, output_file=os.path.join(figs_path, 'ICA%d_left_lateral.png' % i)
                                 )

    plotting.plot_glass_brain(ica_data, threshold=5, colorbar=True, plot_abs=False, display_mode='lyrz',
                               output_file=os.path.join(figs_path,'ICA%d.png' % i)
                               )
