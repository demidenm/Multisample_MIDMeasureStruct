import glob
import sys
import pandas as pd
from nilearn.glm.second_level import SecondLevelModel


def group_onesample(fixed_files,group_out,contrast,smoothing=None):
    """
    This function takes in a list of fixed effect files for a select contrast and
    calculated a group (secondlevel) model by fitting an intercept to length of maps.
    For example, for 10 subject maps of contrast A, the design matrix would include an intercept length 10.

    :param fixed_files: list - a list of files to be used
    :param group_out: tr - output location to save z-statistic for group map
    :param smoothing: int - smoothing to perform at group level, assumed to have been performed at first level, default=None
    :param contrast: str - contrast label that is used in naming output file
    :return: nothing return, z-stat file saved to specified path
    """

    N_maps = len(fixed_files)

    # Create design matrix with intercept (1s) that's length of subjects/length of fixed_files
    design_matrix = pd.DataFrame([1] * N_maps,
                              columns=['Intercept'])

    # Fit secondlevel model
    sec_lvl_model = SecondLevelModel(mask_img=mni_mask, smoothing_fwhm=smoothing)
    sec_lvl_model = sec_lvl_model.fit(second_level_input=fixed_files,
                                      design_matrix=design_matrix)

    # Calculate z-statistic from second lvl map
    zstat_map = sec_lvl_model.compute_contrast(
        second_level_contrast='Intercept',
        second_level_stat_type='t',
        output_type='z_score',
    )

    # group out file, naming subs-N
    zstat_out = '{}/subs-{}_effect-onesample_contrast-{}_stat-zstat.nii.gz'.format(group_out, N_maps, contrast)
    zstat_map.to_filename(zstat_out)

    # Calculate t-statistic from second lvl map
    tstat_map = sec_lvl_model.compute_contrast(
        second_level_contrast='Intercept',
        second_level_stat_type='t',
        output_type='stat',
    )

    # group out file, naming subs-N
    tstat_out = '{}/subs-{}_effect-onesample_contrast-{}_stat-tstat.nii.gz'.format(group_out, N_maps, contrast)
    tstat_map.to_filename(tstat_out)


# Setting path to estimates for eahc run's glm; the fixed precision weighted out

fixed_input = sys.argv[1] # contrast files for subjects for fixed effect (Avg across runs)
group_out = sys.argv[2] # Output for group maps
mni_mask = '/oak/stanford/groups/russpold/data/MLS/code/proj_midinvar/ROI/MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz'

#contrast list
contrasts = ["Lgain-Neut", "LSgain-Neut", "Lgain-Lloss",
             "Lloss-Neut", "LSloss-Neut", "Lloss-Lgain"]


# Running group level act. [uncorrected] for list of contrasts
print('Starting on second level model (One-sample T-test, uncorrected)')
for con in contrasts:
    print(f'\t Working on contrast map: {con}')
    # find all fixed effect contrast maps across subjects
    fix_maps = sorted(glob.glob(f'{fixed_input}/**/*contrast-{con}_stat-effect.nii.gz'))
    group_onesample(fixed_files=fix_maps, group_out=group_out, contrast=con)
