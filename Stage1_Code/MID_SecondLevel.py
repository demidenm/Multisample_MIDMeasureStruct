import glob
import pandas as pd
from nilearn.glm import compute_fixed_effects
from nilearn.glm.second_level import SecondLevelModel


def fixed_effect(sub, task_type, contrasts, input_path, output_path,
                 save_beta=True, save_var=True, save_tstat=True):
    """
    This function takes in a subject, task label, set of computed contrasts using nilearn,
    the path to contrast estimates (beta maps), the output path for fixed effec tmodels and
    specification of types of files to save, the beta estimates, associated variance and t-stat (which is calculated
    based on beta/variance values)
    Several path indices are hard coded, so should update as see fit
    e.g., '{sub}_ses-01_task-{task_type}_effect-fixed_contrast-{c}_stat-effect.nii.gz'

    :param sub: string-Input subject label, BIDS leading label, e.g., sub-01
    :param task_type: string-task type using bids label, e.g., 'mid'
    :param contrasts: list-list of contrast types that are saved from first level
    :param input_path: string-location of first level output files
    :param output_path: string-location to save fixed effects
    :param save_beta: Whether to save 'effects' or beta values, default = True
    :param save_variance: Whether to save 'variance' or beta values, default = True
    :param save_fixed: Whether to save 'tstat', default = True
    :return: nothing return, files are saved
    """
    for c in contrasts:
        betas = sorted(glob.glob(f'{input_path}/{sub}_*{c}_beta*.nii.gz'))
        var = sorted(glob.glob(f'{input_path}/{sub}_*{c}_var*.nii.gz'))

        # conpute_fixed_effects provides
        # (1) contrast map of the effect across runs;
        # (2) var map of between runs effect;
        # (3) t-statistic based on effect of variance;
        fix_effect, fix_var, fix_tstat = compute_fixed_effects(contrast_imgs=betas,
                                                               variance_imgs=var,
                                                               precision_weighted=True)

        if save_beta == True:
            fix_effect_out = f'{output_path}/{sub}_ses-01_task-{task_type}_effect-fixed_contrast-{c}_stat-effect.nii.gz'
            fix_effect.to_filename(fix_effect_out)

        if save_var == True:
            fix_var_out = f'{output_path}/{sub}_ses-01_task-{task_type}_effect-fixed_contrast-{c}_stat-var.nii.gz'
            fix_var.to_filename(fix_var_out)
        if save_tstat == True:
            fix_tstat_out = f'{output_path}/{sub}_ses-01_task-{task_type}_effect-fixed_contrast-{c}_stat-tstat.nii.gz'
            fix_tstat.to_filename(fix_tstat_out)


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
    sec_lvl_model = SecondLevelModel(smoothing_fwhm=smoothing)
    sec_lvl_model = sec_lvl_model.fit(second_level_input=fixed_files,
                                      design_matrix=design_matrix)

    # Calculate z-statistic from second lvl map
    zstat_map = sec_lvl_model.compute_contrast(
        second_level_contrast='Intercept',
        second_level_stat_type='t',
        output_type='z_score',
    )

    # group out file, naming subs-N
    zstat_out = f'{group_out}/subs-{N_maps}_effect-onesample_contrast-{contrast}_stat-zstat.nii.gz'
    zstat_map.to_filename(zstat_out)



# Set paths
data_dir = './Pilot_N3'
in_dir = f'{data_dir}/derivatives/fmri_prep/firstlvl_out/contrast_files'
fix_dir = f'{data_dir}/derivatives/fmri_prep/secondlevel/fixed'


# subjects and contrast list
sub = ['sub-21', 'sub-31']
contrasts = ["Lgain-Neut", "LSgain-Neut", "Lgain-Lloss",
             "Lloss-Neut", "LSloss-Neut", "Lloss-Lgain"]


# running fixed effect model (effect across two runs)
print(f'Step 1. Working on Fixed Effect Model')
for s in sub:
    print(f'\t Working on subject: {s}')
    fixed_effect(sub=s, contrasts=contrasts, task_type='mid',
                 input_path=in_dir, output_path=fix_dir,
                 save_tstat=True, save_beta=True, save_var=True)


# running group model (one sample t-test)
group_dir = f'{data_dir}/derivatives/fmri_prep/secondlevel/group'

# Running group level act. [uncorrected] for list of contrasts
print(f'Step 2. Starting on second level model (One-sample T-test, uncorrected)')
for con in contrasts:
    print(f'\t Working on contrast map: {con}')
    # find all contrast fixed effect contrast maps across subjects
    fix_maps = sorted(glob.glob(f'{fix_dir}/*contrast-{con}_stat-effect.nii.gz'))
    group_onesample(fixed_files=fix_maps,group_out=group_dir,contrast=con)
