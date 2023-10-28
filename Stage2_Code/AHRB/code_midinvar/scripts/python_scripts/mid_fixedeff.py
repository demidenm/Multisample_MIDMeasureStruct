import glob
import sys
import pandas as pd
from nilearn.glm import compute_fixed_effects


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


# Specify subject list
sub = sys.argv[1]

# Setting path to estimates for eahc run's glm; the fixed precision weighted out

first_input = sys.argv[2] # contrast files location
fixed_out = sys.argv[3] # fixed effect output location

#contrast list
contrasts = ["Lgain-Neut", "LSgain-Neut", "Lgain-Lloss",
             "Lloss-Neut", "LSloss-Neut", "Lloss-Lgain"]


# running fixed effect model (effect across two runs)
print(f'Start: Working on Precision Weighted Fixed Effect for {sub}')
fixed_effect(sub=sub, contrasts=contrasts, task_type='mid',
             input_path=first_input, output_path=fixed_out,
             save_tstat=True, save_beta=True, save_var=True)


