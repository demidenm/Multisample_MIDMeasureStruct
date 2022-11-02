import glob
import pandas as pd
import numpy as np
import nibabel as nib
from nilearn.glm.first_level import make_first_level_design_matrix
from nilearn.glm.first_level import FirstLevelModel
from nilearn.plotting import plot_design_matrix
from nilearn import image
import os


class color:
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


def regressors(confound_df, num_compcor=8):
    """

    :param confound_df: pandas dataframe for the *counfounds_timeseries.tsv
    :param num_compcor: the number of components to select from confounds.tsv, default first 8 [0:8]
    :return: list of confound regressors

    """
    # based on analyze_lev1.py from Jeanette mumford, creating confound_regressors
    # First, from confound df, pulling translation, rotation variables and cosine (see fMRIprep doc for using compcor regressors)
    #   from this df creating separate df for column names + values
    base_regressor_df = confound_df.filter(regex='trans|rot|cosine').copy()
    base_regressor_names = list(base_regressor_df.columns)
    base_regressors = base_regressor_df.values

    # Pulling compcor values
    add_regressor_names = []
    add_regressor_names += [i for i in confound_df.columns if
                            'a_comp_cor' in i][:num_compcor]  # if a_comp_cor, taking first 8
    additional_regressors = confound_df.loc[:, add_regressor_names].values

    # combining motion + compcor reg & creating names
    regressor_names = base_regressor_names + add_regressor_names
    regressors = np.hstack((base_regressors, additional_regressors))

    # creating pandas DF of regressors + names
    confound_regressors = pd.DataFrame(regressors, columns=regressor_names)

    return confound_regressors


def create_designmat(events_df, tr, volumes, conf_regressors, hrf_model='glover', stc=True):
    """
    :param events_df: this is the pandas dataframe for the events for the MID task
    :param tr: TR for the BOLD volume,
    :param volumes: volumes in the BOLD
    :param conf_regressors: dataframe of nuisance regressors from def(regressors)
    :param hrf_model: select hrf model for design matrix, default glover
    :param stc: whether slice time correction was done. To adjust the onsets/frame times in design matrix. Default True, alt False
    :return: returns a design matrix for each run with 5 anticipation,
            10 feedback, motion 6 + derivs, compcor + constant and confounds regressor file
    """

    # Adding the the events_df the cue + fixation duration for anticipation duration of model
    events_df["ANTICIPATION_DURATION"] = events_df["CUE_DURATION"] + events_df["FIXATION_DURATION"]

    # concatinating the condition types from events_df, anticipation type [:,0] & feedback type [:,9]
    # concatinating the onsets from events_df for anticipation type [:,1] & feedback type [:,7]
    # concatinating the duration from events_df for anticipation type (see line 18) [:,13] & feedback type [:,8]
    conditions = pd.concat([events_df.loc[:, "TRIAL_TYPE"], events_df.loc[:, "TRIAL_RESULT"]],
                           ignore_index=True)
    onsets = pd.concat([events_df.loc[:, "CUE_ONSET"], events_df.loc[:, "FEEDBACK_ONSET"]],
                       ignore_index=True)
    duration = pd.concat([events_df.loc[:, "ANTICIPATION_DURATION"], events_df.loc[:, "FEEDBACK_DURATION"]],
                         ignore_index=True)

    # create pandas df with events
    design_events = pd.DataFrame({'trial_type': conditions,
                                  'onset': onsets,
                                  'duration': duration})


    # creating design matrix using make_first_level_design_matrix from nilearn
    # Using the BOLD tr and volumes to generate the frame_times: acquisition time in seconds
    tr = tr
    vols = volumes
    frame_times = np.arange(vols) * tr

    if stc == True:
        design_matrix = make_first_level_design_matrix(# default modulation == '1'
            # Offset the times due to slice time correction, see blog post
            # https://reproducibility.stanford.edu/slice-timing-correction-in-fmriprep-and-linear-modeling /
            frame_times=frame_times+(tr/2), events=design_events,
            hrf_model=hrf_model, drift_model=None,
            add_regs=conf_regressors
        )
    elif stc == False:
        design_matrix = make_first_level_design_matrix(  # default modulation == '1'
            # Offset the times due to slice time correction, see blog post
            # https://reproducibility.stanford.edu/slice-timing-correction-in-fmriprep-and-linear-modeling /
            frame_times=frame_times, events=design_events,
            hrf_model=hrf_model, drift_model=None,
            add_regs=conf_regressors
        )

    return design_matrix

def glm_report(model, contrasts, out_dir, alpha=.05):
    """
    :param model: GLM modeled produced using FirstLevel
    :param contrasts: list of contrasts to run
    :param out_dir: str, path to write output with file name ending .html
    :param alpha: threshold activation maps for glass brain plotting, default = .05
    :return: saves html report out_dir path as .html
    """
    """
    
    """
    from nilearn.reporting import make_glm_report
    # https://nilearn.github.io/dev/modules/generated/nilearn.reporting.make_glm_report.html

    report = make_glm_report(model=model, contrasts=contrasts,
                             plot_type='glass', display_mode='lyrz',
                             height_control=None, alpha=alpha
                             # no cluster threshold on individual subj
                             )
    return report.save_as_html(out_dir)



# Setting data (for events.tsv) and derivative (for fmriprep) paths
data_path = '/Users/michaeldemidenko/sherlock/data/AHRB'
deriv_path = '/Users/michaeldemidenko/Desktop/Academia/Stanford/2_F32/Data/Pilot_N1/derivatives/fmri_prep'
firstlvl_output = f'{deriv_path}/firstlvl_out/contrast_files'

# Creating list of subjects and runs to loop over
sub = ['sub-31','sub-21']
runs = ['run-01','run-02']


# define contrasts for MID
contrasts = {
    'Lgain-Neut': 'LargeGain - NoMoneyStake',
    'LSgain-Neut': '.5*LargeGain + .5*SmallGain - 1*NoMoneyStake',
    'Lloss-Neut': 'LargeLoss - NoMoneyStake',
    'LSloss-Neut': '.5*LargeLoss + .5*SmallLoss - 1*NoMoneyStake',
    'Lgain-Lloss': 'LargeGain - LargeLoss',
    'Lloss-Lgain': 'LargeLoss - LargeGain'
}

# provide TR & volumes for associated BOLD
tr = .8
vols = 407


for s in sub:
    print(f'{color.DARKCYAN}Working on subject: {s}{color.END}')
    for r in runs:
        print(f'\tStarting {r} [Eight Steps]')

        print(f'\t\t 1/8 Load Files & set paths')
        # import behavior events .tsv from data path
        beh = pd.read_csv(f'{data_path}/{s}/ses-1/func/{s}_ses-1_task-mid_{r}_events.tsv', sep='\t')
        # get path to confounds from fmriprep & get data in dataframe filling missing values with 0 (e.g., derivatives)
        conf = glob.glob(
            f'{deriv_path}/{s}/ses-1/func/{s}_ses-1_task-mid_{r}_desc-confounds_timeseries.tsv'
        )[0]
        conf_df = pd.read_csv(conf, sep='\t', na_values=['n/a']).fillna(0)
        # Get path to functional mask
        mask_path = glob.glob(
            f'{deriv_path}/{s}/ses-1/func/{s}_ses-1_task-mid_{r}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz'
        )[0]
        # Get path to functional data
        nii_path = glob.glob(
            f'{deriv_path}/{s}/ses-1/func/{s}_ses-1_task-mid_{r}_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz'
        )[0]

        print(f'\t\t 2/8 Create Regressors & Design Matrix for GLM')
        # get list of regressors
        conf_regressors = regressors(confound_df = conf_df, num_compcor=8)

        # run to create design matrix
        design_matrix = create_designmat(events_df=beh, conf_regressors=conf_regressors,
                                         hrf_model='glover', tr=tr, volumes=vols,stc=True)

        print(f'\t\t 3/8 Fit GLM model with 5mm FWHM smoothing, ar1 autocorrelation')
        # fitting glm for sub's run, using associated mask, using ar1 autocorrelation (FSL prewhitening), drift model
        # 'cosine' and .01 highpass (100s filter)
        fmri_glm = FirstLevelModel(subject_label={s}, mask_img=mask_path, t_r=tr, smoothing_fwhm=5,
                                   standardize=False, noise_model='ar1', drift_model=None,high_pass=None
                                   # cosine 0:3 included from fmriprep in desing mat based on 128 s calc
        )

        # Run GLM model using set paths and calculate design matrix
        run_fmri_glm = fmri_glm.fit(nii_path, design_matrices=design_matrix)


        print(f'\t\t 4/8: From GLM model, create z-score contrast maps and save to output path')
        # contrast names and associated contrasts in contrasts defined is looped over
        # contrast name is used in saving file, the contrast is used in deriving z-score
        for con_name, con in contrasts.items():
            # zscore file
            zstat_name = f'{firstlvl_output}/{s}_ses-01_task-mid_{r}_contrast_{con_name}_zstat.nii.gz'
            z_est = run_fmri_glm.compute_contrast(con, output_type='z_score')
            z_est.to_filename(zstat_name)
            # To calculate ROI mean beta, need beta maps and beta maps with variance used in fixed effect calculation
            # fixed effect calculation is used for whole brain analysis inputs.
            # Calc: beta estimate
            beta_name = f'{firstlvl_output}/{s}_ses-01_task-mid_{r}_contrast_{con_name}_beta.nii.gz'
            beta_est = run_fmri_glm.compute_contrast(con, output_type='effect_size')
            beta_est.to_filename(beta_name)
            # Calc: variance
            var_name = f'{firstlvl_output}/{s}_ses-01_task-mid_{r}_contrast_{con_name}_var.nii.gz'
            var_est = run_fmri_glm.compute_contrast(con, output_type='effect_variance')
            var_est.to_filename(var_name)


        # This creates an HTML report for the GLM specified above.
        # The contrasts are printed at a .01 alpha
        print(f'\t\t 5/8: Using GLM model, generate design matrix with model information')
        plot_design_matrix(design_matrix=design_matrix,
                           output_file=f'{firstlvl_output}/{s}_ses-01_task-mid_{r}_model-FirstLevelDesign.png')

        # creating a filtered_bold signal, equivalent to filtered_func.nii.gz from FSL
        print(f'\t\t 6/8: Clean BOLD filter, e.g., detrend and remove signal related confound regressors '
              f'(e.g. consine, compcor, motion)')
        clean_bold = image.clean_img(nii_path, mask_img=mask_path, confounds=conf_regressors, t_r=tr,
                                     detrend=False, standardize=False, low_pass=None,ensure_finite=False
                                     # detrend (cosine 0:3) and regressors produced by fmriprep
        )

        # smoothing the above, cleaned bold image
        clean_bold_fwhm5mm = image.smooth_img(clean_bold, fwhm=5)
        # saving the filered image
        print('\t\t 7/8: Save filtered BOLD image to output path')
        nib.save(img= clean_bold_fwhm5mm,
                 filename=f'{firstlvl_output}/{s}_ses-01_task-mid_{r}_space-MNI152NLin2009cAsym_res-2_desc'
                          f'-filtered_bold.nii.gz'
                 )

        # Calculating SNR
        # to remove afni warning: 3dinfo -DAFNI_NIFTI_TYPE_WARN=NO
        print('\t\t 8/8: Calculate SNR using 3dTstat')
        tsnr_run = f"""3dTstat \
            -prefix {firstlvl_output}/{s}_ses-01_task-mid_{r}_space-MNI152NLin2009cAsym_res-2_desc-filteredSNR_bold.nii \
            -tsnr {firstlvl_output}/{s}_ses-01_task-mid_{r}_space-MNI152NLin2009cAsym_res-2_desc-filtered_bold.nii.gz \
            -overwrite
        """
        os.system(tsnr_run)