import glob
import pandas as pd
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from nilearn.glm.first_level import make_first_level_design_matrix
from nilearn.glm.first_level import FirstLevelModel
from nilearn.plotting import plot_design_matrix
from nilearn import image
import sys

# non-interactive figure to avoid error w/ saving design matrix
plt.switch_backend('Agg')

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
    :param stc: whether to adjust the onsets/frame times in design matrix. Default True, alt False
    :return: returns a design matrix for each run with 5 anticipation,
            10 feedback, motion 6 + derivs, compcor + constant and confounds regressor file
    """

    # Adding the the events_df the cue + fixation duration for anticipation duration of model. Note, ABCD names fixation = anticiation
    events_df["ANTICIPATION_DURATION"] = events_df["Cue.Duration"] + events_df["Anticipation.Duration"]

    # concatinating the condition types from events_df, anticipation type [:,0] & feedback type [:,9]
    # concatinating the onsets from events_df for anticipation type [:,1] & feedback type [:,7]
    # concatinating the duration from events_df for anticipation type (see line 18) [:,13] & feedback type [:,8]
    conditions = pd.concat([events_df.loc[:, "Condition"], events_df.loc[:, "Result"]],
                           ignore_index=True)
    onsets = pd.concat([events_df.loc[:, "Cue.OnsetTime"], events_df.loc[:, "Feedback.OnsetTime"]],
                       ignore_index=True)
    duration = pd.concat([events_df.loc[:, "ANTICIPATION_DURATION"], events_df.loc[:, "FeedbackDuration"]],
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



# Creating list of subjects and runs to loop over
sub = sys.argv[1]

runs = ['run-01','run-02']

# Setting data (for events.tsv) and derivative (for fmriprep) paths
data_path = sys.argv[2]
deriv_path = sys.argv[3]
scratch_out = sys.argv[4]
#firstlvl_out = sys.argv[5]

# define contrasts for MID, or 'Reward' in MLS
contrasts = {
    'Lgain-Neut': 'LgReward - Triangle',
    'LSgain-Neut': '.5*LgReward + .5*SmallReward - 1*Triangle',
    'Lloss-Neut': 'LgPun - Triangle',
    'LSloss-Neut': '.5*LgPun + .5*SmallPun - 1*Triangle',
    'Lgain-Lloss': 'LgReward - LgPun',
    'Lloss-Lgain': 'LgPun - LgReward'
}

# provide TR & volumes for associated BOLD
tr = .800
vols = 403
ses='ses-2YearFollowUpYArm1'
task='MID'

for r in runs:
    print(f'Starting {r} for {sub} [Seven Steps]')

    print('1/7 Load Files & set paths')
    # import behavior events .tsv from data path, rename feedback neutral to diff Miss/Hit
    beh = pd.read_csv(f'{data_path}/{sub}/{ses}/func/{sub}_{ses}_task-{task}_{r}_events.tsv', sep='\t')
    # Neutral condition: if df['Condition'] == 'Triangle' & hit
    beh.loc[(beh['Condition'] == 'Triangle') & ((beh['prbacc'] == 1)), 'Result'] = 'Hit_No Money Stake!'
    # Neutral condition: if df['Condition'] == 'Triangle' & miss
    beh.loc[(beh['Condition'] == 'Triangle') & ((beh['prbacc'] == 0)), 'Result'] = 'Miss_No Money Stake!'
    
    # get path to confounds from fmriprep & get data in dataframe filling missing values with 0 (e.g., derivatives)
    conf = glob.glob(
        f'{deriv_path}/{sub}/{ses}/func/{sub}_{ses}_task-{task}_{r}_desc-confounds_timeseries.tsv'
    )[0]
    
    conf_df = pd.read_csv(conf, sep='\t', na_values=['n/a']).fillna(0)

    # Get path to functional data
    nii_path = glob.glob(
            f'{deriv_path}/{sub}/{ses}/func/{sub}_{ses}_task-{task}_{r}_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz'
    )[0]
    
    # Get path to functional mask
    mask_path = glob.glob(
        f'{deriv_path}/{sub}/{ses}/func/{sub}_{ses}_task-{task}_{r}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz'
    )[0]

    print('2/7 Create Regressors & Design Matrix for GLM')
    # get list of regressors
    conf_regressors = regressors(confound_df = conf_df, num_compcor=8)

    # run to create design matrix
    design_matrix = create_designmat(events_df=beh, conf_regressors=conf_regressors,
                                     hrf_model='SPM', tr=tr, volumes=vols,stc=False)

    print('3/7 Fit GLM model with 5mm FWHM smoothing, ar1 autocorrelation')
    # fitting glm for sub's run, using associated mask, using ar1 autocorrelation (FSL prewhitening), drift model
    # 'cosine' and .01 highpass (100s filter)
    fmri_glm = FirstLevelModel(subject_label={sub}, mask_img=mask_path, t_r=tr, smoothing_fwhm=5,
                               standardize=False, noise_model='ar1', drift_model=None,high_pass=None
                               # cosine 0:3 included from fmriprep in desing mat based on 128 s calc
    )

    # Run GLM model using set paths and calculate design matrix
    run_fmri_glm = fmri_glm.fit(nii_path, design_matrices=design_matrix)


    print('4/7: From GLM model, create z-score contrast maps and save to output path')
    # contrast names and associated contrasts in contrasts defined is looped over
    # contrast name is used in saving file, the contrast is used in deriving z-score
    for con_name, con in contrasts.items():
       # zscore file
       zstat_name = f'{scratch_out}/{sub}_{ses}_task-{task}_{r}_contrast_{con_name}_zstat.nii.gz'
       z_est = run_fmri_glm.compute_contrast(con, output_type='z_score')
       z_est.to_filename(zstat_name)
       # To calculate ROI mean beta, need beta maps and beta maps with variance used in fixed effect calculation
       # fixed effect calculation is used for whole brain analysis inputs.
       # Calc: beta estimate
       beta_name = f'{scratch_out}/{sub}_{ses}_task-{task}_{r}_contrast_{con_name}_beta.nii.gz'
       beta_est = run_fmri_glm.compute_contrast(con, output_type='effect_size')
       beta_est.to_filename(beta_name)
       # Calc: variance
       var_name = f'{scratch_out}/{sub}_{ses}_task-{task}_{r}_contrast_{con_name}_var.nii.gz'
       var_est = run_fmri_glm.compute_contrast(con, output_type='effect_variance')
       var_est.to_filename(var_name)


    # This creates an HTML report for the GLM specified above.
    # The contrasts are printed at a .01 alpha
    print('5/7: Using GLM model, generate design matrix with model information')
    plot_design_matrix(design_matrix=design_matrix,
                       output_file=f'{scratch_out}/{sub}_{ses}_task-{task}_{r}_model-FirstLevelDesign.png')

    # creating a filtered_bold signal, equivalent to filtered_func.nii.gz from FSL
    print('6/7: Clean BOLD filter, e.g., detrend and remove signal related confound regressors '
          '(e.g. consine, compcor, motion)')
    clean_bold = image.clean_img(nii_path, confounds=conf_regressors, t_r=tr,
                                 detrend=False, standardize=False, low_pass=None,ensure_finite=False
                                 # detrend (cosine 0:3) and regressors produced by fmriprep
    )

    # [not using. removing for space reasons] smoothing the above, cleaned bold image
    ##clean_bold_fwhm5mm = image.smooth_img(clean_bold, fwhm=5)
    # saving the filered image
    ##print('\t\t 7/7: Save filtered BOLD image to output path')
    ##nib.save(img= clean_bold_fwhm5mm,
    ##         filename=f'{scratch_out}/{sub}_{ses}_task-{task}_{r}_space-MNI152NLin2009cAsym_res-2_desc-filtered_bold.nii.gz'
    ##         )
    
    print(f'\t\t  .... sub-{sub} run {r} DONE!')

