import pandas as pd
import numpy as np
from nilearn import image
from nilearn.masking import apply_mask


# Set input and output paths
data_dir = './Pilot_N3'
ROIs_dir = f'{data_dir}/ROIs'
fmri_prep = f'{data_dir}/derivatives/fmri_prep'
fstlvl_dir = f'{data_dir}/derivatives/fmri_prep/firstlvl_out/contrast_files'
df_out = f'{data_dir}/data'

# Load in binarized ROIs (these are created using the provided .txt file)
# L/R NAcc is 50thr Harvard-Oxford Subcortical mask
# L/R Ant Insula is the 50thr Harvard-Oxford Cortical mask that is masked by anterior masked downloaded from neurosynth (thresh 8)
#   Thus, remaining part of the Harvard-Oxford Insula is only the anterior portion overlapping with neurosynth mask
r_nacc = image.load_img(f'{ROIs_dir}/Left_NAcc.nii.gz')
L_nacc = image.load_img(f'{ROIs_dir}/Right_NAcc.nii.gz')
r_ains = image.load_img(f'{ROIs_dir}/Right_AnteriorInsula.nii.gz')
l_ains = image.load_img(f'{ROIs_dir}/Left_AnteriorInsula.nii.gz')

# List of Subjects and Runs
subj = ["31","21"]
runs = ["01", "02"]

# Contrasts (beta maps) to work through;
# Differentiated by uppercase, LS = Large+Small; L = Large;
# loss = 'Don't Lose', gain = 'Win', Neut = 'No Money at Stake'
contrasts = ["Lgain-Neut", "LSgain-Neut", "Lgain-Lloss",
             "Lloss-Neut", "LSloss-Neut", "Lloss-Lgain"]

# setting up ROI label and variable name for nifti file
rois = {"L-Ins": l_ains,
        "R-Ins": r_ains,
        "L-NAcc": l_nacc,
        "R-NAcc": r_nacc}


# Creating empty lists where the mean beta values and item names are saved
beta_vals = []
beta_items = []

# looping and extract mean beta estimate by masking img (beta map) by mask (roi).
# the '*beta.nii.gz' image is the 'effect_size' output file from nilearn's firstlevel compute_contrast()
for s in subj:
    for r in runs:
        for c in contrasts:
            for roi,img in rois.items():
                beta_map = image.load_img(f'{fstlvl_dir}/sub-{s}_ses-01_task-mid_run-{r}_contrast_{c}_beta.nii.gz')
                beta_mean = np.mean(apply_mask(imgs=beta_map, mask_img=img))
                beta_vals.append(beta_mean)
                beta_items.append(f'sub-{s}_{r}_{c}_{roi}')

# After loop which creates a row for each subject, run, contrast, roi combination, reshape the data into something that is wide
# wide formate conversation: subjects remain as rows. Columns are run_contrast_roi with associated row values representing beta estimate for contrasts
df = pd.DataFrame({'Items': beta_items,
              'Values': beta_vals})
df[['subj','run', 'contrast','roi']] = df.Items.str.split("_",expand=True)
df = df.drop(columns=['Items'])
df = df.pivot_table(index="subj",
                    columns=["run","contrast","roi"],
                    values ="Values")

# Since the variables are in stacked format, loop through x,y,z (3 stacked variables), to combined and produce {run_contrast_roi} column name
df.columns = [f'{x}_{y}_{z}' for x,y,z in df.columns]
# export signal intensity df into .csv
df.to_csv(f'{df_out}/ROImeansignal.csv')

# After calculating ROI estimates, calculating % signal for each ROI (to exclude for excess dropout)
roi_vals = []
roi_items = []
for s in subj:
    for r in runs:
            for roiname,roi_img in rois.items():
                brain = image.load_img(f'{fmri_prep}/sub-{s}/ses-1/func/sub-{s}_ses-1_task-mid_run-{r}_*brain_mask.nii.gz')
                brain_roimask = apply_mask(imgs=brain, mask_img=roi_img)
                roi = roi_img.get_fdata()
                brain_vox = brain_roimask.sum()
                roi_vox = roi.sum()
                BrRo_perc = brain_vox/roi_vox
                roi_vals.append(BrRo_perc)
                roi_items.append(f'sub-{s}_{r}_{roiname}')

# Do the same as for mean signal intensity df, but for the voxel overlap
ROIvox_df = pd.DataFrame({'Items': roi_items,
              'Values': roi_vals})
ROIvox_df[['subj','run', 'roiname']] = ROIvox_df.Items.str.split("_",expand=True)
ROIvox_df = ROIvox_df.drop(columns=['Items'])
ROIvox_df = ROIvox_df.pivot_table(index="subj",
                                  columns=["run","roiname"],
                                  values ="Values")

ROIvox_df.columns = [f'{x}_{y}' for x,y in ROIvox_df.columns]
# export voxel overlap df into .csv
ROIvox_df.to_csv(f'{df_out}/ROIvMaskVoxels.csv')
