import sys
import os
import pandas as pd
import numpy as np
from nilearn import image
from nilearn.masking import apply_mask


# input/output paths
fmri_prep = sys.argv[1]
ROIs_dir = sys.argv[2]
first_input = sys.argv[3] # contrast files location
roi_out = sys.argv[4] # fixed effect output location

# List of Subjects and Runs
subj = os.listdir(first_input)

runs = ["01", "02"]
task = 'MID'
ses = '2YearFollowUpYArm1' 

#contrast list
contrasts = ["Lgain-Neut", "LSgain-Neut", "Lgain-Lloss",
             "Lloss-Neut", "LSloss-Neut", "Lloss-Lgain"]


# Load in binarized ROIs (these are created using the provided .txt file)
# L/R NAcc is 50thr Harvard-Oxford Subcortical mask
# L/R Ant Insula is the 50thr Harvard-Oxford Cortical mask that is masked by anterior masked downloaded from neurosynth (thresh 8)
#   Thus, remaining part of the Harvard-Oxford Insula is only the anterior portion overlapping with neurosynth mask
l_nacc = image.load_img(f'{ROIs_dir}/Left_NAcc.nii.gz')
r_nacc = image.load_img(f'{ROIs_dir}/Right_NAcc.nii.gz')
r_ains = image.load_img(f'{ROIs_dir}/Right_AnteriorInsula.nii.gz')
l_ains = image.load_img(f'{ROIs_dir}/Left_AnteriorInsula.nii.gz')


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
print("1. Mean ROI Running for N subjects: ", len(subj))

for sub in subj:
    #print(f'\tStart {sub}')
    for r in runs:
         for c in contrasts:
             for roi,img in rois.items():
                 beta_map = image.load_img(f'{first_input}/{sub}/{sub}_ses-{ses}_task-{task}_run-{r}_contrast_{c}_beta.nii.gz')
                 beta_mean = np.mean(apply_mask(imgs=beta_map, mask_img=img))
                 beta_vals.append(beta_mean)
                 beta_items.append(f'{sub}_{r}_{c}_{roi}')

# After loop which creates a row for each subject, run, contrast, roi combination, reshape the data into something that is wide
# wide formate conversation: subjects remain as rows. Columns are run_contrast_roi with associated row values representing beta estimate for contrasts
print('2. Convert beta values and item vectors into a data.frame')
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
df.to_csv(f'{roi_out}/region_roi-meansignal.csv')



# After calculating ROI estimates, calculating % signal for each ROI (to exclude for excess dropout)
roi_vals = []
roi_items = []

print("3.  Voxel in ROI Running for: ", len(subj))
for sub in subj:
    #print(f'\tStart {sub}')
    for r in runs:
        for roiname,roi_img in rois.items():
            brain = image.load_img(f'{fmri_prep}/{sub}/ses-{ses}/func/{sub}_ses-{ses}_task-{task}_run-{r}_*brain_mask.nii.gz')
            brain_roimask = apply_mask(imgs=brain, mask_img=roi_img)
            roi = roi_img.get_fdata()
            brain_vox = brain_roimask.sum()
            roi_vox = roi.sum()
            brainroi_perc = brain_vox/roi_vox
            roi_vals.append(brainroi_perc)
            roi_items.append(f'{sub}_{r}_{roiname}')

# Do the same as for mean signal intensity df, but for the voxel overlap
print('4. Convert % voxel overlap values and item vectors into a data.frame')
roi_vox_df = pd.DataFrame({'Items': roi_items,
              'Values': roi_vals})
roi_vox_df[['subj','run', 'roiname']] = roi_vox_df.Items.str.split("_",expand=True)
roi_vox_df = roi_vox_df.drop(columns=['Items'])
roi_vox_df = roi_vox_df.pivot_table(index="subj",
                                  columns=["run","roiname"],
                                  values ="Values")

roi_vox_df.columns = [f'{x}_{y}' for x,y in roi_vox_df.columns]
# export voxel overlap df into .csv
roi_vox_df.to_csv(f'{roi_out}/region_overlay-roivoxels.csv')
