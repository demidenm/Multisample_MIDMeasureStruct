[![Funded By](https://img.shields.io/badge/NIDA-F32%20DA055334--01A1-yellowgreen?style=plastic)](https://reporter.nih.gov/project-details/10525501)
[![RegReport](https://img.shields.io/badge/Stage_1-Registered_Report-red
)](https://osf.io/a6t8k)
[![DOI](https://img.shields.io/badge/DOI-Publication-blue
)](https://osf.io/a6t8k)

# Multisample_MIDMeasureStruct
Project  Title: 

    A Multi-Sample Evaluation of the Measurement Structure and Function of the Modified MID Task in Adolescents_

Project Authors: 
    
    Michael I. Demidenko, Jeanette A. Mumford, Nilam Ram, Russell A. Poldrack

> Respository & code maintained by: Michael Demidenko
 
 
This github repository contains the Python (.py) and R scripts (.rmd/.html) that are associated with the project. \
This project is being submitted as a Registered Report. As described by the [OSF](https://osf.io/a6t8k), there are two stages to the \
registered report. Stage 1 is the submission and review of introduction, methods and proposed analyses for the publication.\
The Stage 2 is the completion of the analyses (as proposed), results and discussion section.

## Stage 2 Associated Code
Stage 2 contains the final analysis code for the fMRI and Measurement model. In the Stage2_Code folder, there are four subfolders with their own subfolders:
1. ABCD
   1. code_midinvar
      1. ROI
      2. scripts
         1. folders: python_scripts, templates, qc_out
            1. python_script .py files: betamap_roi, mid_firstlevel, mid_fixedeff, mid_group, mid_group runs, tstat_to_cohensd
         2. .sh files: FirstLevel, Group, Group_runs, make_first, roi_batch, view_imgs
   2. config
2. AHRB
   1. code_midinvar
         2. python_scripts
            1. folders: python_scripts, qc_out
               1. python_script .py files: betamap_roi, mid_firstlevel, mid_fixedeff, mid_group, mid_group runs, tstat_to_cohensd
            2. .sh files: FirstLevel, Group, Group_runs, roi_batch
   2. config
3. MLS
   1. code_midinvar
      1. python_scripts
         1. folders: python_scripts, qc_out
            1. python_script .py files: betamap_roi, mid_firstlevel, mid_fixedeff, mid_group, mid_group runs, tstat_to_cohensd
         2. .sh files: FirstLevel, Group, Group_runs, roi_batch
   2. config
4. preprocessng templates

First, the AHRB and MLS analyses are run on Stanford University's Sherlock high-performance computers (HPC) and the ABCD data are run on the \
the University of Minnesota's MSI HPCs. Briefly, each folder (AHRB, MLS and ABCD) has a template script used to run the `fmriprep v23.1.0` preprocessing. \
Then, each folder contains a subfolder code_midinvar with the relevant scripts and ROIs. Note, ROIs are redundant, so that folder is only provided for ABCD. \
the *scripts/ folder contains a collection of python scripts (*scripts/python_scripts/) and bash files. 

**Bash files** (.sh): used to submit jobs to the HPCs for the first level (e.g., FirstLevel.sh), group level \
for average runs (Group.sh) and per run (Group_runs.sh), and the ROI analyses (roi_bath.sh).

**Python scripts** (.py): There are seven separate python scripts used in these analyses. They will be nearly identical for each of the three samples. The scripts and the details briefly:

- betamap_roi.py: Used in ABCD, AHRB and MLS. This script iterates over subjects, runs, and contrasts, extracting the mean beta values by applying the [ROI](./Stage2_Code/ABCD/code_midinvar/ROI) masks to the contrast maps. The results are then structured into a DataFrame that includes beta values, along with subject, run, contrast, and ROI information. Finally, the script calculates the percentage of voxel overlap between the ROI masks and functional brain masks. This information is also formatted and exported into CSV files. 
- extract_values.py: Used in ABCD only. Summarize and compile data related to subjects' performance for each MID task run. Reads data from two TSV files representing runs 1 and 2, calculating the mean framewise displacement for each run + extracts accuracy and mean reaction time (RT) data from a JSON file for each run. The metrics are then appended to a DataFrame and export to CSV.
- mid_firstlevel.py: Used in ABCD, AHRB and MLS. Performs a series of data preprocessing and analysis steps using Nilearn and other libraries. Creates confound regressors (e.g., motion + CompCor) and design matrix include MID task regressors (15) + confounds. Fit GLM for specific design matrix w/ temporal autocorr + prewhitening. Estimate contrast, saving z-scores, beta maps and variances for each. 
- mid_fixedeff.py: Used in ABCD, AHRB and MLS. Compute fixed effects across six contrast for each subjects' MID fMRI data. Uses beta maps and variance maps for each of the six contrasts using Nilearn's compute_fixed_effects function. Resulting fixed effect, variance, and t-statistic maps are saved into the an output directory.
- mid_group.py: Used in ABCD, AHRB and MLS. Performs a one-sample T-test at the group level for a list of six contrasts. It calculates the group mean effect by fitting an intercept to the length of the input fixed-effect maps using Nilearn's SecondLevelModel. Computes z- and t-stats maps for each of the six contrast.
- mid_group_runs.py: Used in ABCD, AHRB and MLS. Same as `mid_group.py` but instead uses run inputs produced by `mid_firstlevel.py`.
- tstat_to_cohensd.py: Used in ABCD, AHRB and MLS. Takes T-statistic maps and calculates Cohen's d effect size based on these statistics and the degrees of freedom (n).

## Stage 1 Associated Code

Stage 1 contains the proposed code for the fMRI and Measurement model.
FMRI analyses are include steps to run the Region of Interest the First and Second Level analyses. These are broken broken down into three scripts:
    
 - First Level (MID_FirstLevel.py): Uses a set of defined functions in combination with Nilearn & 3dTstat to (1) create a design matrix,
    (2) define a GLM model, (3) calculate several pre-defined contrast (beta) maps, (4) generate a cleaned BOLD file and (5) calculate a temporal SNR volume.
 - Group Level (MID_SecondLevel.py): Uses a set of defined functions in combination with Nilearn to (1) calculate a weighted fixed effects model to combine runs, and (2) estimates a group (i.e., second) level activation map (z-score) based on a list of inputed run beta maps.
 - ROI signal intensity (BetaMap_ROI.py): Uses a set of Nilearn functions and Harvord-Oxford + Neurosynth generated ROIs to
    mask run (i.e., first level) beta maps and extract the average signal intensity across the voxels within the specified mask. This data is appended to a pandas 
    dataframe and saved as a .csv file for later use. In addition, the overlap between each subject's brain mask and the MNI ROI mask are used to
    generate the % of voxels from a subjects brain mask are in the MNI ROI mask. These values are saved to a pandas dataframe
    to determine which subjects have <50% voxels within each mask and need to be excluded.
    
The Measurement analyses are performed in R. Files include: R script and R markdown files (+ html file). The initial code was piloted on simulated data. This code chunk (using simsem) is contained within Stage 1 but will not be contained within Stage 2 R files. 
    
 1. R script (MID_CFA-ESEM-EFA-LSEM.r): Contains the core code to import data, define confirmatory model, run CFA model, MG-CFA, ESEM and EFA. Included comparisons 
 2. .html/.rmd (MID_CFA-ESEM-EFA-LSEM.html,.rmd,.pdf): Contains the expand version of the R script. In addition to the code, the script will explain the packages and related notes for the analyses. The .html file is generated to provide sufficient comments and guidance for a user that is less familiar with the packages to draw some insights as to why, what and how is being run during each code block. The associated the .rmd file is provided to reproduce the .html file.
    
    
Revisions (outside of typos & identified coding errors) to the Stage 1 code submitted under this repo will reflect only the reviewers requests/feedback. Thus, revisions may be made to the files within the Stage 1 folder up until in principle acceptance (IPA). Upon receipt of IPA from the editoral team, no further changes will be made to Stage 1. Final R code (i.e., without simulated data but fMRI input files) and relevants figs/output will be added to a Stage 2 folder contained within this repo and details expanded below.

