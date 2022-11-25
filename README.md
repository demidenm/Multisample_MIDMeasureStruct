# Multisample_MIDMeasureStruct
Project  Title: _A Multi-Sample Evaluation of the Measurement Structure and Function of the Modified MID Task in Adolescents_
Project Authors: Michael I. Demidenko, Jeanette A. Mumford, Nilam Ram, Russell A. Poldrack

> Respository & code maintained by: Michael Demidenko
 
 
This github repository contains the Python (.py) and R scripts (.rmd/.html) that are associated with the project. This project is being submitted as a Registered Report. As described by the [OSF](https://osf.io/rr/), there are two stages to the registered report. Stage 1 is the submission and review of introduction, methods and proposed analyses for the publication. The Stage 2 is the completion of the analyses (as proposed), results and discussion section.

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

## Stage 2 Associated Code
