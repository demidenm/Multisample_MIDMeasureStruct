# Multisample_MIDMeasureStruct
 Project  Title: A Multi-Sample Evaluation of the Measurement Structure and Function of the Modified MID Task in Adolescents
 Project Authors: Michael I. Demidenko, Jeanette Mumford, Nilam Ram, Russell A. Poldrack
 Respository & code maintained by: Michael Demidenkoo
 
## Associated Code

This github repository contains the Python (.py) and R (.rmd/.html) that is associated with the publication titiled:
    A Multi-Sample Evaluation of the Measurement Structure and Function of the Modified Monetary Incentive Delay Task in Adolescents.

This publication is being submitted as a Registered Report. 

Stage 1 contains the proposed code for the fMRI and Measurement model.
FMRI analyses are include steps to run the Region of Interest the First and Second Level analyses. These are broken broken down into three scripts:
    
    - First Level (MID_FirstLevel.py): Uses a set of defined functions in combination with Nilearn & 3dTstat to (1) create a design matrix,
    (2) define a GLM model, (3) calculate several pre-defined contrast (beta) maps, (4) generate a cleaned BOLD file and (5) calculate a TSNR volume.
    
    - Group Level (MID_SecondLevel.py): Uses a set of defined functions in combination with Nilearn to (1) calculate a weighted fixed effects model to combined runs,
    and (2) estimates a group (i.e., second) level activation map (z-score) based on a list of inputed run beta maps.
    
    - ROI signal intensity (BetaMap_ROI.py): Uses a set of Nilearn features to use Harvord-Oxford + Neurosynth  generate ROIs to
    mask run (i.e., first level) beta maps and extract the average signal intensity across the voxels within the mask. This data is appended to a pandas 
    dataframe and saved as a .csv file for later use. In addition, the overlap between each subjects brain mask and the MNI ROI mask are used to
    generate the % of voxels from a subjects brain mask are in the MNI ROI mask. These values are saved to a pandas dataframe
    to determine which subjects have <50% voxels within each mask and need to be excluded.
    
The Measurement analyses are included within a single R markdown file. These include an R script file and R markdown (+ html file). The initial code was piloted on
simulated data. This code chunk (using simsem) is contained within Stage 1 but will not be contained within Stage 2 R files. 
    
    - R script (MID_CFA-ESEM-EFA-LSEM.r): Contains the core code to import data, define confirmatory model, run CFA model, MG-CFA, ESEM and EFA. Included comparisons
    
    - .rmd/.html (MID_CFA-ESEM-EFA-LSEM.html,MID_CFA-ESEM-EFA-LSEM.rmd): Contains the expand version of the R script. In addition to the code, the script will explain the packages and related notes for the analyses. 
    The .html file is to help provide sufficient comments and guidance for a user that is less familiar with the packages to draw some insights as to what is being
    use and why during each code block.
    
    
Any revisions (outside of typoes & coding blunders) to the Stage 1 code submitted under this repo will reflect only the reviewers requests/feedback. Thus, revisions may be made to the files within the Stage 1 folder up until in principle acceptance (IPA). 
Upon receipt of IPA from the editoral team, no further changes will be made to Stage 1. Final R code (i.e., without simulated data but fMRI input files) and relevants figs/output will be added to a Stage 2 folder contained within this repo.
