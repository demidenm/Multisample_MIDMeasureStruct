import numpy as np
from nilearn.image import load_img, new_img_like
import sys
import glob
import os

def t_to_cohensd(t_stat, n):
    # calc cohens d as t / sqrt(n)
    cohens_d = t_stat / np.sqrt(n)
    return cohens_d

# in path
inp_dir = sys.argv[1]
img_list = glob.glob(f"{inp_dir}/**/*stat-tstat.nii.gz", recursive=True)

for img_path in img_list :  
    # load data & get data array (fdata)
    f_img = img_path
    img = load_img(f_img)
    data = img.get_fdata()
    
    # While z / sqrt(n) may be same value, to be exact, calculate z-to-p, p-to-t and t-cohens
    n = int(f_img.split('subs-')[1].split('_')[0])
    contrast = f_img.split('contrast')[1].split('_')[0]
    
    # run function
    d_data = t_to_cohensd(data, n)
    
    # create nifti image same as original
    cohens_img = new_img_like(img, d_data)
    
    # Save the result as a new NIfTI file
    new_name = img_path.replace("stat-tstat", "stat-cohensd")
    cohens_img.to_filename(new_name)
