import pandas as pd
import json


# List of subjects and sessions

deriv_dir="/oak/../AHRB/derivatives"
fp_dir = f'{deriv_dir}/fmriprep_23.1.0'
task = 'mid'


# read in list of IDs/ses to summarize
id_ses_df = pd.read_csv('id_ses.tsv', delimiter='\t', header=None)

dfs = []
count_n = 0

# Loop through subjects and sessions
for index, row in id_ses_df.iterrows():
    subject = row.iloc[0]
    ses = row.iloc[1]
    # assign dirs to pull from
    file_path_run1 = f'{fp_dir}/{subject}/ses-{ses}/func/{subject}_ses-{ses}_task-{task}_run-01_desc-confounds_timeseries.tsv'
    file_path_run2 = f'{fp_dir}/{subject}/ses-{ses}/func/{subject}_ses-{ses}_task-{task}_run-02_desc-confounds_timeseries.tsv'
    
    # Read the TSV files into pandas DataFrames
    df_run1 = pd.read_csv(file_path_run1, delimiter='\t')
    df_run2 = pd.read_csv(file_path_run2, delimiter='\t')
    
    # Calculate the mean framewise displacement for run-01/run-02
    mean_fd_run1 = df_run1['framewise_displacement'].mean()
    mean_fd_run2 = df_run2['framewise_displacement'].mean()
    
    # get json info
    j_file = f'{deriv_dir}/taskdescribe_v1.0/{ses}_{task}/{subject}_ses-{ses}_task-{task}_beh-descr.json'
    with open(j_file, 'r') as file:
        data = json.load(file)
    run1_accuracy = data['Run 1']['Overall Accuracy']
    run1_mean_rt = data['Run 1']['Mean RT']
    run2_accuracy = data['Run 2']['Overall Accuracy']
    run2_mean_rt = data['Run 2']['Mean RT']

    
    # Append the results to the DataFrame
    result_df = pd.DataFrame({
        'Subject': [subject],
        'Session': [ses],
        'mFD_run1': [mean_fd_run1],
        'mFD_run2': [mean_fd_run2],
        'acc_run1': [run1_accuracy],
        'acc_run2': [run2_accuracy],
        'mrt_run1': [run1_mean_rt],
        'mrt_run2': [run2_mean_rt]
        })
    
    dfs.append(result_df)
    
    count_n += 1

# Write out the DataFrame to a CSV file
result_df = pd.concat(dfs, ignore_index=True)
result_df.to_csv(f'{deriv_dir}/taskdescribe_v1.0/subs-{count_n}_task-{task}_summ-mot-acc-rt.csv', index=False)
