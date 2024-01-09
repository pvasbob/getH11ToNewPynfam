
import os

directories = {
        'all_jobs_to_cancel' : 'all_jobs_to_cancel.txt'
        }


os.system("squeue | grep qunqun > all_jobs_to_cancel.txt")


parent_dir = os.getcwd()
jobs_cancel_path = os.path.join(parent_dir, directories['all_jobs_to_cancel'])

contents = []

with open(jobs_cancel_path, 'r') as file:
    f = file.readlines()
    lines = [line.rstrip() for line in f]


job_num = []
for line in lines:
    job_num.append(line.rstrip().split()[0])


for num in job_num:
    os.system(f'scancel {num}') 


