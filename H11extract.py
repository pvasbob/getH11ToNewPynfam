import pandas as pd
import os
import shutil


directories = {
     'H11' : 'tests/pynfam_example/000000/hfb_soln/',
     'H11_store' : 'H11_store',
     'work_dir' : 'work_dir',
     'NH11': 'QQ_binary_NH11_X',
     'PH11': 'QQ_binary_NH11_X',
     'output_binfo': 'QQ_binary_output_binfo',
     'opt' : '47_optimize',
     }


parent_dir = os.getcwd()

H11_store_path = os.path.join(parent_dir, directories['H11_store'])
if os.path.exists(H11_store_path):
    H11_store_filelist = os.listdir(H11_store_path)
    if not len(H11_store_filelist) == 0:
        raise Exception(' QQ_binary_NH11 files may already exists in H11_store, quit to check in H11_store')
else:
    os.mkdir(H11_store_path)


work_dir_path = os.path.join(parent_dir, directories['work_dir'])
#copy the odd ones into the even one working_dir.
#os.system('mv ../odd/work_dir/*   ./work_dir')

fintask_list =os.listdir(work_dir_path)
if len(fintask_list) == 0: 
    raise Exception('No finished task, quit')


# check if QQ_binary_NH11 ... in all finished tasks present:
for fintask in fintask_list:
    fintask_path = os.path.join(work_dir_path, fintask)
    H11_path = os.path.join(fintask_path, directories['H11'])
    H11_path_filelist = os.listdir(H11_path)
    if not ((directories['NH11'] in H11_path_filelist) and  \
            (directories['PH11'] in H11_path_filelist) and \
            (directories['output_binfo'] in H11_path_filelist) and 
            (directories['opt'] in H11_path_filelist)): 
        raise Exception(f'H11 files in {fintask} missed, quit')


opt_rec_path = os.path.join(H11_store_path, 'opt_rec')
opt_rec_file=open(opt_rec_path, 'w')
 # copy all Q_binary_NH11 to H11_store directory   
for fintask in fintask_list:
    print(fintask)
    H11op_label = fintask.split('_')[-2]
    H11ph_label = str(int(fintask.split('_')[-1]))
    print(H11op_label, H11ph_label)
    fintask_path = os.path.join(work_dir_path, fintask)
    H11_path = os.path.join(fintask_path, directories['H11'])
    NH11_path = os.path.join(H11_path, directories['NH11'])
    PH11_path = os.path.join(H11_path, directories['PH11'])
    output_binfo_path = os.path.join(H11_path, directories['output_binfo'])
    NH11_new_path = os.path.join(H11_store_path, 'QQ_binary_NH11_' + H11op_label + H11ph_label)
    PH11_new_path = os.path.join(H11_store_path, 'QQ_binary_PH11_' + H11op_label + H11ph_label)
    output_binfo_new_path = os.path.join(H11_store_path, 'QQ_binary_output_binfo_' + H11op_label + H11ph_label)
    output_binfo_new_path = os.path.join(H11_store_path, 'QQ_binary_output_binfo_' + H11op_label + H11ph_label)
    shutil.copy(NH11_path, NH11_new_path)
    shutil.copy(PH11_path, PH11_new_path)
    shutil.copy(output_binfo_path, output_binfo_new_path)

    
    opt_path = os.path.join(H11_path, directories['opt'])
    with open(opt_path) as opt_lines:
        for line in opt_lines:
            line = fintask + ',' + line 
            opt_rec_file.write(line)

opt_rec_file.close()


#opt_rec_file = pd.read_csv(opt_rec_path, sep=',',names=["TLKN", "Omega", "Strength", "v_n"])
#opt_sort_rec_file = opt_rec_file.sort_values(by="TLKN")
opt_rec_file = pd.read_csv(opt_rec_path, sep=',', header=None)
opt_sort_rec_file = opt_rec_file.sort_values(3, ascending=False)

opt_sort_rec_path = os.path.join(H11_store_path, 'opt_sort_rec')
opt_sort_rec_file.to_csv(opt_sort_rec_path, index=False, header= False)


opt_index_sort_rec_path = os.path.join(H11_store_path, 'opt_index_sort_rec')
opt_index_sort_rec_file = open(opt_index_sort_rec_path, 'w')

with open(opt_sort_rec_path, 'r') as lines:
    index = 0
    for line in lines:
        index = index +1
        line = line.strip('\n') + ','+ str(index) 
        print(line)
        opt_index_sort_rec_file.write(line)
        opt_index_sort_rec_file.write('\n')

opt_index_sort_rec_file.close()


 


   









