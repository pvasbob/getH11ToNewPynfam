import pandas as pd
import os
import shutil


cutoff = 60
directories = {
     'H11_store' : 'H11_store',
     'H11_vn_store' : 'H11_vn_store',
     'NH11': 'QQ_binary_NH11',
     'PH11': 'QQ_binary_NH11',
     'output_binfo': 'QQ_binary_output_binfo',
     'opt_index' : 'opt_index_sort_rec',
     }


parent_dir = os.getcwd()

H11_store_path = os.path.join(parent_dir, directories['H11_store'])
if not os.path.exists(H11_store_path):
   raise Exception('H11_store not exists, so can not extract H11 according to v_n')

if cutoff == None:  
    raise Exception('cutoff has to be set.')


H11_vn_store_path = os.path.join(parent_dir, directories['H11_vn_store'])
if os.path.exists(H11_vn_store_path):
    H11_store_filelist = os.listdir(H11_vn_store_path)
    if not len(H11_vn_store_filelist) == 0:
        raise Exception('H11 acccording to vn may already exists, quit to check.')
else:
    os.mkdir(H11_vn_store_path)



opt_index_sort_rec_path = os.path.join(H11_store_path, 'opt_index_sort_rec')
with open(opt_index_sort_rec_path, 'r') as file:
    for line in file:
        if int(line.split(',')[-1]) > cutoff: break
        H11op_label = line.split(',')[0].split('_')[0] 
        H11ph_label = str(int(line.split(',')[0].split('_')[1]))
        print(line.split(',')[-1].strip('\n'),line.split(',')[-2],'\t', H11op_label + H11ph_label)
        NH11_new_path = os.path.join(H11_store_path, 'QQ_binary_NH11_' + H11op_label + H11ph_label)
        PH11_new_path = os.path.join(H11_store_path, 'QQ_binary_PH11_' + H11op_label + H11ph_label)
        output_binfo_new_path = os.path.join(H11_store_path, 'QQ_binary_output_binfo_' + H11op_label + H11ph_label)
        NH11_vn_path = os.path.join(H11_vn_store_path, 'QQ_binary_NH11_' + H11op_label + H11ph_label)
        PH11_vn_path = os.path.join(H11_vn_store_path, 'QQ_binary_PH11_' + H11op_label + H11ph_label)
        output_binfo_vn_path = os.path.join(H11_vn_store_path, 'QQ_binary_output_binfo_' + H11op_label + H11ph_label)
        shutil.copy(NH11_new_path, NH11_vn_path)
        shutil.copy(PH11_new_path, PH11_vn_path)
        shutil.copy(output_binfo_new_path, output_binfo_vn_path)


 

 
        

