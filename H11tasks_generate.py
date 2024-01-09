
import os
import shutil
import pandas as pd


Nucleus = 'Ce'
qrpa_lines = {
      'header' : '================qrpa.inp=============',
      'calc, mode, file io' : [1, 0],
      'O,T,L,K'  : [1, 0, 3, 0],
      'qrpa_eps' : 0.1,
      'max_iter_qrpa' : 100,
      'qrpa_nbroyden' : 40,
      'qrpa_alphamix' : 0.33,
      'line parameters' : [0.0, 0.5, 0.25, 0.0],
      'circle parameters' : [2.22, 0.0, 0.11],
      'half-circle parameters' : [0.0, 200.0],      
      'qrpa_points' : 1, 
     }



directories = {
     'lpfam_run' : 'tests/pynfam_example/000000/hfb_soln/hfb_meta/0000',
     'work_dir' : 'work_dir',
     'L_file' : 'L_file',
     'orig_lpfam_dir' : 'lpynfamH11_std+imag_Leven_width0.5',
     'qrpa_inp': 'qrpa.inp',
     'slurm' : 'slurm_longleaf.sh',
     }



slurm_lines = {
      'env' : '#!/bin/bash',
      'name' : '#SBATCH -J 1031_x.0',
      'partition' : '#SBATCH -p general',
      'ntasks' : '#SBATCH --ntasks=1',
      'time' : '#SBATCH --time=07-00:00:00',
      'output' : '#SBATCH -o out_%A.out',
      'space_line' : '    ',
      'command' : './run_pynfam.py',
      }


parent_dir = os.getcwd()

work_dir_path = os.path.join(parent_dir, directories['work_dir'])
if os.path.exists(work_dir_path):
    raise Exception('work_dir exists, quit to check if job done already.')
else:
    os.mkdir(work_dir_path)

L_file_path = os.path.join(parent_dir, directories['L_file'])
if not os.path.join(parent_dir, L_file_path):
    raise Exception('L_file not exists. quit')


orig_lpfam_path = os.path.join(parent_dir, directories['orig_lpfam_dir'])
if not os.path.exists(orig_lpfam_path):
    raise Exception('original lpynfamH11_std+imag_Leven_width0.5 NOT exist. quit')


amp_list = os.listdir(L_file_path)
if len(amp_list) == 0:
    raise Exception('NO 46_amplitude files in L_file, quit.')

for amp in amp_list:
    # get last 3 numbers for T, L, K
    task = amp.split('_')[-1] 
    if len(list(task)) != 3: 
        raise Exception('Invalide name for 46_amplitude file, quit to check')
    T_op = list(task)[0]
    L_op = list(task)[1]
    K_op = list(task)[2]

    # get lists of local maximum from this amp
    locamax_ind = {}
    amp_path = os.path.join(L_file_path, amp)
    data = pd.read_csv(amp_path, header = None, delimiter = ',')  # header = None  important
    omega = data[0]
    stren = data[3] 
    globmax = max(stren) 
    print(amp)
    ct = 1
    for i in list(range(1, len(stren)-1)):
       #if stren[i] > stren[i-1] and stren[i] > stren[i+1]:        
        if stren[i] > stren[i-1] and stren[i] > stren[i+1] and (stren[i] > globmax/10 or omega[i] < 10.0) :        
            locamax_ind[ct] = omega[i] 
            ct = ct +1

    for ind in locamax_ind:
        task_ind_path = os.path.join(work_dir_path, task + '_' + '{0:03d}'.format(ind))
        shutil.copytree(orig_lpfam_path, task_ind_path)

        task_ind_0000_path = os.path.join(task_ind_path, directories['lpfam_run'])
        qrpa_inp_path = os.path.join(task_ind_0000_path, directories['qrpa_inp'])
        if os.path.exists(qrpa_inp_path): os.remove(qrpa_inp_path)

        qrpa_lines['line parameters'][0] = locamax_ind[ind]
        qrpa_lines['O,T,L,K'][1] = int(T_op) 
        qrpa_lines['O,T,L,K'][2] = int(L_op) 
        qrpa_lines['O,T,L,K'][3] = int(K_op) 

        with open(qrpa_inp_path, "a") as file:
            file.write(qrpa_lines['header'] + '\n') 
            file.write(str(qrpa_lines['calc, mode, file io'])[1:-1] + '\n')
            file.write(str(qrpa_lines['O,T,L,K'])[1:-1] + '\n')
            file.write(str(qrpa_lines['qrpa_eps']) + '\n')
            file.write(str(qrpa_lines['max_iter_qrpa']) + '\n')
            file.write(str(qrpa_lines['qrpa_nbroyden']) + '\n')
            file.write(str(qrpa_lines['qrpa_alphamix']) + '\n')
            file.write(str(qrpa_lines['line parameters'])[1:-1] + '\n')
            file.write(str(qrpa_lines['circle parameters'])[1:-1] + '\n')
            file.write(str(qrpa_lines['half-circle parameters'])[1:-1] + '\n')
            file.write(str(qrpa_lines['qrpa_points']) + '\n')
 

        slurm_path = os.path.join(task_ind_path, directories['slurm'])
        if os.path.exists(slurm_path) : os.remove(slurm_path)
        slurm_lines['name'] = '#SBATCH -J 1031_x.0'
        slurm_lines['name'] = '#SBATCH -J ' + Nucleus  \
                          + str(qrpa_lines['O,T,L,K'][1]) \
                          + str(qrpa_lines['O,T,L,K'][2]) \
                          + str(qrpa_lines['O,T,L,K'][3]) \
                          + str('{0:02d}'.format(int(ind)))



        with open(slurm_path, 'a') as file:
            file.write(slurm_lines['env'] + '\n')
            file.write(slurm_lines['name'] + '\n')
            file.write(slurm_lines['partition'] + '\n')
            file.write(slurm_lines['ntasks'] + '\n')
            file.write(slurm_lines['time'] + '\n')
            file.write(slurm_lines['output'] + '\n')
            file.write(slurm_lines['space_line'] + '\n')
            file.write(slurm_lines['command'] + '\n')

        print(f'job {task + str(ind)} locate at {locamax_ind[ind]} created')


all_tasklist = os.listdir(work_dir_path)
for task_ind in all_tasklist:
    task_ind_path = os.path.join(work_dir_path, task_ind)
    os.chdir(task_ind_path)
    os.system('sbatch slurm_longleaf.sh') 


print('job assignment on longleaf done')
