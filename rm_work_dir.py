
import os
import shutil 

parent_dir = os.getcwd()

work_dir = os.path.join(parent_dir, 'work_dir')
shutil.rmtree(work_dir)

