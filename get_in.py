from os import scandir, chdir, mkdir
from shutil import copy2 as cp
from sys import argv

if len(argv) > 1:
    config_dir = argv[1]
else:
    config_dir = '.'
    
ins_path = '{}\ins'.format(config_dir)
try:
    mkdir(ins_path)
except FileExistsError:
    ins_path += '1'
    mkdir(ins_path)

for d in scandir():
    if '.gid' in d.name:
        in_path = '{}\{}\{}.in'.format(config_dir, d.name, d.name.split('.gid')[0])
        cp(in_path, ins_path)

