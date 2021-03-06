from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast

parser= argparse.ArgumentParser()

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-MIND1","--mind-atom", default = '-3.5', help = "input minimum distance between atoms", type = str)
parser.add_argument("-MIND2","--mind-symm", default = '2.2',help = "input minimum distance between symmetry", type = str)
parser.add_argument("-resl","--resolution", default = '3.0', help = "input resolution", type = str)
parser.add_argument("-ESEL", "--minimum-e", default = '1.2', help = 'input Minimum E', type = str)
parser.add_argument("-TEST", "--test-min-del", nargs='+', help = 'after find define the new starting atom')
parser.add_argument("-NTRY", "--number-of-try", default = 1000, help = 'enter the number of try', type = int)
parser.add_argument("-FIND", "--number-of-atoms", default = 4, help='enter the number of atoms', type = str)
parser.add_argument("-SFAC", "--type-of-atoms", help='enter the type of atom', type = str)

parser.add_argument("-thre","--likelihood-threshold",help = "setup likelihood threshold", type = float)

args = parser.parse_args()

if args.reflection_mtz:
    reflection_file = args.reflection_mtz

min_dist_atom = '-3.5'

min_dist_symm = '2.2'

NOATOM = 4
NOTRY = 1000


if args.number_of_try:
    NOTRY = args.number_of_try

if args.number_of_atoms:
    NOATOM = args.number_of_atoms



################################ Deploying SHELXC AND GNENERATE FILES FOR SHELXD ###################################################

# Input xxx_xxx_xx_x.mtz 
# Generate xxx_xxx_xx_x.sca 
# Output: 1) xxxxxxxxx.hkl 
#         2) xxxxxxxxx_fa.hkl
#         3) xxxxxxxxx_fa.ins
#         4) xxxxxxxxx_sad.cif
# The whole process is idempotent

########################################## Extract MTZ Information#################################
process = subprocess.Popen('phenix.mtz.dump '+reflection_file, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()

split_out=out.splitlines()

mylist = []
for i in range(len(split_out)):
    item = split_out[i].decode("utf-8")
    if item != '':
        mylist.append(split_out[i].decode("utf-8"))

my_mtz={}
for each_line in mylist:
    #each_line=each_line.replace(' ','')
    if len(each_line.split(":"))== 2:
        key, value = each_line.split(":")[0],each_line.split(":")[1]
        my_mtz[key] = value

#find the unit cell and symmetry parameter
unit_cell_search = 'unit cell'

for key in my_mtz:
      if unit_cell_search in key.lower():
            my_mtz[key]=my_mtz[key].replace('(','').replace(')','').replace(',',' ')
            CELL= my_mtz[key]
        #print('unit_cell='+my_mtz[key])


space_group_search = 'space group symbol'

for key in my_mtz:
      if space_group_search in key.lower():
            
            SPACE_GROUP= my_mtz[key]
        #print('space_group='+my_mtz[key])

print(CELL)

print(SPACE_GROUP)
#############################################################################################################

########################### Convert mtz to sca ################################

if os.path.isfile("shelxc.inp"):
    os.remove("shelxc.inp")

os.system('mtz2sca '+reflection_file)

reflection_file=reflection_file.replace('.mtz','.sca')

molecule_inp = {'CELL':CELL, 'SPAG':SPACE_GROUP, 'FIND': NOATOM, 'NTRY':NOTRY, 'SFAC':args.type_of_atoms, 'SAD':reflection_file}

for i in molecule_inp:
    my_string = i + ' ' +json.dumps(molecule_inp[i])
    my_string = my_string.replace('"','')
    print(my_string,file=open("shelxc.inp", "a"))

title_name = reflection_file.replace('.sca','')#.replace('_','')


##print(job_name)

#print('shelxc '+job_name+' < shelxc.inp')

os.system('shelxc '+title_name+' < shelxc.inp')

##################################################################################
###### Modify ins file######
if title_name+'_fa.ins' in os.listdir(os.getcwd()):
    content=open(title_name+'_fa.ins', "r").read()

    ins_list = []
    for i in content.split('\n'):
        ins_list.append(i)
else:
    print('.ins file has not been generated')

for i in range(len(ins_list)):
    #Change MIND
    if 'MIND' in ins_list[i]:
        dist_atoms = ins_list[i].split(' ')[1]
        dist_symm = ins_list[i].split(' ')[2]
        ins_list[i] = ins_list[i].replace(dist_atoms, args.mind_atom).replace(dist_symm, args.mind_symm)
        
    #Change RESOLUTION
    if 'SHEL' in ins_list[i]:
        resolution_origin = ins_list[i].split(' ')[2]
        ins_list[i] = ins_list[i].replace(resolution_origin, args.resolution)

    elif 'END' in ins_list[i]:
        for j in range (i, len(ins_list)):
            ins_list[j] = ' '

# for i in ins_list:
#     if i == ' ':
#         del(ins_list[ins_list.index(i)])

if args.minimum_e:
    ins_list.append('ESEL '+args.minimum_e)

if args.test_min_del:
    # test_min_del_list = ast.literal_eval(args.test_min_del)
    ins_list.append('TEST '+str(args.test_min_del[0])+' '+str(args.test_min_del[1]))


ins_list.append('END')

for j in range(2):
    for i in ins_list:
        if i == ' ':
            del(ins_list[ins_list.index(i)])

os.remove(title_name+"_fa.ins")
for i in ins_list:
    print(i, file=open(title_name+"_fa.ins", "a"))

#tell users the change
print('Modified '+title_name+'_fa.ins file according to the input parameter:')
for i in ins_list:
    print(i)

#############################################################################################################

#################################################### DEPLOYING SHELXD##############################################
ins_file = title_name+'_fa.ins'
hkl_file = title_name+'_fa.hkl'
os.system('shelxd ' + title_name + '_fa')

#################################################### modify pdb file ######################################
if title_name+'_fa.pdb' in os.listdir(os.getcwd()):
    pdbFileIn = title_name+'_fa.pdb'
else:
    print('pdb file has not been correctly generated')

pdbFileOut = title_name+'_fa_cleaned.pdb'

THRESHOLD = 0.3

if args.likelihood_threshold:
    THRESHOLD = args.likelihood_threshold

my_pdb = []
with open(pdbFileIn,'r') as f:
    for line in f:
        my_pdb.append(line.rstrip('\n'))


exp_num_atoms = 0
for line in range(len(my_pdb)):
    if 'HETATM' in my_pdb[line]:
        likelihood = float(my_pdb[line].split()[8])
	exp_num_atoms += 1
        if ' S ' in my_pdb[line]:
            my_pdb[line] = my_pdb[line].replace(' S ','SE')
        if likelihood < THRESHOLD:
            my_pdb[line] = 'to_be_deleted'  
	    exp_num_atoms -= 1
my_pdb.remove('to_be_deleted')

print('exp_num_atoms = '+str(exp_num_atoms),file = open("Guessed_atom_number.txt","a"))
if os.path.isfile(pdbFileOut):
    os.remove(pdbFileOut)

for pdb_line in my_pdb:
    print(pdb_line,file=open(pdbFileOut, "a"))
