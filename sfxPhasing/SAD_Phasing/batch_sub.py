# This is the top file of the SAD pipeline. 
# It is used to carry out the batch submission
# This .py will do several automation and grid definition, and parse
# these information to SAD-automation file.
#
from __future__ import print_function
import sys
import subprocess
import os
import time
import datetime
import argparse
import re
import shlex
import json
import ast
import numpy as np
import shutil
import random

############ This is the original directory which is the 'root' directory of your results ##########################
original_path = os.getcwd()
if os.path.isfile("final_result.txt"):
    os.remove("final_result.txt")
    
parser= argparse.ArgumentParser()

parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection", type = str)
parser.add_argument("-seq","--sequence-file",help = "input the sequence file", type = str)
parser.add_argument("-SFAC","--atom-type", help = "input the name of atom of this SAD (case insensitive)" , type = str)
#removing the queue arg because is passed from sbatch submission script
#parser.add_argument("-q", "--queue", help = "input the computing queue you want to use", type = str)
parser.add_argument("-n", "--number-of-cores", help = "input the number of core you want to use", type = str)
parser.add_argument("-na", "--number-of-atoms", help = "input number of heavy anomalous scatterers", type = int)
parser.add_argument("-DSUL_R","--disulfide-range", nargs = '+', help='input the range of the disulfide range (optional)',type = str)
parser.add_argument("-RESOL_R","--resolution-range", nargs = '+', help='input the range of the resolution range (optional)',type = str)
parser.add_argument("-THRE_R","--threshold-range", nargs = '+', help='input the range of the threshold range (optional)',type = str)
parser.add_argument("-ATOM_R","--atom-range", nargs = '+', help='input the range of the atom number range (optional)',type = str)
parser.add_argument("-AutoBuild","--AutoBuild-polish", help = 'type N if you do not want it and type Y if you like it',type = str)
parser.add_argument("-shf", "--shifter", help = 'set this if you want to use shifter', type = str)
args = parser.parse_args()
#pr = cProfile.Profile()
#pr.enable()
#parse reflection file
if args.reflection_mtz:
    #reflectionFile = args.reflection_mtz
    #reflectionFile = os.path.abspath(args.reflection_mtz)
    job_name = args.reflection_mtz.replace('.mtz','')
    reflectionFile = os.path.abspath(args.reflection_mtz)
    print(reflectionFile)
else:
    print('Please input the mtz file')
    sys.exit()

if args.sequence_file:
    sequenceFile = os.path.abspath(args.sequence_file)
else:
    print('Please input sequence')
    sys.exit()
    
if args.atom_type:
    atomType = (args.atom_type).upper()
else:
    print('Please input the scattering atoms')
    sys.exit()

#if args.queue:
#    computeQueue = args.queue
#else:
#    print ('Please input the computing queue you want to use')
#    sys.exit()
    
if args.number_of_cores:
    coreNumber = args.number_of_cores
else:
    print ('Please address the core number you want to use')
    sys.exit()

if args.shifter:
    print('You have chosen to use shifter to run Se_SAD part')
    command_prefix = 'srun -n 1 shifter /img/load_everything.sh'
else:
    print('You have chosen to run using a conda env')
    command_prefix = 'srun -n 1'
###################### Setup the Grid Range #########################

print("This program utilize two softwares, CCP4i2 SHELXC/D and Crank2. Further refinement can be done by initiating Autobuild")

# Define the minimum distance between atoms and symmetry units. Users can change the definition in the parameter.json file in the 
# same directory.
with open('parameter.json') as json_file:
    data = json.load(json_file)
    for p in data['MIND1']:
        if p == atomType:
            mind_atom = str(data['MIND1'][atomType])
            
    for p in data['MIND2']:
        if p == atomType:
            mind_symm = str(data['MIND2'][atomType])
            
    for p in data["Low Resolution CutOff"]:
        if p == atomType:
            low_resolution_cut = str(data['Low Resolution CutOff'][atomType])
            
#Define the SE, S, and DSUL number from the sequence file
sequence=[]
with open(sequenceFile,'r') as f:
    for line in f:
        sequence.append(line.rstrip('\n'))
        
for line in sequence:
    if '>' in line:
        sequence.remove(line)

protein = ''.join(sequence)

single_S_or_SE_number=0
double_sulfur_number=0

for i in protein:
    if i=='M':
        single_S_or_SE_number+=1 #For S-SAD, M = S; For Se_SAD, M=SE
    elif i=='C':
        double_sulfur_number+=1
        
max_DSUL = 0
SE_num = 0
if single_S_or_SE_number == 0:
    if args.number_of_atoms:
        SE_num = args.number_of_atoms
        max_DSUL = double_sulfur_number/2
        max_S = single_S_or_SE_number+double_sulfur_number
    else:
        print("There is no methionine in this protein. Please enter the number of heavy atoms using -ATOM_R")
        #sys.exit()

else :
    max_DSUL = double_sulfur_number//2 #changing for python3
    max_S = single_S_or_SE_number+double_sulfur_number
    max_SE = single_S_or_SE_number

#######################################################
########calling shifter for phenix this step is executed only once at the beginning of the workflow

process = subprocess.Popen(command_prefix + ' phenix.mtz.dump '+reflectionFile, 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()
print(err)
print(out)
split_out=out.splitlines()

for line in split_out:
    if b'Resolution range' in line:
        print(line)
        resolution = round(float(line.split(b' ')[-1]),1)    #changing for python3
        

DSUL_range = range(1,max_DSUL+1)
resolution_range = np.arange(resolution, resolution+1.0, 0.1) 

if atomType == 'S':
    atom_find = range(max_S//2, max_S+1) #changing for python3
else:
    if single_S_or_SE_number == 0:
        atom_find = range(SE_num-1, SE_num+2)
    else:
        atom_find = range(max_SE/2, max_SE+1)

    
    
thre_range = np.linspace(0.2,0.4,3) #(0.2,0.5,4)


################################# Update the grid range if asked by users ##############################

def get_range(x,y):
    if x == y:
        return np.array([x])
    else:
        return np.arange(x,y+0.1,0.1)

if args.disulfide_range:
    if atomType == 'S':
        DSUL_range = range(int(args.disulfide_range[0]),int(args.disulfide_range[1])+1)
    else:
        print('The atom type is not S. Input disulfide range will not be used')
else:
    print("Defualt DSUL search range will be used")
    
    
if args.resolution_range:
    resolution_range = get_range(round(float(args.resolution_range[0]),1),round(float(args.resolution_range[1]),1))
else:
    print("Defualt resolution search range will be used")
    
    
if args.atom_range:
    atom_find = range(int(args.atom_range[0]),int(args.atom_range[1])+1)
else:
    print("Defualt atom number search range will be used")
    
if args.threshold_range:
    thre_range = get_range(round(float(args.threshold_range[0]),1),round(float(args.threshold_range[1]),1))
else:
    print("Defualt threshold search range will be used")

################################# Randomize the job submission ########################################
print ("Disulfide search range: "+str(DSUL_range))
print ("Threshold cut-off range: "+str(thre_range))
print ("Resolution search range: "+str(resolution_range))
print ("Atom number search range:"+str(atom_find))

#sys.exit() #to be deleted



directory_list = []
command_list = []
se_SAD = os.path.abspath(os.getcwd()+'/Se_SAD_automation.py')
crank = os.path.abspath(os.getcwd()+'/crank2_script.py')
shelx = os.path.abspath(os.getcwd()+'/SHELX_script.py')
#creat a list of directory and corresponding command
if max_DSUL > 0 and atomType == 'S':
    for dsul in DSUL_range:
        for thre in thre_range:
            for resolution in resolution_range:
                for number in atom_find:
                    directory = 'DSUL'+str(dsul)+'/threshold'+str(thre)+'/resolution'+str(resolution)+'/atom_number'+str(number)
                    directory_list.append(directory)
                    automation_cl = 'python '+se_SAD+' -shelx '+shelx+ ' -crank '+crank+' -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.5 -thre '+str(thre)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -P '+original_path
                    #automation_cl = 'python ' +se_SAD+ '-crank '+crank+' -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.3 -thre '+str(thre)+' -DSUL '+str(dsul)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -lresl '+low_resolution_cut+' -P '+original_path
                    #automation_cl = 'strace -tt -f -etrace=file,read,write,close,lseek,ioctl -o out.log python Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.3 -thre '+str(thre)+' -DSUL '+str(dsul)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -lresl '+low_resolution_cut+' -P '+original_path
                    command_list.append(automation_cl)


#Do not consider DSUl parameter
elif max_DSUL == 0 or atomType != 'S':
    print ('0 disulfide bond found or the atom type you are looking for is not Sulfur. Do not consider the grid search for disulfied number')
    for thre in thre_range:
        for resolution in resolution_range:
            for number in atom_find:
                directory = 'threshold'+str(thre)+'/resolution'+str(resolution)+'/atom_number'+str(number)
                directory_list.append(directory)
                #automation_cl = 'python Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.5 -thre '+str(thre)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -P '+original_path
                automation_cl = 'python '+se_SAD+' -shelx '+shelx+ ' -crank '+crank+' -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.5 -thre '+str(thre)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -P '+original_path
                #print(automation_cl)
                #print(sequenceFile)
                #automation_cl = 'strace -tt -f -etrace=file,read,write,close,lseek,ioctl -o out.log  python Se_SAD_automation.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -resl '+str(resolution)+' -FIND '+str(number)+' -ESEL 1.5 -thre '+str(thre)+' -SFAC '+atomType+' -MIND1 '+mind_atom+' -MIND2 '+mind_symm+' -P '+original_path
                command_list.append(automation_cl)

# Shuffle the jobs
matching = list(zip(directory_list,command_list))
          
random.shuffle(matching) 

directory_list,command_list = zip(*matching)

# run jobs
start_time_cp = datetime.datetime.now()
for i in range(5):     #set to 5 for testing/debug purposes
#for i in range(len(directory_list)):
    print(len(directory_list))
    start_time = datetime.datetime.now()
    os.system('mkdir -p -m777 '+directory_list[i])
    #os.system('mkdir -p '+directory_list[i])
    # removing redundant copy operations
    #os.system('cp Se_SAD_automation.py SHELX_script.py crank2_script.py autobuild.py '+sequenceFile+' '+reflectionFile+' '+directory_list[i])  
    #os.system('cp Se_SAD_automation.py SHELX_script.py crank2_script.py autobuild.py '+directory_list[i])  
    end_time = datetime.datetime.now()
    delta = end_time - start_time
    print('Time spent for copying before each job')
    print(str(float(delta.seconds) + float(delta.microseconds) / 1000000))
    os.chdir("./"+directory_list[i])
    print(os.getcwd())
    ##this will block replacing with Popen
    #os.system('srun -A m3506 -C haswell -q '+computeQueue+' --cores-per-socket '+coreNumber+' -o %J.log '+command_list[i])
    print(command_list[i])
    f = open('output.log', 'w') 
    #process = subprocess.Popen(command_list[i], stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell= True)
    #process = subprocess.Popen('srun --cores-per-socket '+coreNumber+' -o %J.log '+command_list[i], stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    ########run shifter with memory per job plus craynetwork settings so multiple jobs can ran on the same node
    process = subprocess.Popen('srun -n 1 --mem=10000 --gres=craynetwork:0 --cpus-per-task='+coreNumber+' -o %J.log shifter /img/load_everything.sh '+command_list[i] , stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    #introducing some waiting time between job submissions
    time.sleep (5)
    #d = dict(os.environ)
    #process = subprocess.Popen('srun -n 1 -o %J.log shifter /img/load_everything.sh '+command_list[i],env = d, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    print('Job submitted')

    ## only uncomment this for debug otherwise job submission will block
    #out,err = process.communicate()
    #print(out)
    #process
    #print(err)
    os.chdir(original_path)    
end_time_cp = datetime.datetime.now()
delta = end_time_cp - start_time_cp
# this is for measurement only 
with open('times500_2_jobs_optimization.log', 'w') as f:
    f.write('Total time for copy operation:' + str(float(delta.seconds) + float(delta.microseconds) / 1000000)+'\n')
#################### Wait for to start Autobuild to polish the model ###########################
##helper method

##Jobs count method 
def job_count():
    finished_jobs = []
    try:
        with open('final_result.txt','r') as f:
            for line in f:
                finished_jobs.append(line.replace('\n',''))

        return len(finished_jobs)
    except:
        return 0

##Job selected to do autobuild
def case_select ():
    results = []
    with open('final_result.txt','r') as f:
        for line in f:
            results.append(line.replace('\n',''))

    R_score = []
    for result in results:
        R_score.append(round(float(result.split('R_free:')[-1].split('/')[0]),3))

    R_score = np.asarray(R_score)
    polish_cases_Rfree = np.where(R_score == R_score.min())[0]

    polish_cases_0 = []
    for i in polish_cases_Rfree:
        polish_cases_0.append(results[i])
    print(polish_cases_0)
    Residue_score = []
    for i in polish_cases_0:
        Residue_score.append(int(i.split('Residue:')[-1]))

    Residue_score = np.asarray(Residue_score)
    polish_cases_residue = np.where(Residue_score == Residue_score.max())[0]

    polish_cases_1 = []
    for i in polish_cases_residue:
        polish_cases_1.append(polish_cases_0[i])

    #return polish_cases_1
    resolution_filter = []
    if len(polish_cases_1) > 1:
        for i in polish_cases_1:
            resolution_filter.append(float(i.split('resolution')[-1].split('/atom')[0]))
        resolution_filter = np.asarray(resolution_filter)
        polish_cases_resolution = np.where(resolution_filter == resolution_filter.min())[0]
    
    
        polish_cases_2 = []
        for i in polish_cases_resolution:
            polish_cases_2.append(polish_cases_1[i])
            #print("have to pick randomly")
        if len(polish_cases_2) > 1:
            print("have to pick randomly")
            select_case = polish_cases_2[random.randint(0,len(polish_cases_2)-1)]

        else: select_case = polish_cases_2[0]
            
    else:
        select_case = polish_cases_1[0]

    print ("The best grid right now is "+select_case.split('R:')[0])
    print ("It has the score R:"+select_case.split('R:')[-1])
    return select_case.split('R:')[0]

###set to 5 for testing/debug
Total_jobs = 5
#Total_jobs = len(directory_list)
Half_total_jobs = Total_jobs // 2
percent97_total_jobs = Total_jobs * 97 //100
half_finished = False
percent97_finished = False

while half_finished == False:
    
    if job_count() < max(1,Half_total_jobs):
        pass
    else:
	
        selected_job_directory1 = case_select()
        if args.AutoBuild_polish != 'N':
            print('Hello2')
            os.system("mkdir Autobuild1")
            os.system("cp autobuild.py "+sequenceFile+" "+reflectionFile+" "+selected_job_directory1+"result.pdb Autobuild1")
            os.chdir("Autobuild1")
            autobuild_cl = 'python autobuild.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -rfff 0.05 -nproc '+str(coreNumber)+' -pdb result.pdb'
            ### this will block changing it with Popen
            #os.system('srun -A m3506  -C haswell -q '+computeQueue+' --cores-per-socket= '+coreNumber+' -o %J.log '+autobuild_cl)
            os.chdir(original_path)
            half_finished = True


while percent97_finished == False:
    
    if job_count() < max(1,percent97_total_jobs):
        pass
    else:
        percent97_total_jobs = True
        selected_job_directory2 = case_select()
        if selected_job_directory2 == selected_job_directory1:
            break
        else:
            if args.AutoBuild_polish != 'N':
                print('Hello')
                os.system("mkdir Autobuild2")
                os.system("cp autobuild.py "+sequenceFile+" "+reflectionFile+" "+selected_job_directory2+"result.pdb Autobuild2")
                os.chdir("Autobuild2")
                autobuild_cl = 'python autobuild.py -rfl '+reflectionFile+' -seq '+sequenceFile+' -rfff 0.05 -nproc '+coreNumber+' -pdb result.pdb'
                #os.system('srun -A m3506 -C haswell -q '+computeQueue+' --cores-per-socket '+coreNumber+' -o %J.log '+autobuild_cl)
                process = subprocess.Popen('srun -A m3506 -C haswell -q '+computeQueue+' --cores-per-socket '+coreNumber+' -o %J.log '+autobuild_cl, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out,err = process.communicate()
                print(out)
                print(err)
                os.chdir(original_path)
                percent97_finished = True


pr.disable()
with open( 'cpu_out.txt', 'w') as output_file:
    sys.stdout = output_file
    pr.print_stats( sort='time' )
    sys.stdout = sys.__stdout__



   
        
