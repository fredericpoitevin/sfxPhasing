from __future__ import print_function
import sys
import subprocess
import os
import argparse
import re
import json
import ast

if os.path.isfile("output.eff"):
    os.remove("output.eff")

if os.path.isfile("final_result.txt"):
    os.remove("final_result.txt")


current_path = os.getcwd()

with open('FILE_SETUP.json') as json_file:
    parameter = json.load(json_file)
    component_num = int(parameter['component number'])

parser= argparse.ArgumentParser()

#########################################CREATE DICTIONARY ###############################################
my_temp={"hklin":"None", 
      
         "labin":"FP,SIGFP",
      
         "resolution_cutoff": "2.5",
      
         "no_tncs":"False",
      
      "crystal_symmetry":{
          "unit_cell":"None",
          "space_group":"None",
          
      },
      
      "mode" : "*quick full",
      
      "symmetry_exploration":"pointgroup enantiomorph *dataset",
      
      "output":{
          "root":"mrage",
          "max_solutions_to_write" : "1",
          "gui_base_dir": "None",
          "save_aniso_data": "False",
          "job_title":"None"
        },
      
      "composition":{
          "count":"None",
          "component": {
                "sequence" : "None",
                "mtype" : "*protein dna rna",
                "stoichiometry" : "1",
              
                "ensemble": {
                      "coordinates": {
                            "pdb" : "None",
                            "identity" :"None",
                            "error_translation" : "*default floored_chothia_and_lesk rmsd"
                        },
                    
                      "alignment":"None",
                },
              
                "model_collection":{
                      "trim" : "False",
                      "coordinates" :{
                            "pdb" : "None",
                            "identity": "None",
                            "error_translation" : "*default floored_chothia_and_lesk rmsd"
                      }
                },
              
                "template": {
                      "pdb":"None",
                      "alignment" : "None"
                },
              
                "template_list":"None",
              
                "homology":{
                      "file_name":"None",
                      "max_hits": "3"
                },
              
                "search": {
                      "services":"local ncbi",
                      "max_hits": "3"
                },
            }
        },
      
        "assembly": {
              "use_assembled" : "False",
              "local_search": {
                    "sweep" : "5",
                    "extent" : "5"
              },
              "component": {
                    "file_name" : "None",
                    "operation": {
                          "euler" : "0 0 0",
                          "displacement" : "0 0 0"
                    }
              }
        },
        "queue": {
            "technology" : "lsf *threading sge multiprocessing pbs",
            "cpus":"12",
            "submission_command": "None",
            "qslot_cpus":"1"
        },  
     
        "packing_pool": "20",
        "rotation_peaks_cutoff" : "0.75",
        "post_refinement_cutoff" : "0.75",
        "final_selection_cutoff" : "0.75",
        "b_factor_refinement" : "True",
        "significant_peak_threshold": "7.0",
        "sculptor_protocols" : "11 10 13 *12 1 3 2 5 4 7 6 9 8 all minimal",
        "template_equivalence": "True",
        "assembly_acceptance_policy" : "observed *always never",
        "assembly_ensemble_creation_policy":"*observed always never",
        "pdb_mirror": "*rcsb pdbe pdbj",
        "exclude_pdb_ids": "None",
        "simple_run": {
            "enable": "False",
            "component_sequence" : "None",
            "template_model": "None",
            "alignment":"None",
            "blast_services":"local ncbi",
            "max_blast_hits": "3"
        }
 
     }
########################################################################################################################################################




######################################################## ADD ALL ARGUMENT BELOW ######################################################################
#parse the mtz file
parser.add_argument("-rfl","--reflection-mtz", help="input the mtz file of reflection",type = str)

#parse the data label
parser.add_argument("-labin","--data-labels", help="enter the data lables shown in the reflection file", type = str)

#parse number of components
parser.add_argument("-n","--component-count", default = 1, help="input the number of chains",type = int)

#parse count of chains
parser.add_argument("-c","--chains-count", default = 1, help="input the number of chains",type = int)


#create resolution cutoff parse
parser.add_argument("-res","--resolution", default = 2.5, help="input the resolution cut off", type = float)

# input the root directory
parser.add_argument("-P", "--path", help = "input the orginal path", type = str)

parser.add_argument("-cpus", "--number-of-cores", help = "input the number of core you want to use", type = str)

for i in range(1,component_num+1):#10:
    #create pdb argument parse
    j=str(i)
    parser.add_argument("-pdbE"+j,"--component"+j+"-ensemble-pdb",help="input the pdb"+j+" file of ensemble",type = str)
    parser.add_argument("-idenE"+j,"--component"+j+"-ensemble-identity", help="input the identity"+j+" file of ensemble",type = float)
    parser.add_argument("-errtE"+j,"--component"+j+"-ensemble-error-translation",help="input the error translation"+j+" file of ensemble",type = str)

    parser.add_argument("-pdbM"+j,"--component"+j+"-model-pdb", help="input the pdb"+j+" file of model",type = str)
    parser.add_argument("-idenM"+j,"--component"+j+"-model-identity", help="input the identity"+j+" file of model",type = float)
    parser.add_argument("-errtM"+j,"--component"+j+"-model-error-translation",help="input the error translation"+j+" file of model",type = str)

    parser.add_argument("-pdbT"+j,"--component"+j+"-template-pdb",help="input the pdb"+j+" file of template",type = str)
    parser.add_argument("-idenT"+j,"--component"+j+"-template-identity",help="input the identity"+j+" file of template",type = float)
    parser.add_argument("-errtT"+j,"--component"+j+"-template-error-translation",help="input the error translation"+j+" file of template",type = str)

    #create seq argument parse
    parser.add_argument("-seq"+j,"--component"+j+"-sequence",help="input the sequence"+j+" file of ensemble",type = str) 

    #create homology argument parse
    parser.add_argument("-hom"+j,"--component"+j+"-homology",help="input the pdb"+j+" file of ensemble",type = str)


#########################################################################################################



######################################### PASS ARGUMENT FROM COMMAND LINE ##########################################################
args = parser.parse_args()


if args.reflection_mtz:
    my_temp['hklin'] = args.reflection_mtz
if args.data_labels:
    my_temp['labin'] = args.data_labels

number_of_chains = 1
component_number = 1

if args.component_count:
    number_of_chains = component_num
    for i in range(number_of_chains):
        component_number = str(i+1)
        dictionary_string = str(my_temp['composition']['component'])
        my_temp['composition']['component'+component_number] = ast.literal_eval(dictionary_string)
    my_temp['composition'].pop('component', None)

if args.resolution:
    my_temp['resolution_cutoff'] = args.resolution

if args.chains_count:
    my_temp['composition']['count'] = args.chains_count

if args.number_of_cores:
    my_temp['queue']['cpus'] = args.number_of_cores

for i in range(1, component_num +1):

    #get attribute of pdb
    get_arg_ensemble_pdb = getattr(args, 'component'+str(i)+'_ensemble_pdb')
    get_arg_ensemble_id = getattr(args, 'component'+str(i)+'_ensemble_identity')
    get_arg_ensemble_et = getattr(args, 'component'+str(i)+'_ensemble_error_translation')
    get_arg_model_pdb = getattr(args, 'component'+str(i)+'_model_pdb')
    get_arg_model_id = getattr(args, 'component'+str(i)+'_model_identity')
    get_arg_model_et = getattr(args, 'component'+str(i)+'_model_error_translation')
    get_arg_template_pdb = getattr(args, 'component'+str(i)+'_template_pdb')
    get_arg_template_id = getattr(args, 'component'+str(i)+'_template_identity')
    get_arg_template_et = getattr(args, 'component'+str(i)+'_template_error_translation')

    #get attribute of sequence
    get_arg_sequence = getattr(args,'component'+str(i)+'_sequence')

    #get attribute of homology
    get_arg_homology = getattr(args,'component'+str(i)+'_homology')

    #pass ensemble parameter
    if get_arg_ensemble_pdb:
        my_temp['composition']['component'+str(i)]['ensemble']['coordinates']['pdb'] = get_arg_ensemble_pdb

    if get_arg_ensemble_id:
        my_temp['composition']['component'+str(i)]['ensemble']['coordinates']['identity'] = get_arg_ensemble_id

    if get_arg_ensemble_et:
        my_temp['composition']['component'+str(i)]['ensemble']['coordinates']['error_translation'] = get_arg_ensemble_et

    #pass model parameter
    if get_arg_model_pdb:
        my_temp['composition']['component'+str(i)]['model_collection']['coordinates']['pdb'] = get_arg_model_pdb

    if get_arg_model_id:
        my_temp['composition']['component'+str(i)]['model_collection']['coordinates']['identity'] = get_arg_model_id

    if get_arg_model_et:
        my_temp['composition']['component'+str(i)]['model_collection']['coordinates']['error_translation'] = get_arg_model_et

    #pass template
    if get_arg_template_pdb:
        my_temp['composition']['component'+str(i)]['template']['pdb'] = get_arg_template_pdb

    if get_arg_template_id:
        my_temp['composition']['component'+str(i)]['template']['identity'] = get_arg_template_id

    if get_arg_template_et:
        my_temp['composition']['component'+str(i)]['template']['error_translation'] = get_arg_template_et

    #pass sequence
    if get_arg_sequence:
        my_temp['composition']['component'+str(i)]['sequence'] = get_arg_sequence

    #pass homology
    if get_arg_homology:
        my_temp['composition']['component'+str(i)]['homology']['file_name']=get_arg_homology

################################################################################################################################




########################################### GET INFORMATION FROM MTZ FILE ######################################

process = subprocess.Popen('phenix.mtz.dump '+my_temp['hklin'], 
                          stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,shell=True)

out,err = process.communicate()
print(out)
print(err)
split_out=out.splitlines()

mylist = []
for i in range(len(split_out)):
    item = split_out[i].decode("utf-8")
    if item != '':
        mylist.append(split_out[i].decode("utf-8"))

my_mtz={}
for each_line in mylist:
    if len(each_line.split(":"))== 2:
        key, value = each_line.split(":")[0],each_line.split(":")[1]
        my_mtz[key] = value

#find the unit cell and symmetry parameter
unit_cell_search = 'unit cell'

for key in my_mtz:
      if unit_cell_search in key.lower():
            my_mtz[key]=my_mtz[key].replace('(','').replace(')','').replace(',',' ')
            my_temp['crystal_symmetry']['unit_cell']= my_mtz[key]

space_group_search = 'space group symbol'

for key in my_mtz:
      if space_group_search in key.lower():
            my_temp['crystal_symmetry']['space_group']= my_mtz[key]


####################################### CREATING OUTPUT.EFF####################################################
for i in my_temp:
    my_string = i + '=' +json.dumps(my_temp[i])

    for j in range (1,component_num+1):
        k = str(j)
        if 'component'+k in my_string:
            my_string = my_string.replace('component'+k,'component')

    my_string=my_string.replace('}','\n}')
    my_string=my_string.replace(':','=')
    my_string=my_string.replace('={','{\n')
    my_string=my_string.replace('= {','{\n')
    my_string=my_string.replace('"','')
    my_string=my_string.replace('mrage','"mrage"')

    if i != 'labin':
        my_string=my_string.replace (',', '\n')
    print(my_string,file=open("output.eff", "a"))
    print('',file=open("output.eff", "a"))

##############################################################################################################




####################################### DEPLOY PHENIX COMMAND LINE#############################################
process = subprocess.Popen(' phaser.MRage --verbosity=VERBOSE output.eff', 
                         stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,shell=True)

#os.system('phaser.MRage --verbosity=VERBOSE output.eff')
#redirecting phaser output to file 
out,err = process.communicate()
print(out)
print(err)
split_out=out.splitlines()
with open('test.txt', 'w') as f:
    for i in range(len(split_out)):
     f.write(split_out[i].decode("utf-8")+'\n')
####################################################################################################################

#######################Output result to the root directory#########################################
for item in os.listdir(os.getcwd()):
    #if '.log' in item:
    if '.txt' in item:
        logFile = item

mylog = []

with open(logFile) as f:
    for lines in f:
        mylog.append(lines)
print('CHECK'+mylog[0])
print(len(mylog))
try:
    print('Here I am')
    resultTitle_index = mylog.index('Evaluation for probability of solution being correct:\n')
    print(resultTitle_index)
except:
    print ('The case is not successful')

for line in mylog:
    if 'P(total)=' in line:
        print(line)
        end_of_result = mylog.index(line)

result = []
P_total = mylog[end_of_result].replace('\n',' ')

result.append(P_total)

for i in mylog[resultTitle_index+1:end_of_result]:
    try:
        prob = i.split("=>")[-1].split('=')[-1].replace('\n',' ')
        TFZ = i.split("=>")[0].split('=')[-1].replace('\n','')
        result.append('TFZ = '+TFZ+' Probability = '+prob)
    
    except:
        pass

result_statement = current_path+' '
for i in result:
    result_statement += ''+i+''
    
os.chdir(args.path)
print(result_statement.replace(args.path,''),file=open('final_result.txt','a'))
