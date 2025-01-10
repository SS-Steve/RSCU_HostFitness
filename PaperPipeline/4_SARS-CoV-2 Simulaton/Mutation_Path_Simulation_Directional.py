### System requirement ###

import numpy as np
import pandas as pd
from scipy import stats
import sklearn
import math
import os
import multiprocessing
import tqdm
import pickle

import src.ss_utils as ss
import src.ss_database as db
import src.mutation_utils as mu

# Linux commands
# nohup python Mutation_Path_Simulation_Directional.py > Mutation_Path_Simulation_Directional.log 2>&1 &

Virus_Host_Dict = {'Bat Hp-betacoronavirus/Zhejiang2013'                 : 'not_human',
                   'Bat coronavirus BM48-31/BGR/2008'                    : 'not_human',
                   'Betacoronavirus England 1'                           : 'human',
                   'Middle East respiratory syndrome-related coronavirus': 'human',
                   'Betacoronavirus Erinaceus/VMC/DEU/2012'              : 'not_human',
                   'Betacoronavirus HKU24'                               : 'not_human',
                   'Bovine coronavirus'                                  : 'human',
                   'Human coronavirus OC43'                              : 'human',
                   'Human coronavirus HKU1'                              : 'human',
                   'Murine hepatitis virus'                              : 'not_human',
                   'Rat coronavirus Parker'                              : 'not_human',
                   'Pipistrellus bat coronavirus HKU5'                   : 'not_human',
                   'Rabbit coronavirus HKU14'                            : 'not_human',
                   'Rousettus bat coronavirus'                           : 'not_human',
                   'Rousettus bat coronavirus HKU9'                      : 'not_human',
                   'SARS coronavirus Tor2'                               : 'human',
                   'Severe acute respiratory syndrome coronavirus 2'     : 'human',
                   'Tylonycteris bat coronavirus HKU4'                   : 'not_human'}

###============================================================================================================###
### Input parameters setting ###
###============================================================================================================###

Target_Virus_Name = 'Rousettus bat coronavirus HKU9'

Covid_Virus_Name = 'Severe acute respiratory syndrome coronavirus 2'

MUTATION_PATH_DIRECTION = 'HIGH2LOW'   # 'LOW2HIGH' or 'HIGH2LOW'

Dict_HIP_Trace = {}
TERMINATION_ACT = False
Termination_Iteration_n = 200 # The simulation is aborted if no improvement in those numbers of iterations
Reverse_Filtering_Decision = True

CPU_USED = os.cpu_count() - 5
# CPU_USED = int(os.cpu_count() / 2)
# CPU_USED = 5

Model_Path = f"Trained_Model//RF_model_RSCUpTAXONOMYpCDS_LENGTH_human.pickle"

SAVE_FOLDER = 'Mutation_Path_Simulation' # In 'Output' folder 


###============================================================================================================###
### Program start ###
###============================================================================================================###


# Determine 'LOW2HIGH' or 'HIGH2LOW' mode
if isinstance(MUTATION_PATH_DIRECTION, str) and MUTATION_PATH_DIRECTION.upper() in ['LOW2HIGH']:
    HIP_Indicator = 1
elif isinstance(MUTATION_PATH_DIRECTION, str) and MUTATION_PATH_DIRECTION.upper() in ['HIGH2LOW']:
    HIP_Indicator = 0
else:
    raise ValueError(f"'MUTATION_PATH_DIRECTION' cannot only be 'LOW2HIGH' or 'HIGH2LOW'.")

# Import data

RSCU = ss.load_RTG_Data("RSCU")
CDS_LENGTH = ss.load_RTG_Data("CDS_LENGTH")
INFO = ss.load_RTG_Data("INFO")
TAXONOMY = ss.load_RTG_Data("TAXONOMY")

# Load pre-defined list and dictionary
multibox_codon_table_list = db.multibox_codon_table_list[:]
no_stop_mut_table_list = db.no_stop_mut_table_list.copy()
empty_Codon_Change_Dict = mu.empty_Codon_Change_Dict.copy()
Virus_Accession_Dict = mu.Virus_Accession_Dict.copy()

# Load model 
with open(Model_Path, 'rb') as file:
    Model = pickle.load(file)

### Target virus ###

# Target_Virus_Name = 'Bat coronavirus BM48-31/BGR/2008'
Target_NCBI_ID = Virus_Accession_Dict[Target_Virus_Name]

# Calculate Target Virus HIP
Target_Virus_Index = (INFO['id'] == Target_Virus_Name)
Target_RSCU = RSCU[Target_Virus_Index]
X_Target_Virus = pd.concat([Target_RSCU, 
                            TAXONOMY[Target_Virus_Index], 
                            CDS_LENGTH[Target_Virus_Index]], axis=1)
X_Target_Virus = np.array(X_Target_Virus)

Target_HIP = Model.predict_proba(X_Target_Virus)[0][HIP_Indicator]

TARGET_Info_Dict = {'Virus_Name':Target_Virus_Name,
                    'NCBI_ID':Target_NCBI_ID,
                    'RSCU':Target_RSCU,
                    'HIP':Target_HIP}


### SARS-CoV-2 ###

# Covid_Virus_Name = 'Severe acute respiratory syndrome coronavirus 2'
Covid_NCBI_ID = Virus_Accession_Dict[Covid_Virus_Name]

# Calculate SARS-CoV-2 HIP
Covid_Virus_Index = (INFO['id'] == Covid_Virus_Name)
Covid_RSCU = RSCU[Covid_Virus_Index]
Betacoronavirus_TAXONOMY = TAXONOMY[Covid_Virus_Index]

X_Covid = pd.concat([Covid_RSCU, 
                     Betacoronavirus_TAXONOMY, 
                     CDS_LENGTH[Covid_Virus_Index]], axis=1)
X_Covid = np.array(X_Covid)

Covid_HIP = Model.predict_proba(X_Covid)[0][HIP_Indicator]

COVID_Info_Dict = {'Virus_Name':Covid_Virus_Name,
                   'NCBI_ID'  :Covid_NCBI_ID,
                   'RSCU'     :Covid_RSCU,
                   'HIP'      :Covid_HIP}

# Reset index for future use
Betacoronavirus_TAXONOMY = Betacoronavirus_TAXONOMY.reset_index(drop=True)


#############################################################################
# Setup Start and Stop Virus (Determine Forward or Reverse Mutation Path) ###
#############################################################################

if HIP_Indicator == 1:
    Start_Virus_Info_Dict, Stop_Virus_Info_Dict = (TARGET_Info_Dict, COVID_Info_Dict)
elif HIP_Indicator == 0:
    Start_Virus_Info_Dict, Stop_Virus_Info_Dict = (COVID_Info_Dict, TARGET_Info_Dict)


Start_Virus_Name = Start_Virus_Info_Dict['Virus_Name']
Start_NCBI_ID = Start_Virus_Info_Dict['NCBI_ID']
Start_RSCU = Start_Virus_Info_Dict['RSCU']
Start_HIP = Start_Virus_Info_Dict['HIP']

Stop_Virus_Name = Stop_Virus_Info_Dict['Virus_Name']
Stop_RSCU = Stop_Virus_Info_Dict['RSCU']
Stop_HIP = Stop_Virus_Info_Dict['HIP']

Target_RSCU_dict = mu.X_RSCU_to_RSCU_dict(Stop_RSCU)

# Download newest CDS from NCBI for Target virus
if isinstance(Start_NCBI_ID, list) is True: # for Segmented viruses
    start_whole_seq = ''
    start_all_cds_dict = {}
    for NCBI_ID in Start_NCBI_ID:
        all_cds_dict_t = ss.get_all_cds(NCBI_ID)
        for pid, cds_seq_t in all_cds_dict_t.items():
            if pid not in start_all_cds_dict:
                start_all_cds_dict[pid] = cds_seq_t
            else:
                while pid in start_all_cds_dict:
                    pid += '_'
                start_all_cds_dict[pid] = cds_seq_t
        whole_seq_t = ss.get_whole_seq(NCBI_ID)[NCBI_ID]
        start_whole_seq += whole_seq_t
else:
    start_all_cds_dict = ss.get_all_cds(Start_NCBI_ID)
    start_whole_seq = ss.get_whole_seq(Start_NCBI_ID)[Start_NCBI_ID]


### Get Contecated CDS and removed the start and stop codon ###
# Because the start and stop codon cannot be mutated in the simulation
n_cds = len(start_all_cds_dict)
f = ''
START_COVID_CDS_LEN, START_CND_SINGLE = {}, {}
for _id_, cds_seq in start_all_cds_dict.items():
    START_COVID_CDS_LEN[_id_] = len(cds_seq)
    cds = cds_seq[3:len(cds_seq)-3]
    START_CND_SINGLE[_id_] = ss.codon_n_dict(cds)
    f += cds
f = f + (n_cds * ('___' + '___'))


# Generate Codon_N_Dict and RSCU_Dict
START_CND  = ss.codon_n_dict(f)
START_RSCU = ss.codon_n_dict_to_RSCU(START_CND)
START_RSCU = ss.dict_extract_list(START_RSCU, multibox_codon_table_list)

# Generate CDS-Length data
START_COVID_CDS_LEN_Dict = mu.single_cds_len_dict_to_COVID_LEN_X(START_COVID_CDS_LEN, start_whole_seq, start_all_cds_dict)

# Generate the Start X array 
START_X = pd.concat([pd.DataFrame.from_dict({0:START_RSCU}, 'index'), 
                     Betacoronavirus_TAXONOMY,
                     pd.DataFrame.from_dict({0:START_COVID_CDS_LEN_Dict}, 'index')], 
                    axis=1)
START_X = np.array(START_X)


Current_Codon_N_Dict        = START_CND.copy()
Current_Codon_N_Dict_Single = START_CND_SINGLE.copy()
#Current_RSCU_dict           = START_RSCU.copy()
Current_LEN_dict            = START_COVID_CDS_LEN_Dict.copy()
Current_CDS_LEN_dict        = START_COVID_CDS_LEN.copy()
#Current_TAX_dict            = COVID_TAX_dict.copy()
Current_Human_PP            = Start_HIP
Second_Lowest_Human_PP      = 2

Stage_Lowest_Data = [{
    'mutation'            : 'None',
    'Codon_N_Dict'        : Current_Codon_N_Dict,
    'COVID_LEN_dict'      : Current_LEN_dict,
    'Codon_N_Dict_Single' : Current_Codon_N_Dict_Single, 
    'CDS_LEN_dict'        : Current_CDS_LEN_dict,
    'Codon_Change_Dict'   : empty_Codon_Change_Dict,
    'pre_mutation'        : [],
}]
 
# Create Mutation Pool containing all the possible mutations
Mutation_Pool = []
for k, mutation_list in no_stop_mut_table_list.items():
    if k in ['-', '+']:
        for mut in mutation_list:
            for cds_loc in list(Current_CDS_LEN_dict.keys()): # Consider mutating different CDS
                mutploc = f"{mut}@{cds_loc}"
                if mutploc not in Mutation_Pool:
                    Mutation_Pool.append(mutploc)
    else:
        for mut in mutation_list:
            if mut not in Mutation_Pool:
                Mutation_Pool.append(mut)
                

# SAVE_FOLDER = 'Mutation_Path_Simulation'
ss.create_folder(SAVE_FOLDER)

Virus_Folder = f"{Start_Virus_Name.replace('/', '_')}-to-{Stop_Virus_Name.replace('/', '_')}"

print (f"### Start_Virus: {Start_Virus_Name} - Start_Virus_HIP = {Start_HIP}")
print (f"### Target_Virus: {Stop_Virus_Name} - Target_Virus_HIP = {Stop_HIP}")
print (f"### Output Folder: {SAVE_FOLDER}//{Virus_Folder}")
#============================================#


ss.create_folder(f"{SAVE_FOLDER}//{Virus_Folder}")
ss.create_folder(f"{SAVE_FOLDER}//{Virus_Folder}//Raw_data")
Before_Clean_Stage_Lowest_Data_n = len(Stage_Lowest_Data)


# Dict_HIP_Trace = {}
# TERMINATION_ACT = False
# Termination_Iteration_n = 500 # The simulation is aborted if no improvement in those numbers of iterations
# Reverse_Filtering_Decision = True

Pre_Human_PP = Start_HIP
Target_PP_Threshold = Stop_HIP
Current_Max_Human_PP, Recorded_Max_Correlation, Max_Correlation, Pre_Max_Correlation = 0, 0, 0, 0
# Target_PP_Threshold = 0.5 # Set stop HIP if you want 

Current_Mutation = 1
Pre_set_Stage_Data_Count = 1

while Current_Human_PP < Target_PP_Threshold and TERMINATION_ACT is False:
# while Current_Human_PP < Target_PP_Threshold and TERMINATION_ACT is False and Current_Mutation < 5: # Limitations on cycles


    # Create iteration folder
    Current_Mutation_Folder = f"Current_Mutation_[{Current_Mutation}]"
    ss.create_folder(f"{SAVE_FOLDER}//{Virus_Folder}//Raw_data//{Current_Mutation_Folder}")

    print (f"# Current_Mutation: {Current_Mutation} - Total Candidate Sequence: {len(Stage_Lowest_Data)} - Before_Target_Corr_Select: {Before_Clean_Stage_Lowest_Data_n}")


    # === Main Program === #

    All_Lowest_Data = []
    Found_Lower_PP = False
    Current_Lowest_Data = None
    res, c = {}, 0
    for loc_in_Stage, Current_Single_Lowest_Data in enumerate(Stage_Lowest_Data):

        print (f"Current_Human_Pred_Proba: {Current_Human_PP}")
        print (f"Current Candidate Sequence: {loc_in_Stage+1} / {len(Stage_Lowest_Data)} - Previous mutation: {Current_Single_Lowest_Data['mutation']}")

        Current_Codon_N_Dict          = Current_Single_Lowest_Data['Codon_N_Dict']
        Current_Codon_N_Each_CDS_Dict = Current_Single_Lowest_Data['Codon_N_Dict_Single']
        Current_LEN_dict              = Current_Single_Lowest_Data['COVID_LEN_dict']
        Current_CDS_LEN_dict          = Current_Single_Lowest_Data['CDS_LEN_dict']
        Current_Codon_Change_Dict     = Current_Single_Lowest_Data['Codon_Change_Dict']
        Previous_Mutation             = Current_Single_Lowest_Data['pre_mutation']

        def Process(ind):

            # Get mutation from the Mutation Pool
            mut = Mutation_Pool[ind]
            Update_Previous_Mutation = Previous_Mutation + [mut]

            Reverse_Filtering = mu.Reverse_Mut_Filtering(Current_Codon_Change_Dict, mut, Reverse_Filtering_Decision)
            codon_change_dict = mu.mutation_dict_add_mutation(Current_Codon_Change_Dict, mut)

            # Update RSCU and LEN data with Mutation
            cds_loc = None
            if '@' in mut:
                cnd_mut, cds_loc = mut.split('@')
                cds_nd = mu.mutate_codon_n_dict_single(Current_Codon_N_Each_CDS_Dict, cnd_mut, cds_loc) # Check if the codon runs out
            else:
                cnd_mut = mut
                cds_nd  = Current_Codon_N_Each_CDS_Dict.copy()
            cnd = ss.codon_n_dict_mutation(Current_Codon_N_Dict, cnd_mut, False)

            if cds_nd is False or cnd is False or codon_change_dict is False or Reverse_Filtering is True:
                output = False
            else:
                rscu = ss.codon_n_dict_to_RSCU(cnd)
                x_rscu = ss.dict_extract_list(rscu, multibox_codon_table_list)
                cld, cds_ld = mu.mutate_len_dict(Current_LEN_dict, mut, Current_CDS_LEN_dict)
            
                # Generate new X and calculate predicted probability
                new_X = pd.concat([pd.DataFrame.from_dict({0:x_rscu}, 'index'), 
                                   Betacoronavirus_TAXONOMY,
                                   pd.DataFrame.from_dict({0:cld}, 'index')], 
                                  axis=1)
                new_X = np.array(new_X)

                pred_proba = Model.predict_proba(new_X)[0][HIP_Indicator]

                # Organise the output dictionary
                data_dict = {'mutation'            : mut,
                             'pre_mutation'        : Update_Previous_Mutation,
                             'Codon_N_Dict'        : cnd,
                             'Codon_N_Dict_Single' : cds_nd,
                             'COVID_LEN_dict'      : cld,
                             'CDS_LEN_dict'        : cds_ld,
                             'Codon_Change_Dict'   : codon_change_dict,
                             'pred_proba'          : pred_proba}
                output = [ind, pred_proba, data_dict]

            return output


        # Multi-processing main program
        pool = multiprocessing.Pool(processes = CPU_USED)
        Iter_List = list(range(len(Mutation_Pool)))
        for _OUTPUT_RESULT_ in tqdm.tqdm(pool.imap_unordered(Process, Iter_List), 
                                         total   = len(Iter_List),
                                         disable = True):
            if _OUTPUT_RESULT_ is not False:
                ind, pred_proba, data_dict = _OUTPUT_RESULT_

                # Get List of Lower PP
                All_Lowest_Data.append(data_dict)

                # Organise all the mutation data into dictionary for dataframe 
                res[c] = {'index'        : ind,
                          'all_data_ind' : c,
                          'loc_in_Stage' : loc_in_Stage,
                          'mutation'     : data_dict['mutation'],
                          'pre_mutation' : ss.list_to_str(data_dict['pre_mutation']),
                          'pred_proba'   : pred_proba}
                c += 1
                
            else:
                pass

        pool.close()

        # Organise into dataframe 
        df_res = pd.DataFrame.from_dict(res, 'index')
        df_res = df_res.sort_values('index', ignore_index=True)

    # Sort dataframe and save 
    df_show = df_res.sort_values('pred_proba', ascending=False, ignore_index=True) 
    ss.save_file(df_show, f"{SAVE_FOLDER}//{Virus_Folder}//Raw_data//{Current_Mutation_Folder}//[{Current_Mutation}]_MUTATION_RES")

    # ==================== #


    # Get the Lowest data in all the results

    df_lowest_above_pp  = ss.extract_df_row_top_n(df_show, 'pred_proba', Pre_set_Stage_Data_Count, ASCENDING=False)

    Stage_Lowest_Data = [All_Lowest_Data[i] for i in df_lowest_above_pp['all_data_ind'].tolist()]

    # Record previous HIP
    Pre_Human_PP = Current_Human_PP
    # Update HIP
    Current_Human_PP = df_show.loc[0, 'pred_proba']
    # Record Max HIP
    if Current_Human_PP > Current_Max_Human_PP:
        Current_Max_Human_PP = Current_Human_PP

    Before_Clean_Stage_Lowest_Data_n = len(Stage_Lowest_Data)


    # Get the most correlate to Target_Virus (RSCU)
    if len(Stage_Lowest_Data) > Pre_set_Stage_Data_Count:

        Recoding_Involve = False
        cc_dict = {}
        cc = 0 
        for loc, data_dict in enumerate(Stage_Lowest_Data):
            if   data_dict['mutation'] in no_stop_mut_table_list[1]:
                Recoding_Involve = True
                Recoding = '1nt'
            elif data_dict['mutation'] in no_stop_mut_table_list[2]:
                Recoding_Involve = True
                Recoding = '2nt'
            else:
                Recoding = 'False'

            Codon_N_Dict = data_dict['Codon_N_Dict']
            t_rscu = ss.codon_n_dict_to_RSCU(Codon_N_Dict)
            t_rscu = ss.dict_extract_list(t_rscu, multibox_codon_table_list)
            corr_coef = ss.corr_coef_pearson(
                list(t_rscu.values()), 
                list(Target_RSCU_dict.values())) # Calculate Pearson Correlation
            cc_dict[cc] = {'ind':loc, 'corr_coef':corr_coef, 'recoding':Recoding}
            cc += 1

        df_cc_t = pd.DataFrame.from_dict(cc_dict, 'index')
        ss.save_file(df_cc_t, f"{SAVE_FOLDER}//{Virus_Folder}//Raw_data//{Current_Mutation_Folder}//[{Current_Mutation}]_Target_Corr_Select_RES")

        if Recoding_Involve is True:
            if '1nt' in df_cc_t['recoding'].tolist():
                df_cc_t = df_cc_t[df_cc_t['recoding'] == '1nt']
            else:
                if '2nt' in df_cc_t['recoding'].tolist():
                    df_cc_t = df_cc_t[df_cc_t['recoding'] == '2nt']
                else:
                    pass
            df_cc_t = df_cc_t.reset_index(drop=True)
        else:
            pass

        df_cc_t, _ = ss.extract_df_row_top_n(df_cc_t, 'corr_coef', Pre_set_Stage_Data_Count, ASCENDING=False, remain=True)
        # Record previous Max_Correlation
        Pre_Max_Correlation = Max_Correlation
        # Update Max_Correlation
        Max_Correlation = df_cc_t.loc[0, 'corr_coef']
        # Record Max Correlation Recorded
        if Max_Correlation > Recorded_Max_Correlation:
            Recorded_Max_Correlation = Max_Correlation

        print (f"Maximal Correlation Coeffient : {Max_Correlation}")


        cc_ind_list       = df_cc_t['ind'].tolist()
        Stage_Lowest_Data = [Stage_Lowest_Data[i] for i in cc_ind_list]

    # Extract dominant genome
    if len(Stage_Lowest_Data) > Pre_set_Stage_Data_Count:
        Stage_Lowest_Data, _ = mu.Stage_Lowest_Data_find_dominant(Stage_Lowest_Data)
    else:
        pass

    # Save HIP Trace information

    Dict_HIP_Trace[Current_Mutation-1] = {'Current_Mutation':Current_Mutation,
                                          'Current_Human_PP':Current_Human_PP,
                                          'Pre_Human_PP':Pre_Human_PP,
                                          'Changed_PP':Current_Human_PP-Pre_Human_PP,
                                          'Current_Max_Human_PP':Current_Max_Human_PP,
                                          'Gap_to_Max_HIP':Pre_Human_PP - Current_Max_Human_PP,
                                          'Current_Max_Correlation':Max_Correlation,
                                          'Max_Correlation_Recorded':Recorded_Max_Correlation,
                                          'Pre_Max_Correlation':Pre_Max_Correlation,
                                          'Gap_to_Max_Correlation_Recorded':Recorded_Max_Correlation-Pre_Max_Correlation,
                                          'Target_PP':Target_PP_Threshold,
                                          'Start_PP':Start_HIP,
                                          'Total_Mutation_Run':len(Iter_List)}
    DF_HIP_Trace = pd.DataFrame.from_dict(Dict_HIP_Trace, 'index')
    ss.save_file(DF_HIP_Trace, f"{SAVE_FOLDER}//{Virus_Folder}//HIP_Trace")
    
    # Determine if terminate the simulation - terminate if there are improvements in HIP
    tpp_list = DF_HIP_Trace['Gap_to_Max_HIP'].tolist()[(-Termination_Iteration_n):]
    tpp_list.sort(reverse=True)
    if len(tpp_list) >= Termination_Iteration_n:
        if tpp_list[0] < 0:
            TERMINATION_ACT = True
            print (f"TERMINATION_ACT = {TERMINATION_ACT}")
    # Don't terminate it if there are improvements in correlation to COVID
    if TERMINATION_ACT is True:
        tmc_list = DF_HIP_Trace['Gap_to_Max_Correlation_Recorded'].tolist()[(-Termination_Iteration_n):]
        tmc_list.sort(reverse=True) 
        if tmc_list[0] <= 0:
            TERMINATION_ACT = False

    # Increase Current_Mutation
    Current_Mutation += 1


# Create iteration folder
Current_Mutation_Folder = f"Current_Mutation_[{Current_Mutation}]"
ss.create_folder(f"{SAVE_FOLDER}//{Virus_Folder}//Raw_data//{Current_Mutation_Folder}")

ss.save_file(df_show, f"{SAVE_FOLDER}//{Virus_Folder}//Raw_data//{Current_Mutation_Folder}//[{Current_Mutation}]_MUTATION_RES")
ss.save_file(df_cc_t, f"{SAVE_FOLDER}//{Virus_Folder}//Raw_data//{Current_Mutation_Folder}//[{Current_Mutation}]_Target_Corr_Select_RES")

Dict_HIP_Trace[Current_Mutation-1] = {'Current_Mutation':Current_Mutation,
                                        'Current_Human_PP':Current_Human_PP,
                                        'Pre_Human_PP':Pre_Human_PP,
                                        'Changed_PP':Current_Human_PP-Pre_Human_PP,
                                        'Current_Max_Human_PP':Current_Max_Human_PP,
                                        'Gap_to_Max_HIP':Pre_Human_PP - Current_Max_Human_PP,
                                        'Current_Max_Correlation':Max_Correlation,
                                        'Max_Correlation_Recorded':Recorded_Max_Correlation,
                                        'Pre_Max_Correlation':Pre_Max_Correlation,
                                        'Gap_to_Max_Correlation_Recorded':Recorded_Max_Correlation-Pre_Max_Correlation,
                                        'Target_PP':Target_PP_Threshold,
                                        'Start_PP':Start_HIP,
                                        'Total_Mutation_Run':len(Iter_List)}
DF_HIP_Trace = pd.DataFrame.from_dict(Dict_HIP_Trace, 'index')
ss.save_file(DF_HIP_Trace, f"{SAVE_FOLDER}//{Virus_Folder}//HIP_Trace")

if TERMINATION_ACT is True:
    print (f"### Simulation Terminated! - Target_Virus: {Stop_Virus_Name} - Final_Human_Pred_Proba: {Current_Human_PP} - Start_Human_Pred_Proba: {Start_HIP}")
    print (f"### No improvement in #{Termination_Iteration_n} iterations.")
else:
    print (f"### Simulation Done! - Target_Virus: {Stop_Virus_Name} - Final_Human_Pred_Proba: {Current_Human_PP} - Start_Human_Pred_Proba: {Start_HIP}")
  


# Linux commands
# nohup python Mutation_Path_Simulation_Directional.py > Mutation_Path_Simulation_Directional.log 2>&1 &


