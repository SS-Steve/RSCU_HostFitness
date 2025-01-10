import numpy as np
import pandas as pd
import src.ss_utils as ss
import src.ss_database as db


def mutate_len_dict(Current_LEN_dict, Mutation, Current_CDS_LEN_dict=None):
    '''
    ### Change the len dict regarding the mutations
    '''
    OUTPUT_LEN_dict = Current_LEN_dict.copy()
    if Current_CDS_LEN_dict is not None:
        OUTPUT_CDS_LEN_dict = Current_CDS_LEN_dict.copy()
    if '@' in Mutation:
        Mutation, Mutation_Loc = Mutation.split('@')
    else:
        Mutation_Loc = None
    if '->' in Mutation:
        pass
    else:
        if   '-' in Mutation:
            OUTPUT_LEN_dict['genome_length']    -= 3
            OUTPUT_LEN_dict['contenCDS_length'] -= 3
            OUTPUT_LEN_dict['cds_length_mean']  = OUTPUT_LEN_dict['contenCDS_length'] / OUTPUT_LEN_dict['cds_count']
            if Mutation_Loc is not None: # Calculate the CDS len std according to mutation location
                OUTPUT_CDS_LEN_dict[Mutation_Loc] -= 3
                OUTPUT_LEN_dict['cds_length_std'] = np.std(list(OUTPUT_CDS_LEN_dict.values()))
        elif '+' in Mutation:
            OUTPUT_LEN_dict['genome_length']    += 3
            OUTPUT_LEN_dict['contenCDS_length'] += 3
            OUTPUT_LEN_dict['cds_length_mean']  = OUTPUT_LEN_dict['contenCDS_length'] / OUTPUT_LEN_dict['cds_count']
            if Mutation_Loc is not None: # Calculate the CDS len std according to mutation location
                OUTPUT_CDS_LEN_dict[Mutation_Loc] += 3
                OUTPUT_LEN_dict['cds_length_std'] = np.std(list(OUTPUT_CDS_LEN_dict.values()))
        else:
            pass
    if Current_CDS_LEN_dict is not None:
        OUTPUT = [OUTPUT_LEN_dict, OUTPUT_CDS_LEN_dict]
    else:
        OUTPUT = OUTPUT_LEN_dict.copy()
    return OUTPUT


def mutate_codon_n_dict_single(Current_Codon_N_Dict_Single, Mut, Loc):
    '''
    ### Change the codon_n_dict according to mutations
    '''
    new_Codon_N_Dict_Single = Current_Codon_N_Dict_Single.copy()
    if '-' in Mut and '>' not in Mut:
        pre = Mut.replace('-', '')
        if new_Codon_N_Dict_Single[Loc][pre] == 0:
            OUTPUT = False
        else:
            new_Codon_N_Dict_Single[Loc][pre] -= 1
            OUTPUT = new_Codon_N_Dict_Single.copy()
    else:
        OUTPUT = new_Codon_N_Dict_Single.copy()
    return OUTPUT


def mutation_dict_add_mutation(Codon_Change_Dict, Mut):
    '''
    ### Change the codon_change_dict according to mutations
    '''
    output = Codon_Change_Dict.copy()
    if '>' in Mut:
        pre, post = Mut.split('->')
        if pre in output:
            output[pre] -= 1
        else:
            output[pre] = -1
        if post in output:
            output[post] += 1
        else:
            output[post] = 1
    else:
        if '@' in Mut:
            mut = Mut.split('@')[0]
        else:
            mut = Mut
        if   '-' in mut:
            pre = mut.replace('-', '')
            if pre in output:
                output[pre] -= 1
            else:
                output[pre] = -1
        elif '+' in mut:
            post = mut.replace('+', '')
            if post in output:
                output[post] += 1
            else:
                output[post] = 1
    return output


def Reverse_Mut_Filtering(Codon_Change_Dict, Mut, Reverse_Filtering_Decision=True):
    '''
    ### Filter the reverse mutations
    '''
    Reverse_Filtering = False
    if Reverse_Filtering_Decision is False:
        pass
    else:
        if '@' in Mut:
            mut = Mut.split('@')[0]
        else:
            mut = Mut
        if '->' in mut:
            pre, post = mut.split('->')
        else:
            if '-' in mut:
                pre = mut.split('-')[1]
                post = None
            elif '+' in mut:
                pre = None
                post = mut.split('+')[1]
        if pre is not None and Codon_Change_Dict[pre] > 0:
            Reverse_Filtering = True
        if post is not None and Codon_Change_Dict[post] < 0:
            Reverse_Filtering = True
    return Reverse_Filtering


def single_cds_len_dict_to_COVID_LEN_X(Single_CDS_Len_Dict, start_whole_seq, start_all_cds_dict):
    '''
    ### Convert the Single_CDS_Len_Dict to CDS_Len_Dict for X (dataset)
    '''
    mean_std = ss.mean_std(list(Single_CDS_Len_Dict.values()))
    COVID_LEN_X_dict = {'genome_length'   :len(start_whole_seq),
                        'contenCDS_length':sum(list(Single_CDS_Len_Dict.values())),
                        'cds_count'       :len(start_all_cds_dict),
                        'cds_length_mean' :mean_std['mean'],
                        'cds_length_std'  :mean_std['std']}
    return COVID_LEN_X_dict


def X_RSCU_to_RSCU_dict(X_RSCU):
    '''
    ### Convert a X dataset of RSCU for a dictionary - later used to calculate correlation coefficient
    '''
    d_t = X_RSCU.to_dict('list')
    Output_Dict = {}
    for codon in db.multibox_codon_table_list:
        Output_Dict[codon] = d_t[f"[RSCU]_{codon}"][0]
    return Output_Dict


def Stage_Lowest_Data_find_dominant(Stage_Lowest_Data):
    list_of_search = []
    data_dict_pool = {}
    for data_dict in Stage_Lowest_Data:
        codon_n_dict_str = str(data_dict['Codon_N_Dict'])
        list_of_search.append(codon_n_dict_str)
        if   codon_n_dict_str not in data_dict_pool:
            data_dict_pool[codon_n_dict_str] = [data_dict]
        else:
            data_dict_pool[codon_n_dict_str].append(data_dict)
    lec = ss.list_element_count(list_of_search)
    max_count = max(list(lec.values()))
    select_list = []
    for k, v in lec.items():
        if v == max_count:
            select_list.append(k)
    OUTPUT_LIST = []
    d_mut_path = {}
    for dict_ind, i_str in enumerate(select_list):
        OUTPUT_LIST.append(data_dict_pool[i_str][0]) # Take the first one of the same sequence
        list_t = []
        for data_dict in data_dict_pool[i_str]:
            list_t.append(data_dict['pre_mutation'])
        d_mut_path[dict_ind] = {'Codon_N_Dict':data_dict['Codon_N_Dict'],
                                'pre_mutation':list_t}
    return OUTPUT_LIST, d_mut_path


# def Find_first_mutation_path(Mutation_Num):
#     df_in = ss.input_file(f'33//Each_Step_Data//[{Mutation_Num}]_MUTATION_RES', loc='Output')
#     mut_list = df_in.loc[0, 'pre_mutation']
#     mut_list = mut_list[1:-1].split(',')
#     res = {}
#     c = 0
#     for loc, mut in enumerate(mut_list):
#         if '->' in mut:
#             pre, post = mut.split('->')
#         else:
#             if '-' in mut:
#                 pre = mut[1:]
#                 post = None
#             elif '+' in mut:
#                 pre = None
#                 post = mut[1:]
#         res_t = {'index':loc+1,
#                  'pre':pre,
#                  'post':post}
#         res[c] = res_t
#         c += 1
#     df_res = pd.DataFrame.from_dict(res, 'index')
#     pre_count_dict  = ss.list_element_count(df_res['pre'].tolist())
#     post_count_dict = ss.list_element_count(df_res['post'].tolist())
#     for k, v in pre_count_dict.items():
#         if k in post_count_dict:
#             post_count_dict[k] -= v
#         else:
#             post_count_dict[k] = -v
#     return post_count_dict


empty_Codon_Change_Dict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'TCT': 0, 
                           'TCC': 0, 'TCA': 0, 'TCG': 0, 'TAT': 0, 'TAC': 0, 
                           'TAA': 0, 'TAG': 0, 'TGT': 0, 'TGC': 0, 'TGA': 0, 
                           'TGG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 
                           'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'CAT': 0, 
                           'CAC': 0, 'CAA': 0, 'CAG': 0, 'CGT': 0, 'CGC': 0, 
                           'CGA': 0, 'CGG': 0, 'ATT': 0, 'ATC': 0, 'ATA': 0, 
                           'ATG': 0, 'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 
                           'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'AGT': 0, 
                           'AGC': 0, 'AGA': 0, 'AGG': 0, 'GTT': 0, 'GTC': 0, 
                           'GTA': 0, 'GTG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 
                           'GCG': 0, 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0, 
                           'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

Virus_Accession_Dict = {'Severe acute respiratory syndrome coronavirus 2'      : 'NC_045512',
                        'Bovine coronavirus'                                   : 'NC_003045',
                        'Human coronavirus OC43'                               : 'NC_006213',
                        'Middle East respiratory syndrome-related coronavirus' : 'NC_019843',
                        'Betacoronavirus England 1'                            : 'NC_038294',
                        'SARS coronavirus Tor2'                                : 'NC_004718',
                        'Human coronavirus HKU1'                               : 'NC_006577',
                        'Tylonycteris bat coronavirus HKU4'                    : 'NC_009019',
                        'Rabbit coronavirus HKU14'                             : 'NC_017083',
                        'Bat coronavirus BM48-31/BGR/2008'                     : 'NC_014470',
                        'Betacoronavirus Erinaceus/VMC/DEU/2012'               : 'NC_039207',
                        'Bat Hp-betacoronavirus/Zhejiang2013'                  : 'NC_025217',
                        'Pipistrellus bat coronavirus HKU5'                    : 'NC_009020',
                        'Rousettus bat coronavirus'                            : 'NC_030886',
                        'Betacoronavirus HKU24'                                : 'NC_026011',
                        'Rat coronavirus Parker'                               : 'NC_012936',
                        'Murine hepatitis virus'                               :['NC_001846', 'NC_048217'],
                        'Rousettus bat coronavirus HKU9'                       : 'NC_009021'}


Max_n_mutation_Ref_Dict = {'Bat coronavirus BM48-31_BGR_2008'   : 979,
                           'Bat Hp-betacoronavirus_Zhejiang2013': 1004,
                           'Betacoronavirus HKU24'              : 858,
                           'Pipistrellus bat coronavirus HKU5'  : 1270,
                           'Rabbit coronavirus HKU14'           : 1153,
                           'Rat coronavirus Parker'             : 2360,
                           'Rousettus bat coronavirus'          : 1659,
                           'Rousettus bat coronavirus HKU9'     : 1693}

