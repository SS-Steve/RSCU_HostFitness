import pandas as pd
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
import os 
import matplotlib.pyplot as plt

import src.ss_database as db

author_name  = 'Shuquan(Steve)Su'                       
Entrez.email = "stevensu1991@gmail.com"
NCBI_API_KEY = '8ad0809ddce384eb810757af76eb9abdce09'
slash_mark = '//'

def PARAMS_CAP(**params):
    '''
    ### Capitalise all the input params 
    '''
    PARAMS_OUT = {}
    for param in params:
        PARAMS_OUT[param.upper()] = params[param]
    return PARAMS_OUT


def PARAMS(KEYWORD_LIST, DEFAULT=None, **params):
    '''
    ### Set input params
    # KEYWORD_LIST: a list of keyword searching within input params
    # DEFAULT: Default value if no keywords are found in parmas
    '''
    OUTPUT = DEFAULT
    for keyword in KEYWORD_LIST:
        if keyword in params:
            OUTPUT = params[keyword]
            break
    return OUTPUT


def create_folder(folder_name,option='Output',folder_path=''):
### Create a folder or dictionary
    system_required_folder, output_folder= '', 'Output'
    if option == 'Output':
        folder = output_folder
    elif option == 'System':
        folder = system_required_folder
    elif option == 'Other':
        folder = folder_path
    if isinstance(folder_name,str):
        create_foldername = folder + slash_mark + folder_name
        try:
            os.mkdir(create_foldername)
        except:
            pass
    elif isinstance(folder_name,list):
        for i in folder_name:
            create_foldername = folder + slash_mark + i
            try:
                os.mkdir(create_foldername)
            except:
                pass


def save_file(file, name=None, type0='csv', **params):
### save file, name = output file name, type0 = output file type
    params      = PARAMS_CAP(**params)
    dpi_value   = PARAMS(['DPI'], 1000, **params) 
    option      = PARAMS(['OPTION', 'OPT'], None, **params) 
    system_required_folder, output_folder, database_folder  = '', 'Output', ''
    if name is None:
        name = r'1_newly_saved'
    save_type_list = ['csv','fig','excel','txt']
    if type0 in save_type_list:
        # time_check = time.strftime("%Y-%m-%d[%Hh%Mm%Ss]", time.localtime())
        if option == 'time':
            # filename = output_folder + slash_mark + name + '__' + time_check + '_' + author_name
            pass
        else:
            filename = output_folder + slash_mark + name
        if type0 == 'csv':
            filename = filename + '.csv'
            file.to_csv(filename,index = False)
        elif type0 == 'excel':
            filename = filename + '.xlsx'
            if type(file) == dict:
                writer = pd.ExcelWriter(filename)
                for l in file:
                    file[l].to_excel(writer,sheet_name=str(l))
                writer.save()
            else:
                file.to_excel(writer,sheet_name='output')
        elif type0 == 'fig':
            try:
                fig_save = fig.get_figure()
            except:
                try: 
                    fig_save = file.get_figure()
                except:
                    fig_save = file
            try:
                fig_save.savefig(filename,dpi=dpi_value,bbox_inches='tight',transparent=True)
            except:
                plt.savefig(filename,dpi=dpi_value,bbox_inches='tight',transparent=True)
        elif type0 == 'txt':
            if isinstance(file, list) == True:
                file_t = ''
                for i in file:
                    file_t = file_t + i + '\n'
                file = file_t
            filename = filename + '.txt'
            with open(filename, 'w') as f:
                f.write(file)
                f.close()
    else:
        print ('wrong save type.')
        print ('support format: ' + str(save_type_list))


def database_load(FILENAME='', OPTION=None):
    '''
    ### Load different types of variables in database
    # i.e. model,dict, list, str...
    '''
    system_required_folder, output_folder, database_folder  = '', 'Output', ''
    import pickle 
    if   OPTION == 'System':
        filename = system_required_folder + slash_mark + FILENAME + '.pickle'
    elif OPTION == 'Output':
        filename = output_folder + slash_mark + FILENAME + '.pickle'
    elif OPTION == 'Custom':
        filename = FILENAME + '.pickle'
    elif OPTION is None:
        filename = database_folder + slash_mark + FILENAME + '.pickle'
    else:
        if isinstance(OPTION, str) == True:
            if OPTION[-1] == slash_mark:
                filename = OPTION + FILENAME + '.pickle'
            else:
                filename = OPTION + slash_mark + FILENAME + '.pickle'
        else:
            filename = None
    if filename is not None:
        with open(filename,'rb') as f:
            OUTPUT = pickle.load(f)
        f.close()
    else:
        OUTPUT = None
    return OUTPUT     


def database_save(file, name='', option=''):
    '''
    ### Save different types of variables in database
    # i.e. model,dict, list, str...
    '''
    system_required_folder, output_folder, database_folder  = '', 'Output', ''
    import pickle
    if option == 'System':
        filename = system_required_folder + slash_mark + name + '.pickle'
    elif option == 'Output':
        filename = output_folder + slash_mark + name + '.pickle'
    elif option == 'Custom':
        filename = name + '.pickle'
    else:
        filename = database_folder + slash_mark + name + '.pickle'
    with open(filename,'wb') as f:
        pickle.dump(file,f)
    f.close()


def mean_std(in0, raw_output=False):
    '''
    ### Calculate multiple descriptive data
    '''
    list_data = []
    dict_describe = {}
    if isinstance(in0, dict) == True:
        list_data = list(in0.values())
    elif isinstance(in0, list) == True:
        list_data = in0[:]
    else:
        try:
            list_data = list(in0)
        except:
            print ('not supported input.')
            pass
    if len(list_data) == 0:
        dict_describe = {'mean':None, 'std':None}
    else:
        dict_describe['mean']             = np.mean(list_data)
        dict_describe['std']              = np.std(list_data)
    if raw_output is True:
        dict_describe['raw_data'] = list_data
    return dict_describe


def list_to_str(List_In=[], Linking_Element=','):
    '''
    ### Convert a list to a string
    '''
    if List_In == []:
        OUTPUT = '[]'
    else:
        OUTPUT = '['
        for element in List_In:
            OUTPUT += f"{element}{Linking_Element}"
        OUTPUT = OUTPUT[:-1]
        OUTPUT += ']'
    return OUTPUT


def list_element_count(LIST_IN=[]):
    '''
    ### Count the number of every element in a list
    '''
    OUTPUT = {}
    for element in LIST_IN:
        if element in OUTPUT:
            OUTPUT[element] += 1
        else:
            OUTPUT[element]  = 1
    return OUTPUT


def dict_extract_list(dict0, list0='', option='keys'):
    '''
    ### Extract dict items according to list of keys
    # dict0: input dictionary
    # list0: list of keys
    '''
    if list0 == '':
        list0 = list(dict0.keys())
    dict_out = {}
    if option == 'keys':
        for i in list0:
            dict_out[i] = dict0[i]
    else:
        dict_out = dict0.copy()
    return dict_out


def load_RTG_Data(Select_Keyword_Input=None):
    '''
    ### Load the RTG data in folder 'Input//RTG_codon_data//'
    '''
    if isinstance(Select_Keyword_Input, str):
        select = Select_Keyword_Input.upper()
    else:
        raise ValueError('Input should be a string.')
    Input_Folder = 'Input//RTG_Data//ALL_VIRUS//'
    if   select in ['CU']:
        Output = pd.read_csv(f"{Input_Folder}//CU.csv")
    elif select in ['AAU']:
        Output = pd.read_csv(f"{Input_Folder}//AAU.csv")
    elif select in ['RSCU']:
        Output = pd.read_csv(f"{Input_Folder}//RSCU.csv")
    elif select in ['CDS_LENGTH', 'CDSL', 'CDS_L', 'CDS']:
        Output = pd.read_csv(f"{Input_Folder}//CDS_LENGTH.csv")
    elif select in ['ATGC']:
        Output = pd.read_csv(f"{Input_Folder}//ATGC.csv")
    elif select in ['HOST']:
        Output = pd.read_csv(f"{Input_Folder}//HOST.csv")
    elif select in ['INFO']:
        Output = pd.read_csv(f"{Input_Folder}//INFO.csv")
    elif select in ['PARTITE']:
        Output = pd.read_csv(f"{Input_Folder}//PARTITE.csv")
    elif select in ['TAXONOMY', 'TAX']:
        Output = pd.read_csv(f"{Input_Folder}//TAXONOMY.csv")
    elif select in ['START_STOP']:
        Output = pd.read_csv(f"{Input_Folder}//START_STOP.csv")
    elif select in ['HS_CORR']:
        Output = pd.read_csv(f"{Input_Folder}//HS_CORR.csv")
    elif select in ['HS_CORR_AA']:
        Output = pd.read_csv(f"{Input_Folder}//HS_CORR_AA.csv")
    elif select in ['TAXONOMY_RAW', 'TAX_RAW']:
        Output = pd.read_csv(f"{Input_Folder}//TAXONOMY_RAW.csv")
    return Output


def translate(seq=''):
    from Bio.Seq import Seq
    if (len(seq) % 3) == 0:
        seq = Seq(seq)
        product = str(seq.translate())
    else:
        logging.warning('Cannot divide by 3.')
        product = None
    return product


def get_all_cds_detail(virus_genome_id = '', option = '', **params):  
    '''
    ### Obtain a combined CDS from NCBI according to accession#
    # INPUT: single NCBI ID 
    # OUTPUT: String of combined all CDS of virus genome
    '''
    params  = PARAMS_CAP(**params)
    WARNING = PARAMS(['WARN', 'WARNING', 'W'], True, **params)
    dict_output = {}
    try:
        hd1 = Entrez.efetch(db      = 'nucleotide', 
                            id      = virus_genome_id, 
                            rettype = 'gb', 
                            api_key = NCBI_API_KEY)
        seq      = SeqIO.read(hd1,'gb')
        features = seq.features
        for feature in features:
            if feature.type == "CDS":
                cds_info_dict = {'type'             : 'CDS', 
                                 'location'         : str(feature.location),
                                 'location_extract' : feature.location}
                for i in feature.qualifiers:
                    if len(feature.qualifiers[i]) == 1:
                        cds_info_dict[i] = str(feature.qualifiers[i][0])
                    else:
                        cds_info_dict[i] = feature.qualifiers[i]    
                cds = SeqRecord(feature.location.extract(seq).seq)
                cds_info_dict['mRNA_sequence'] = str(cds.seq)
                try:
                    dict_output[cds_info_dict['gene']] = cds_info_dict
                except:
                    dict_output['protein:' + cds_info_dict['protein_id']] = cds_info_dict
    except:
        if WARNING == True:
            print ('-Error occured on ' + virus_genome_id)
            print ('--CANNOT download sequence from internet.')
        pass
    if option == 'df':
        df0 = pd.DataFrame()
        c   = 0
        for i in dict_output:
            df0.loc[c,'ID'] = i
            for ii in dict_output[i]:
                if type(dict_output[i][ii]) == str:
                    df0.loc[c,ii] = dict_output[i][ii]
                else:
                    df0.loc[c,ii] = str(dict_output[i][ii])
            c += 1
        dict_output = df0
    return dict_output


def get_all_cds(virus_genome_id): 
    dict_output  = {}
    dict_input   = get_all_cds_detail(virus_genome_id)
    for i in dict_input:
        dict_output[i] = dict_input[i]['mRNA_sequence']
    return dict_output


def get_whole_seq(ncbi_id = '', **params):
    '''
    ### Download a complete sequence according to NCBI Accession ID
    '''
    params  = PARAMS_CAP(**params)
    WARNING = PARAMS(['WARNING', 'WARN', 'W'], True, **params)
    dict_out = {}
    try:
        hd1 = Entrez.efetch(db = 'nucleotide',
                            id = ncbi_id,
                            rettype = 'gb', 
                            api_key=NCBI_API_KEY)
        seq = SeqIO.read(hd1,'gb')
        seq_out = str(seq.seq)
        dict_out[ncbi_id] = seq_out
    except:
        if WARNING == True:
            print ('-Error occured on ' + ncbi_id)
            print ('--CANNOT download sequence from internet.')
        pass 
    return dict_out



def codon_dict(GS='',choice=''):
    codon_dict = {}
    if choice == 'pair':
        optional_length = 3
    elif type(choice) == int:
        optional_length = choice
    else:
        optional_length = 0
    GSL = len (GS)
    GSLtest0 = GSL / 3
    GSLtest = GSL / 3
    position = 1
    if GSLtest - int(GSLtest) == 0:
        while GSLtest > 0:
            codon = GS[(position - 1) * 3:position * 3 + optional_length]
            codon_dict[position] = codon
            GSLtest = GSLtest - 1 
            position = position + 1
    else:
        print ('Cannot divide by 3.')
    if choice == 'pair':
        del codon_dict[len(codon_dict)]
    return codon_dict


def codon_list(GS='',choice=''):
    codon_d = codon_dict(GS,choice)
    codon_list = list(codon_d.values())
    return codon_list


def codon_n_dict(CDS_Seq=''):
    '''
    ### Calculate codon number dict from a CDS Sequence
    '''
    codon_table_list = db.codon_table_list[:]
    cl = codon_list(CDS_Seq)
    Codon_N_Dict = {}
    for codon in codon_table_list:
        Codon_N_Dict[codon] = cl.count(codon)
    return Codon_N_Dict


def codon_n_dict_to_aa_n_dict(Codon_N_Dict={}):
    '''
    ### Generate a AA_N_Dict from a Codon_N_Dict
    '''
    translation_aa_codon_dict = db.translation_aa_codon_dict.copy()
    AA_N_Dict = {}
    for aa, aa_codon_list in translation_aa_codon_dict.items():
        sum_n = 0
        for codon in aa_codon_list:
            sum_n += Codon_N_Dict[codon]
        AA_N_Dict[aa] = sum_n
    return AA_N_Dict


def codon_n_dict_to_RSCU(Codon_N_Dict):
    '''
    ### Convert a dicitonary of Codon_N_Dict into a RSCU dictionary
    '''
    codon_codon_box_dict_numeric = db.codon_codon_box_dict_numeric.copy()
    AA_N_Dict = codon_n_dict_to_aa_n_dict(Codon_N_Dict)
    RSCU = {}
    for codon, codon_count in Codon_N_Dict.items():
        aa_count = AA_N_Dict[translate(codon)]
        if aa_count == 0:
            RSCU[codon] = 0
        else:
            rscu = (codon_count / aa_count) * codon_codon_box_dict_numeric[codon]
            RSCU[codon] = rscu
    return RSCU


def codon_n_dict_mutation(Codon_N_Dict, Mutation, Warning_Opt=True):
    '''
    ### Alter the dictionary of the Codon_N_Dict with Mutation
    '''
    if   isinstance(Mutation, str) is True:
        Mutation_List = [Mutation]
    elif isinstance(Mutation, list) is True:
        Mutation_List = Mutation
    else:
        if Warning_Opt is True:
            print("Wrong type of Mutation.")
        Mutation_List = []
    new_Codon_N_Dict = Codon_N_Dict.copy()
    for Mut in Mutation_List:
        if '->' in Mut:
            Pre, Post = Mut.split('->')
            if new_Codon_N_Dict[Pre] > 0:
                new_Codon_N_Dict[Pre]  -= 1
                new_Codon_N_Dict[Post] += 1
            else:
                new_Codon_N_Dict = False
                if Warning_Opt is True:
                    print(f"Not sufficient codon '{Pre}' for mutation {Mut}.")
                pass
        else:
            if   '-' in Mut:
                Pre = Mut.replace('-', '')
                if new_Codon_N_Dict[Pre] > 0:
                    new_Codon_N_Dict[Pre]  -= 1
                else:
                    new_Codon_N_Dict = False
                    if Warning_Opt is True:
                        print(f"Not sufficient codon '{Pre}' for mutation {Mut}.")
                    pass
            elif '+' in Mut:
                Post = Mut.replace('+', '')
                new_Codon_N_Dict[Post] += 1
            else:
                pass
    return new_Codon_N_Dict


def extract_df_row_top_n(DF_In, Column='', Top_N=100, **params):
    '''
    ### Extract dataframe rows according to a Columns and Top-N value
    '''
    params            = PARAMS_CAP(**params)
    ASCEND_OPT        = PARAMS(['ASCENDING', 'ASCEND', 'A'], True, **params)
    REMAIN_OUTPUT_OPT = PARAMS(['REMAIN', 'REMAIN_OUT'], False, **params)
    df_in = DF_In.sort_values(Column, ascending=ASCEND_OPT, ignore_index=True)
    Threshold = df_in.loc[Top_N-1, Column]
    if   ASCEND_OPT is True:
        OUTPUT_T = df_in[df_in[Column] <= Threshold]
        OUTPUT_R = df_in[df_in[Column] > Threshold]
    elif ASCEND_OPT is False:
        OUTPUT_T = df_in[df_in[Column] >= Threshold]
        OUTPUT_R = df_in[df_in[Column] < Threshold]
    else:
        OUTPUT_T, OUTPUT_R = None, None
    if REMAIN_OUTPUT_OPT is True:
        OUTPUT_T = OUTPUT_T.reset_index(drop=True)
        OUTPUT_R = OUTPUT_R.reset_index(drop=True)
        OUTPUT = [OUTPUT_T, OUTPUT_R]
    else:
        OUTPUT = OUTPUT_T.reset_index(drop=True)
    return OUTPUT


def in_range(num, range_in):
    '''
    ### Determine if the number is in a range
    # num: input number
    # range_in: a list contain two numbers (The first number should be smaller) - string of number means exclude that number
    '''
    ans = ''
    if len(range_in) == 2:
        if type(range_in[0]) == str and type(range_in[1]) != str: 
            if num > float(range_in[0]) and num <= range_in[1]:
                ans = True
            else:
                ans = False
        elif type(range_in[0]) != str and type(range_in[1]) == str: 
            if num >= range_in[0] and num < float(range_in[1]):
                ans = True
            else:
                ans = False
        elif type(range_in[0]) == str and type(range_in[1]) == str: 
            if num > float(range_in[0]) and num < float(range_in[1]):
                ans = True
            else:
                ans = False
        else:
            if num >= range_in[0] and num <= range_in[1]:
                ans = True
            else: 
                ans = False                
    else:
        print('wrong input')
    return ans


def determine_in_range(in0, determine_dict={}):
    '''
    ### Determine which range is the input number in
    # in0: input number
    # determine_dict: a dictionary of different input range
    '''
    ans = ''
    for i in determine_dict:
        if in_range(in0, determine_dict[i]) == True:
            ans = i
    return ans


def corr_coef_pearson(array1,array2,option=''):
    '''
    ### calculate Pearson correlation coeffecient - sample must follow normal distribution
    # array1: a list or a numpy array
    # array2: another list or numpy array
    '''
    correlation_coefficient_determine = {'Significant':[0.95 ,1     ], 
                                        'High':       [0.5  ,'0.95'], 
                                        'Moderate':   [0.3  ,'0.5' ], 
                                        'Low':        [0.05 ,'0.3' ],
                                        'No':         [0    ,'0.05']}
    if len(array1) == len(array2):
        corr_coef = np.corrcoef(array1, array2)[1,0]
        corr_coef = float(corr_coef)
        if option == 'all':
            description = determine_in_range(corr_coef,
                                             correlation_coefficient_determine)
            output = {'corr_coef'   :corr_coef,
                      'description' :description}
        else:
            output = corr_coef
    else:
        if option == 'all':
            output = {'corr_coef'   :np.nan,
                      'description' :'unequal_length'}
        else:
            output = np.nan
    return output




