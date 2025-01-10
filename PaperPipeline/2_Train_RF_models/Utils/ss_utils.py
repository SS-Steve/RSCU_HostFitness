### Import basic modules

import numpy as np
import pandas as pd
from scipy import stats
import time
import os
import matplotlib.pyplot as plt

from Utils import ss_database as db

slash_mark = "//"
input_folder = ""
output_folder = ""



########################################################################################################################
### General functions ###
########################################################################################################################


def PARAMS_CAP(**params):
    """
    ### Capitalise all the input params
    """
    PARAMS_OUT = {}
    for param in params:
        PARAMS_OUT[param.upper()] = params[param]
    return PARAMS_OUT


def PARAMS(KEYWORD_LIST, DEFAULT=None, **params):
    """
    ### Set input params
    # KEYWORD_LIST: a list of keyword searching within input params
    # DEFAULT: Default value if no keywords are found in parmas
    """
    OUTPUT = DEFAULT
    for keyword in KEYWORD_LIST:
        if keyword in params:
            OUTPUT = params[keyword]
            break
    return OUTPUT


def input_file(FILENAME="", FILE_TYPE=None, **params):
    """
    ### input file, name = input file name (no suffix), FILE_TYPE = input file type
    """
    params = PARAMS_CAP(**params)
    FILE_TYPE_IN = PARAMS(
        ["TYPE", "FILE_TYPE", "FILETYPE", "EXTENSION", "EX"], "PD_CSV", **params
    )
    LOCATION_IN = PARAMS(["LOCATION", "LOC"], "", **params)
    OPTION_IN = PARAMS(["OPTION", "OPT"], "DF", **params)
    WARNING = PARAMS(["WARNING", "WARN", "W"], True, **params)
    if isinstance(FILE_TYPE_IN, str) == True:
        FILE_TYPE_IN = FILE_TYPE_IN.upper()
    if FILE_TYPE is None:
        FILE_TYPE = FILE_TYPE_IN
    OUTPUT = ""
    if LOCATION_IN is None:
        filename = input_folder + slash_mark + FILENAME
    elif LOCATION_IN in ["path", "Path", "PATH", "p", "P"]:
        filename = FILENAME
    else:
        if LOCATION_IN[-1] == slash_mark:
            filename = LOCATION_IN + FILENAME
        else:
            filename = LOCATION_IN + slash_mark + FILENAME
    if FILE_TYPE in ["CSV"]:
        import csv

        filename = filename + ".csv"
        try:
            file0 = open(filename)
            OUTPUT = csv.reader(file0)
        except:
            if WARNING is True:
                print("not found.")
            pass
    elif FILE_TYPE in ["PD_CSV", "DF", "DATAFRAME"]:
        filename = filename + ".csv"
        try:
            OUTPUT = pd.read_csv(filename)
        except:
            if WARNING is True:
                print("not found.")
            pass
    elif FILE_TYPE in ["TXT", "TEXT"]:
        filename = filename + ".txt"
        try:
            with open(filename) as f:
                OUTPUT = f.readlines()
                f.close()
        except:
            if WARNING is True:
                print("not found.")
            pass
    elif FILE_TYPE in ["FCS", "FLOWJO", "FJ"]:
        import flowkit

        filename = filename + ".fcs"
        try:
            sample_load = flowkit.Sample(filename)
        except:
            pass
        if OPTION_IN is None or OPTION_IN in ["DF", "DATAFRAME"]:
            OUTPUT = sample_load.as_dataframe(source="raw")
        elif OPTION_IN in ["ARRAY", "A", "NP"]:
            OUTPUT = sample_load.get_events()
    elif FILE_TYPE in ["FNA", "SEQ"]:
        from Bio import SeqIO

        filename = filename + ".fna"
        fasta_sequences = SeqIO.parse(open(filename), "fasta")
        OUTPUT = []
        for fasta in fasta_sequences:
            OUTPUT.append(fasta)
    else:
        OUTPUT = None
    return OUTPUT


def save_file(file, name=None, type0="csv", **params):
    """
    ### save file, name = output file name, type0 = output file type
    """
    params = PARAMS_CAP(**params)
    dpi_value = PARAMS(["DPI"], 1000, **params)
    option = PARAMS(["OPTION", "OPT"], None, **params)
    if name is None:
        name = r"1_newly_saved"
    save_type_list = ["csv", "fig", "excel", "txt"]
    if type0 in save_type_list:
        time_check = time.strftime("%Y-%m-%d[%Hh%Mm%Ss]", time.localtime())
        if option == "time":
            filename = (
                output_folder
                + slash_mark
                + name
                + "__"
                + time_check)
        else:
            filename = name
            # filename = output_folder + slash_mark + name
        if type0 == "csv":
            filename = filename + ".csv"
            file.to_csv(filename, index=False)
        elif type0 == "excel":
            filename = filename + ".xlsx"
            if type(file) == dict:
                writer = pd.ExcelWriter(filename)
                for l in file:
                    file[l].to_excel(writer, sheet_name=str(l))
                writer.save()
            else:
                file.to_excel(writer, sheet_name="output")
        elif type0 == "fig":
            try:
                fig_save = fig.get_figure()
            except:
                try:
                    fig_save = file.get_figure()
                except:
                    fig_save = file
            try:
                fig_save.savefig(
                    filename, dpi=dpi_value, bbox_inches="tight", transparent=True
                )
            except:
                plt.savefig(
                    filename, dpi=dpi_value, bbox_inches="tight", transparent=True
                )
        elif type0 == "txt":
            if isinstance(file, list) == True:
                file_t = ""
                for i in file:
                    file_t = file_t + i + "\n"
                file = file_t
            filename = filename + ".txt"
            with open(filename, "w") as f:
                f.write(file)
                f.close()
    else:
        print("wrong save type.")
        print("support format: " + str(save_type_list))


def database_save(file, filename):
    """
    ### Save different types of variables in database
    # i.e. model, dict, list, str...
    """
    import pickle

    with open(f"{filename}.pickle", "wb") as f:
        pickle.dump(file, f)
    f.close()


def database_load(filename=""):
    """
    ### Load different types of variables in database
    # i.e. model, dict, list, str...
    """
    import pickle

    with open(f"{filename}.pickle", "rb") as f:
        OUTPUT = pickle.load(f)
    f.close()
    return OUTPUT


def create_folder(folder_name,option='Output',folder_path=''):
    '''
    ### Try to create a folder
    '''
    if option == 'Output':
        folder = output_folder
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


def dict_save_as_txt(DICT_IN={}, FILENAME='', **params):
    '''
    ### Save the dictionary as txt file (i.e. parameters dictionary)
    '''
    params = PARAMS_CAP(**params)
    LOC    = PARAMS(['LOC', 'LOCATION', 'L'], 'Output', **params)
    LOC    = LOC.upper()
    import logging
    if   isinstance(DICT_IN, dict) is True:
        list_run = [DICT_IN]
    elif isinstance(DICT_IN, list) is True:
        list_run = DICT_IN[:]
    else:
        list_run = []
        logging.warning('No Dictionary.')
    write_string = ''
    for Dictonary in list_run:
        # Find longest key lenght
        key_lenght = len(max(list(Dictonary.keys()), key=len))
        for k, v in Dictonary.items():
            write_string += f"{' '*(key_lenght-len(k))}{k} : {v} \n"
        write_string += '\n'
    if   LOC in ['OUTPUT', 'OUT']:
        FOLDERNAME = output_folder
    elif LOC in ['INPUT', 'IN']:
        FOLDERNAME = input_folder
    filename = f"{FOLDERNAME}{slash_mark}{FILENAME}" 
    if filename[:-4] != '.txt':
        filename += '.txt'
    if len(list_run) > 0:
        fo = open(filename, "w")
        fo.write(write_string)
        fo.close()


def del_head_space(str0, target=" "):
    """
    ### Delele space or other character at the begining of a string
    # str0: input string
    # target: characters you want to delete
    """
    str_out = ""
    if type(str0) == str:
        c = 0
        while c < len(str0):
            str_t = str0[c : c + 1]
            if str_t == target:
                c = c + 1
            else:
                break
        str_out = str0[c : len(str0)]
    return str_out


def list_remove_list(list0=[], list1=[], list2=[]):
    """
    ### Remove values in a list
    """
    list_main = list0[:]
    list_sub_1 = list1[:]
    list_sub_2 = list2[:]
    for i in list_sub_1:
        if i in list_main:
            list_main.remove(i)
    for ii in list_sub_2:
        if ii in list_main:
            list_main.remove(ii)
    return list_main


def merge_dict(dict1={}, dict2={}):
    """
    ### Merge 2 dict into 1 dict
    """
    merge_dict = dict1.copy()
    merge_dict.update(dict2)
    return merge_dict


def dict_head(dict0, num=3):
    """
    ### Output first few items of a dictionary
    # dict0: input dictionary
    # num: number of items you want to extract
    """
    dict_out = {}
    c = 0
    for i in dict0:
        dict_out[i] = dict0[i]
        c = c + 1
        if c >= num:
            break
    return dict_out


def dict_extract_list(dict0, list0=[]):
    """
    ### Extract dict items according to list of keys
    # dict0: input dictionary
    # list0: list of keys
    """
    dict_out = {}
    for i in list0:
        dict_out[i] = dict0[i]
    return dict_out


def dict_to_list(input_dict={}):
    """
    ### Convert dict into list - combine key and value into a str
    """
    output_list = []
    for i in input_dict:
        output_list.append(str(i) + ":" + str(input_dict[i]))
    return output_list


def transpose_df(df0, direction="h", choice=0):
    """
    ### Transpose dataframe
    """
    output_df = pd.DataFrame()
    if direction == "h":
        column_list = df0.columns.tolist()
        output_df["columns"] = column_list
        count = 0
        while count < df0.shape[0]:
            output_column_name = "row_" + str(count)
            output_df[output_column_name] = df0.loc[count, :].tolist()
            count = count + 1
    elif direction == "v":
        if choice == 0:
            column_list = df0.iloc[:, 0].tolist()
            count = 0
        elif choice == 1:
            column_list = list(np.arange(df0.shape[0]))
            count = 0
        for i in column_list:
            output_df[i] = df0.loc[count, :].tolist()
            count = count + 1
        if choice == 0:
            output_df = output_df.drop([0], axis=0)
    return output_df


def merge_df(df1, df2, option="h", column=""):
    """
    ### Merge two dataframe into one dataframes
    # df1: first input dataframe
    # df2: second input dataframe
    # option: option of how do you want to merge two daaframes - ['v','h','h_on']
    #       'v': 'vertically'
    #       'h': 'horizontally'
    #       'h_on': 'horizontally' on a column both dataframes contian
    # column: name string of a shared column both dataframes contain
    """
    if option == "v":
        df_out = df1.append(df2, ignore_index=True)
    elif option == "h_on":
        df_out = pd.merge(df1, df2, on=column)
    elif option == "h":
        df_out = pd.concat([df1, df2], join="outer", axis=1)
    return df_out


def df_to_dict(DF_IN, DF_DIRECTION="V", **params):
    """
    ### Convert a datafarme into a dict
    # DF_IN: input dataframe
    # DF_DIRECTION: Direction of converting ('Vertical' or 'Horizontal')
    # ID_COL= or ID=: column name or row index for dictionary keys (Default: None)
    """
    params = PARAMS_CAP(**params)
    ID_COL = PARAMS(
        ["ID", "ID_COL", "ID_COLUMN", "IDCOL", "IDCOLUMN", "ID_ROW", "IDROW"],
        None,
        **params,
    )
    df_in = DF_IN.copy()
    DICT_OUTPUT = {}
    DF_DIRECTION = DF_DIRECTION.upper()
    if DF_DIRECTION in ["V", "VE", "VER", "VERTICAL"]:
        column_list = df_in.columns.tolist()
        if ID_COL == None or ID_COL not in column_list:
            df_in = df_in.reset_index(drop=True)
            index_list = df_in.index.tolist()
            for index in index_list:
                dict_t = {}
                for column in column_list:
                    dict_t[column] = df_in.loc[index, column]
                DICT_OUTPUT[index + 1] = dict_t
        else:
            key_list = df_in[ID_COL].tolist()
            df_in = df_in.reset_index(drop=True)
            index_list = df_in.index.tolist()
            column_list.remove(ID_COL)
            del df_in[ID_COL]
            for index in index_list:
                dict_t = {}
                for column in column_list:
                    dict_t[column] = df_in.loc[index, column]
                DICT_OUTPUT[key_list[index]] = dict_t
    elif DF_DIRECTION in ["H", "HO", "HOR", "HORIZONTAL"]:
        index_list = df_in.index.tolist()
        column_list = df_in.columns.tolist()
        if ID_COL == None or ID_COL not in index_list:
            for column in column_list:
                dict_t = {}
                for index in index_list:
                    dict_t[index] = df_in.loc[index, column]
                DICT_OUTPUT[column] = dict_t
        else:
            key_list = df_in.iloc[ID_COL].tolist()
            index_list.remove(ID_COL)
            c = 0
            for column in column_list:
                dict_t = {}
                for index in index_list:
                    dict_t[index] = df_in.loc[index, column]
                DICT_OUTPUT[key_list[c]] = dict_t
                c += 1
    return DICT_OUTPUT


def dict_to_df(DICT_INPUT={}, DIRECTION=None, **params):
    """
    ### Convert dict into dataframe
    # DIRECTION or D=: indicator for direction of output dataframe ['H' or 'V']
    # COL=: column option for output dataframe
    """
    params = PARAMS_CAP(**params)
    DIRECTION_IN = PARAMS(["DIRECTION", "D", "DIRECT", "DIR"], None, **params)
    COL_OPTION_IN = PARAMS(
        ["OPTION", "O", "OPT", "COL", "COLUMN", "COLUMNS"], "KEYS", **params
    )
    if DIRECTION == None:
        if DIRECTION_IN == None:
            DIRECTION = "H"
        else:
            DIRECTION = DIRECTION_IN
    if isinstance(DIRECTION, str) == True:
        DIRECTION = DIRECTION.upper()
    if isinstance(COL_OPTION_IN, str) == True:
        COL_OPTION_IN = COL_OPTION_IN.upper()
    DICT_OUTPUT = pd.DataFrame()
    key_list = list(DICT_INPUT.keys())
    value_list = list(DICT_INPUT.values())
    if DIRECTION in ["H", "HORIZONTAL", "LATERAL", "L"]:
        if COL_OPTION_IN in ["KEYS", "KEY", "K"]:
            for k in key_list:
                v = DICT_INPUT[k]
                if type(v) == dict:
                    v = dict_to_list(v)
                if type(v) == list:
                    count = 0
                    while count < len(v):
                        DICT_OUTPUT.loc[count, k] = v[count]
                        count = count + 1
                else:
                    DICT_OUTPUT.loc[0, k] = DICT_INPUT[k]
        elif COL_OPTION_IN in ["INDEX", "IN", "NUMBER"]:
            DICT_OUTPUT["keys"] = DICT_INPUT.keys()
            DICT_OUTPUT["values"] = DICT_INPUT.values()
            DICT_OUTPUT = transpose_df(DICT_OUTPUT, "v", 1)
    elif DIRECTION in ["V", "VERTICAL"]:
        if isinstance(COL_OPTION_IN, dict) == True:
            for key in COL_OPTION_IN:
                if key.upper() in ["K", "KEYS", "KEY"]:
                    key_col = COL_OPTION_IN[key]
                if key.upper() in ["V", "VALUES", "VALUE"]:
                    value_col = COL_OPTION_IN[key]
        elif (
            isinstance(COL_OPTION_IN, list) == True
            or isinstance(COL_OPTION_IN, tuple) == True
        ):
            key_col = COL_OPTION_IN[0]
            value_col = COL_OPTION_IN[1]
        else:
            key_col = "keys"
            value_col = "values"
        DICT_OUTPUT[key_col] = DICT_INPUT.keys()
        DICT_OUTPUT[value_col] = DICT_INPUT.values()
    return DICT_OUTPUT


def input2list(Input_Argument=None):
    """
    ### Converts an input argument to a list.
    """
    if (
        isinstance(Input_Argument, float)
        or isinstance(Input_Argument, int)
        or isinstance(Input_Argument, str)
    ):
        Output = [Input_Argument]
    elif isinstance(Input_Argument, list):
        Output = Input_Argument[:]
    elif isinstance(Input_Argument, tuple):
        Output = list(Input_Argument)
    elif isinstance(Input_Argument, dict):
        Output = list(Input_Argument.keys())
    else:
        Output = None
    return Output


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


def cross_dict_value_list(DICT_IN={}):
    '''
    ### Cross all the lists (all combinations) in a dictionary (i.e. parameters grid)
    '''
    import itertools
    keys, values = zip(*DICT_IN.items())
    OUTPUT       = [dict(zip(keys, v)) for v in itertools.product(*values)]
    return OUTPUT


def runtime(start_time_point=None, out=False):
    from datetime import datetime
    if start_time_point == None:
        time = datetime.now()
    else:
        time = str(datetime.now() - start_time_point)
        if out == False:
            print (time)
            time = None
        else:
            time = time.split(".")[0]
    return time


def df_col_relocate(DF_IN, COL, LOCATION=0):
    '''
    ### Relocate a column of a dataframe to asigned location
    # DF_IN: input dataframe
    # COL:   a string of column name
    # LOCATION: an int for index location
    '''
    if isinstance(COL, str) == True and isinstance(LOCATION, int) == True:
        col_list = DF_IN.columns.tolist()
        if   LOCATION > len(col_list):
            LOCATION = len(col_list)
        elif LOCATION < 0:
            LOCATION = 0
        col_list.remove(COL)
        col_list = col_list[:LOCATION] + [COL] + col_list[LOCATION:]
        OUTPUT = DF_IN.reindex(columns = col_list)
    else:
        OUTPUT = None
    return OUTPUT


########################################################################################################################
### Model training functions ###
########################################################################################################################


def classifier_and_parameter_grid_select(KEYWORD=''):
    '''
    ### Choose a classifier and its corresponding pre-set parameter grid
    '''
    if isinstance(KEYWORD, str) == True:
        KEYWORD = KEYWORD.upper()
        if   KEYWORD in ['LOGREG', 'LOG_REG', 'LOGISTIC']:
            from sklearn.linear_model import LogisticRegression
            MODEL      = LogisticRegression()
            PARAM_GRID = db.Logistic_Regression_parameter_grid.copy()
            MODEL_NAME = 'LogReg'
        elif KEYWORD in ['RIDGE']:
            from sklearn.linear_model import RidgeClassifier
            MODEL      = RidgeClassifier()
            PARAM_GRID = db.Ridge_parameter_grid.copy()
            MODEL_NAME = 'Ridge'
        elif KEYWORD in ['DT', 'DECISION_TREE', 'TREE']:
            from sklearn.tree import DecisionTreeClassifier
            MODEL      = DecisionTreeClassifier()
            PARAM_GRID = db.Decision_Tree_parameter_grid.copy()
            MODEL_NAME = 'DT'
        elif KEYWORD in ['RF', 'RANDON_FOREST', 'FOREST']:
            from sklearn.ensemble import RandomForestClassifier
            MODEL      = RandomForestClassifier()
            PARAM_GRID = db.Random_Forest_parameter_grid.copy()
            MODEL_NAME = 'RF'
        elif KEYWORD in ['ET', 'EXTRA_TREE']:
            from sklearn.tree import ExtraTreeClassifier
            MODEL      = ExtraTreeClassifier()
            PARAM_GRID = db.Extra_Tree_parameter_grid.copy()
            MODEL_NAME = 'ET'
        elif KEYWORD in ['KNN', 'KNEAREST']:
            from sklearn.neighbors import KNeighborsClassifier
            MODEL      = KNeighborsClassifier()
            PARAM_GRID = db.KNearestNeighbors_parameter_grid.copy()
            MODEL_NAME = 'KNN'
        elif KEYWORD in ['GAUSSIAN', 'GAU']:
            from sklearn.gaussian_process import GaussianProcessClassifier
            MODEL      = GaussianProcessClassifier()
            PARAM_GRID = db.Gaussian_parameter_grid.copy()
            MODEL_NAME = 'Gaussian'
        elif KEYWORD in ['NB', 'Gaussian_NB', 'Naive_Bayes']:
            from sklearn.naive_bayes import GaussianNB
            MODEL      = GaussianNB()
            PARAM_GRID = db.Gaussian_NaiveBayes_parameter_grid.copy()
            MODEL_NAME = 'NB'
        elif KEYWORD in ['MNB', 'Multinomial_NB']:
            from sklearn.naive_bayes import MultinomialNB
            MODEL      = MultinomialNB()
            PARAM_GRID = db.Multinomial_NaiveBayes_parameter_grid.copy()
            MODEL_NAME = 'MNB'
        elif KEYWORD in ['LinearSVC', 'LSVC', 'LinearSVM', 'Linear_SVC', 'Linear_SVM']:
            from sklearn.svm import LinearSVC
            MODEL      = LinearSVC()
            PARAM_GRID = db.LinearSVC_parameter_grid.copy()
            MODEL_NAME = 'LSVC'
        elif KEYWORD in ['SVC', 'SVM']:
            from sklearn.svm import SVC
            MODEL      = SVC()
            PARAM_GRID = db.SVC_parameter_grid.copy()
            MODEL_NAME = 'SVC'
        elif KEYWORD in ['LDA']:
            from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
            MODEL      = LinearDiscriminantAnalysis()
            PARAM_GRID = db.LDA_parameter_grid.copy()
            MODEL_NAME = 'LDA'
        elif KEYWORD in ['QDA']:
            from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
            MODEL      = QuadraticDiscriminantAnalysis()
            PARAM_GRID = db.QDA_parameter_grid.copy()
            MODEL_NAME = 'QDA'
        elif KEYWORD in ['SGD', 'Gradient_Descent']:
            from sklearn.linear_model import SGDClassifier
            MODEL      = SGDClassifier()
            PARAM_GRID = db.Stochastic_Gradient_Descent_parameter_grid.copy()
            MODEL_NAME = 'SGD'
        elif KEYWORD in ['AB', 'ADABOOST']:
            from sklearn.ensemble import AdaBoostClassifier
            MODEL      = AdaBoostClassifier()
            PARAM_GRID = db.AdaBoost_parameter_grid.copy()
            MODEL_NAME = 'AB'
        elif KEYWORD in ['LGBM']:
            from lightgbm import LGBMClassifier
            MODEL      = LGBMClassifier()
            PARAM_GRID = db.LGBM_parameter_grid.copy()
            MODEL_NAME = 'LGBM'
        elif KEYWORD in ['XGB', 'XGBOOST']:
            from xgboost.sklearn import XGBClassifier
            MODEL      = XGBClassifier()
            PARAM_GRID = db.XGB_parameter_grid.copy()
            MODEL_NAME = 'XGB'
        elif KEYWORD in ['MLC', 'NN']:
            from sklearn.neural_network import MLPClassifier
            MODEL      = MLPClassifier()
            PARAM_GRID = db.MLC_parameter_grid.copy()
            MODEL_NAME = 'MLC'
        else:
            MODEL, PARAM_GRID, MODEL_NAME = None
    else:
        MODEL, PARAM_GRID, MODEL_NAME = None
    return MODEL, PARAM_GRID, MODEL_NAME


def cross_val_score(MODEL, X, y, **params):
    '''
    ### Perfrom cross validataion in sklearn
    '''
    from sklearn.model_selection import cross_val_score
    OUTPUT = cross_val_score(MODEL, X, y, **params)
    return list(OUTPUT)


def accuracy_score(model, X, y=None, **params):
    '''
    ### Calculate the accuracy_score of the model and input data
    '''
    from sklearn.metrics import accuracy_score
    try:
        y_pred = model.predict(X)
        y_true = y
    except:
        y_pred = X
        y_true = model
    OUTPUT = accuracy_score(y_true, y_pred, **params)
    return OUTPUT



def balanced_accuracy_score(model, X, y=None, **params):
    '''
    ### Calculate the balanced_accuracy_score of the model and input data
    '''
    from sklearn.metrics import balanced_accuracy_score
    try:
        y_pred = model.predict(X)
        y_true = y
    except:
        y_pred = X
        y_true = model
    OUTPUT = balanced_accuracy_score(y_true, y_pred, **params)
    return OUTPUT


def roc_auc_score(model, X, y): 
    '''
    ### Calculate the Area-Under-Curve of Receiver-Operating-Characteristic curve
    '''
    from sklearn.metrics import roc_auc_score
    y_pred = model.predict_proba(X)
    y_pred = np.transpose([pred[1] for pred in y_pred])
    OUTPUT = roc_auc_score(y, y_pred)
    return OUTPUT
    
    
def precision_score(model, X, y=1, POS_LABEL=1, **params):
    '''
    ### Calculate the precision_score of the model and input data
    '''
    from sklearn.metrics import precision_score
    try:
        y_pred = model.predict(X)
        y_true = y
        LABELS = model.classes_
    except:
        y_pred    = X
        y_true    = model
        POS_LABEL = y
        LABELS    = None
    OUTPUT = precision_score(y_true, y_pred, 
                             labels    = LABELS, 
                             pos_label = POS_LABEL,
                             **params)
    return OUTPUT


def recall_score(model, X, y, POS_LABEL=1, **params):
    '''
    ### Calculate the recall_score of the model and input data
    '''
    from sklearn.metrics import recall_score
    try:
        y_pred = model.predict(X)
        y_true = y
        LABELS = model.classes_
    except:
        y_pred    = X
        y_true    = model
        POS_LABEL = y
        LABELS    = None
    OUTPUT = recall_score(y_true, y_pred, 
                          labels    = LABELS, 
                          pos_label = POS_LABEL,
                          **params)
    return OUTPUT


def f1_score(model, X, y, POS_LABEL=1, **params):
    '''
    ### Calculate the f1_score of the model and input data
    '''
    from sklearn.metrics import f1_score
    try:
        y_pred = model.predict(X)
        y_true = y
        LABELS = model.classes_
    except:
        y_pred    = X
        y_true    = model
        POS_LABEL = y
        LABELS    = None
    OUTPUT = f1_score(y_true, y_pred, 
                      labels    = LABELS, 
                      pos_label = POS_LABEL,
                      **params)
    return OUTPUT


def pr_auc_score(model, X_test, y_test, Pos_Label=1):
    '''
    ### Calculate the Area-Under-Curve of the Precision-Recall curve
    '''
    from sklearn.metrics import precision_recall_curve
    from sklearn.metrics import auc
    y_pred = model.predict_proba(X_test)
    y_true = y_test
    for loc, value in enumerate(model.classes_):
        if value == Pos_Label:
            Pos_Ind = loc
    prec, recall, thresholds = precision_recall_curve(y_true, y_pred[:, Pos_Ind], pos_label=Pos_Label)
    auc_score = auc(recall, prec)
    return auc_score
    
    
def matthews_corrcoef(model, X, y=None, **params):
    '''
    ### Calculate the matthews_corrcoef of the model and input data
    '''
    from sklearn.metrics import matthews_corrcoef
    try:
        y_pred = model.predict(X)
        y_true = y
    except:
        y_pred = X
        y_true = model
    OUTPUT = matthews_corrcoef(y_true, y_pred, **params)
    return OUTPUT
    
    
def ap_score(model, X_test, y_test, Pos_Label=1):
    '''
    ### Calculate the Average Precision score (AP) of the Precision-Recall curve
    '''
    y_pred = model.predict_proba(X_test)
    y_true = y_test
    for loc, value in enumerate(model.classes_):
        if value == Pos_Label:
            Pos_Ind = loc
    from sklearn.metrics import average_precision_score
    aps = average_precision_score(y_true, y_pred[:, Pos_Ind], pos_label=Pos_Label)
    return aps
    

def classifier_true_false_positive_negative(model, X_test, y_test, Pos_Label=None, **params):
    '''
    ### Calculate the the True/False Positive/Negative and their rates in a 2-class classifier model
    '''
    params       = PARAMS_CAP(**params)
    POS_LABEL_IN = PARAMS(['POS_LABEL', 'POS'], Pos_Label, **params)
    OUTPUT_OPT   = PARAMS(['OUTPUT', 'OPTION', 'OPT', 'OUT'], None, **params)
    y_true = list(y_test)
    y_pred = list(model.predict(X_test))
    TP,FP,TN,FN = 0, 0, 0, 0
    class_list = model.classes_.tolist()
    if Pos_Label is None:
        neg_label, pos_label = class_list
    else:
        pos_label = Pos_Label
        class_list.remove(Pos_Label)
        neg_label = class_list[0]
    for i in range(len(y_pred)): 
        if y_true[i]==y_pred[i]==pos_label:
            TP += 1
        elif y_pred[i]==pos_label and y_true[i]!=y_pred[i]:
            FP += 1
        elif y_true[i]==y_pred[i]==neg_label:
            TN += 1
        elif y_pred[i]==neg_label and y_true[i]!=y_pred[i]:
            FN += 1
        else:
            pass
    try:
        TPR = TP / (TP + FN) 
    except ZeroDivisionError:
        TPR = np.nan
    try:
        TNR = TN / (TN + FP) 
    except ZeroDivisionError:
        TNR = np.nan
        
    try:
        FPR = FP / (FP + TN) 
    except ZeroDivisionError:
        FPR = np.nan
    try:
        FNR = FN / (FN + TP) 
    except ZeroDivisionError:
        FNR = np.nan
    try:
        Precision = TP / (TP + FP)
    except ZeroDivisionError:
        Precision = np.nan
    try:
        Recall = TP / (TP + FN)
    except ZeroDivisionError:
        Recall = np.nan
    try:
        Sensitivity = TP / (TP + FN) 
    except ZeroDivisionError:
        Sensitivity = np.nan
    try:
        Specificity = TN / (TN + FP) 
    except ZeroDivisionError:
        Specificity = np.nan
    Missed_Label = TP+FP+TN+FN-len(y_true)
    if OUTPUT_OPT in ['dict', 'DICT', dict]:
        OUTPUT = {
            'TP'           : TP,
            'TN'           : TN,
            'FP'           : FP,
            'FN'           : FN,
            'TPR'          : TPR,
            'TNR'          : TNR,
            'FPR'          : FPR,
            'FNR'          : FNR,
            'Precision'    : Precision,
            'Recall'       : Recall,
            'Sensitivity'  : Sensitivity,
            'Specificity'  : Specificity,
            'PosLabel'     : pos_label,
            'Missed_Label' : Missed_Label}
    else:
        OUTPUT = None
        print (f'          TP : {TP}')
        print (f'          TN : {TN}')
        print (f'          FP : {FP}')
        print (f'          FN : {FN}')
        print (f'         TPR : {TPR}')
        print (f'         TNR : {TNR}')
        print (f'         FPR : {FPR}')
        print (f'         FNR : {FNR}')
        print (f'   Precision : {Precision}')
        print (f'      Recall : {Recall}')
        print (f' Sensitivity : {Sensitivity}')
        print (f' Specificity : {Specificity}')
        print (f'    PosLabel : {pos_label}')
        print (f'Missed_Label : {Missed_Label}')
    return OUTPUT


def classifier_test_scores(model, X_test, y_test, POS_LABEL=None, **params):
    '''
    ### Calculate mulitple scores for evaluation of classifier models
    '''
    params       = PARAMS_CAP(**params)
    POS_LABEL_IN = PARAMS(['POS_LABEL', 'POS'], POS_LABEL, **params)
    OUTPUT_OPT   = PARAMS(['OUTPUT', 'OPTION', 'OPT', 'OUT'], None, **params)
    AS  = accuracy_score(model, X_test, y_test)
    BAS = balanced_accuracy_score(model, X_test, y_test)
    try:
        RAS = roc_auc_score(model, X_test, y_test)
    except:
        RAS = None
    if POS_LABEL_IN is None:
        PS, RS, F1S, AP, PAS, MCC = None, None, None, None, None, None
    else:
        PS  = precision_score(model, X_test, y_test, POS_LABEL)
        RS  = recall_score(model, X_test, y_test, POS_LABEL)
        F1S = f1_score(model, X_test, y_test, POS_LABEL)
        AP  = ap_score(model, X_test, y_test, POS_LABEL)
        PAS = pr_auc_score(model, X_test, y_test, POS_LABEL)
        MCC = matthews_corrcoef(model, X_test, y_test)
    if OUTPUT_OPT in ['dict', 'DICT', dict]:
        OUTPUT = {'Accuracy'          : AS,
                  'Balanced_Accuracy' : BAS,
                  'Roc_Auc_Score'     : RAS,
                  'Average_Precision' : AP,
                  'PR_Auc_Score'      : PAS,
                  'Positive_Label'    : POS_LABEL_IN,
                  'Precision_Score'   : PS,
                  'Recall_Score'      : RS,
                  'F1_score'          : F1S,
                  'Matthews_CorrCoef' : MCC}
    else:
        OUTPUT = None
        print (f'         Accuracy : {AS}')
        print (f'Balanced_Accuracy : {BAS}')
        print (f'    Roc_Auc_Score : {RAS}')
        print (f'Average_Precision : {AP}')
        print (f'     PR_Auc_Score : {PAS}')
        print (f'   Positive_Label : {POS_LABEL_IN}')
        print (f'  Precision_Score : {PS}')
        print (f'     Recall_Score : {RS}')
        print (f'         F1_score : {F1S}')
        print (f'Matthews_CorrCoef : {MCC}')
    return OUTPUT
    

########################################################################################################################
### Other functions ###
########################################################################################################################


def verify_id_col(*arg):
    """
    ### Verify all the dataframe column 'id' are identical
    """
    df0 = arg[0]
    for df in arg:
        if df0["id"].tolist() != df["id"].tolist():
            raise ValueError("Virus ID order is wrong!")
            break
    print("Virus ID order is correct.")


def rename_dataset_columns(Dataset_In, Note=""):
    """
    ### Add a noting of dataset at the head of every columns of dataset
    """
    Dataset_Out = Dataset_In.copy()
    col_list = []
    for col in Dataset_In.columns.tolist():
        col_list.append(f"[{Note}]_{col}")
    Dataset_Out.columns = col_list
    Dataset_Out = Dataset_Out.reset_index(drop=True)
    return Dataset_Out


def df_col_numerise(SERIES_IN, REF_DICT=None):
    '''
    ### Numerise (int) a series or a column of dataframe according to the REF_DICT
    '''
    OUTPUT = {}
    list_ = SERIES_IN.tolist()
    if REF_DICT is None:
        REF_DICT = {}
        for i, j in enumerate(list(set(set(list_)))):
            REF_DICT[j] = i
    for loc, element in enumerate(list_):
        OUTPUT[loc] = REF_DICT[element]
    OUTPUT = pd.Series(OUTPUT)
    return OUTPUT


def load_RTG_Data(Keyword=None):
    """
    ### Load RTG data with keyword
    """
    if isinstance(Keyword, str):
        keyword = Keyword.upper()
    else:
        keyword = None
    if keyword is None:
        Output_Data = None
    elif keyword in ["CU"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//CU.csv")
    elif keyword in ["AAU"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//AAU.csv")
    elif keyword in ["RSCU"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//RSCU.csv")
    elif keyword in ["INFO"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//INFO.csv")
    elif keyword in ["HOST"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//HOST.csv")
    elif keyword in ["PARTITE", "P", "PART"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//PARTITE.csv")
    elif keyword in [
        "GENOME_LENGTH",
        "GEN_LEN",
        "LENGTH",
        "LEN",
        "CDS_LENGTH",
        "CDS_LEN",
    ]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//CDS_LENGTH.csv")
    elif keyword in ["ATGC"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//ATGC.csv")
    elif keyword in ["START_STOP", "SS", "STST"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//START_STOP.csv")
    elif keyword in ["HS_CORR", "CORR"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//HS_CORR.csv")
    elif keyword in ["HS_CORR_AA", "CORR_AA", "CORRAA"]:
        Output_Data = pd.read_csv("RTG_Data//ALL_VIRUS//HS_CORR_AA.csv")
    else:
        Output_Data = None
    return Output_Data
