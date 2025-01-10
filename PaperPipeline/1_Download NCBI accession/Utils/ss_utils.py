### Import basic modules

import numpy as np
import pandas as pd
from scipy import stats
import time
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord

from Utils import ss_database as db

slash_mark = "//"
input_folder = ""
output_folder = ""

# API key for NCBI - Input your NCBI API
NCBI_API_KEY = "8ad0809ddce384eb810757af76eb9abdce09"


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
                + time_check
                + "_"
                + author_name
            )
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


def discriptive_data(in0, raw_output=False):
    """
    ### Calculate multiple descriptive data
    """
    import scipy

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
            print("not supported input.")
            pass
    if len(list_data) == 0:
        dict_describe = {
            "count": None,
            "sum": None,
            "mean": None,
            "median": None,
            "variance": None,
            "std": None,
            "min": None,
            "max": None,
            "25%": None,
            "50%": None,
            "75%": None,
            "skewness": None,
            "kurtosis": None,
            "weighted_average": None,
            "gmean_0": None,
        }
    else:
        dict_describe["count"] = len(list_data)
        dict_describe["sum"] = sum(list_data)
        dict_describe["mean"] = np.mean(list_data)
        dict_describe["median"] = np.median(list_data)
        dict_describe["variance"] = np.var(list_data)
        dict_describe["std"] = np.std(list_data)
        dict_describe["min"] = min(list_data)
        dict_describe["max"] = max(list_data)
        dict_describe["25%"] = np.quantile(list_data, 0.25)
        dict_describe["50%"] = np.quantile(list_data, 0.5)
        dict_describe["75%"] = np.quantile(list_data, 0.75)
        dict_describe["skewness"] = scipy.stats.skew(list_data)
        dict_describe["kurtosis"] = scipy.stats.kurtosis(list_data)
        dict_describe["weighted_average"] = np.average(list_data)
        if list_data.count(0) > 0:
            dict_describe["gmean_0%"] = list_data.count(0) / len(list_data)
            list_data_gmean = list_data[:]
            while 0 in list_data_gmean:
                list_data_gmean.remove(0)
            gmean = stats.gmean(list_data_gmean)
        else:
            gmean = stats.gmean(list_data)
            dict_describe["gmean_0"] = 0
        dict_describe["gmean"] = gmean
    if raw_output is True:
        dict_describe["raw_data"] = list_data
    return dict_describe


def one_hot_encode_with_unknown(
    labels, true_labels_in=True, false_labels_in=False, unknown_labels_in=""
):
    """
    ### One-hot encodes a list of labels with true, false, and unknown encodings.
    """
    true_labels = input2list(true_labels_in)
    false_labels = input2list(false_labels_in)
    unknown_labels = input2list(unknown_labels_in)

    unique_labels = list(set(labels))
    if unknown_labels in unique_labels:
        unknown_index = unique_labels.index(unknown_labels)
        unique_labels.pop(unknown_index)
        unique_labels.append(unknown_labels)
    else:
        unknown_index = -1
        unique_labels.append(unknown_labels)

    encoded_labels = []
    for label in labels:
        if label in true_labels:
            encoded_label = 1
        elif label in false_labels:
            encoded_label = -1
        elif label in unknown_labels:
            encoded_label = 0
        else:
            raise ValueError(
                "Label '{}' not found in true, false, or unknown labels".format(label)
            )
        encoded_labels.append(encoded_label)

    return encoded_labels


########################################################################################################################
### Biology related functions ###
########################################################################################################################


def get_all_cds_detail(virus_genome_id="", option="", **params):
    """
    ### Obtain a combined CDS from NCBI according to accession
    # INPUT: single NCBI ID
    # OUTPUT: String of combined all CDS of virus genome
    """
    params = PARAMS_CAP(**params)
    WARNING = PARAMS(["WARN", "WARNING", "W"], True, **params)
    dict_output = {}
    try:
        hd1 = Entrez.efetch(
            db="nucleotide", id=virus_genome_id, rettype="gb", api_key=NCBI_API_KEY
        )
        seq = SeqIO.read(hd1, "gb")
        features = seq.features
        for feature in features:
            if feature.type == "CDS":
                cds_info_dict = {
                    "type": "CDS",
                    "location": str(feature.location),
                    "location_extract": feature.location,
                }
                for i in feature.qualifiers:
                    if len(feature.qualifiers[i]) == 1:
                        cds_info_dict[i] = str(feature.qualifiers[i][0])
                    else:
                        cds_info_dict[i] = feature.qualifiers[i]
                cds = SeqRecord(feature.location.extract(seq).seq)
                cds_info_dict["mRNA_sequence"] = str(cds.seq)
                try:
                    dict_output[cds_info_dict["gene"]] = cds_info_dict
                except:
                    dict_output[
                        "protein:" + cds_info_dict["protein_id"]
                    ] = cds_info_dict
    except:
        if WARNING == True:
            print("-Error occured on " + virus_genome_id)
            print("--CANNOT download sequence from internet.")
        pass
    if option == "df":
        df0 = empty_df()
        c = 0
        for i in dict_output:
            df0.loc[c, "ID"] = i
            for ii in dict_output[i]:
                if type(dict_output[i][ii]) == str:
                    df0.loc[c, ii] = dict_output[i][ii]
                else:
                    df0.loc[c, ii] = str(dict_output[i][ii])
            c += 1
        dict_output = df0
    return dict_output


def get_combined_cds(virus_genome_id):
    """
    ### Combine all the cds into catencated sequence
    """
    dict_output = {}
    dict_input = get_all_cds_detail(virus_genome_id)
    combined_CDS = ""
    for i in dict_input:
        combined_CDS = combined_CDS + dict_input[i]["mRNA_sequence"]
    dict_output[virus_genome_id] = combined_CDS
    return dict_output


def segmented_get_combined_cds(virus_genome_id_dict={}):
    """
    ### Combine all the cds of a dictionary of segmented virus into catencated sequence
    """
    dict_output = {}
    for i in virus_genome_id_dict:
        for ii in virus_genome_id_dict[i]:
            combined_cds = list(get_combined_cds(ii).values())[0]
            if i not in dict_output.keys():
                dict_output[i] = combined_cds
            else:
                dict_output[i] = dict_output[i] + combined_cds
    return dict_output


def ATGC_cal(seq=""):
    """
    ### GC% and total codons number calculation
    # INPUT: mRNA sequence
    # OUTPUT: List contains Codon number, GC%, Error nucleotide numer
    """
    ATGC_dict = {}
    GS = ""
    if isinstance(seq, str):
        GS = seq
    elif isinstance(seq, list):
        for i in seq:
            GS = GS + i
    GSL = len(GS)
    if GSL == 0:
        ATGC_dict = {
            "bp#": 0,
            "/3?": np.nan,
            "A#": np.nan,
            "A%": np.nan,
            "T#": np.nan,
            "T%": np.nan,
            "G#": np.nan,
            "G%": np.nan,
            "C#": np.nan,
            "C%": np.nan,
            "AT#": np.nan,
            "AT%": np.nan,
            "GC#": np.nan,
            "GC%": np.nan,
            "not_ATGC#": np.nan,
        }
    else:
        ATGC_dict["bp#"] = GSL
        GSLtest = GSL / 3
        if GSLtest - int(GSLtest) == 0:
            ATGC_dict["/3?"] = "Y"
        else:
            ATGC_dict["/3?"] = "N"
        ATGC_dict["A#"] = GS.count("A")
        ATGC_dict["A%"] = GS.count("A") / GSL
        ATGC_dict["T#"] = GS.count("T")
        ATGC_dict["T%"] = GS.count("T") / GSL
        ATGC_dict["G#"] = GS.count("G")
        ATGC_dict["G%"] = GS.count("G") / GSL
        ATGC_dict["C#"] = GS.count("C")
        ATGC_dict["C%"] = GS.count("C") / GSL
        ATGC_dict["AT#"] = ATGC_dict["A#"] + ATGC_dict["T#"]
        ATGC_dict["AT%"] = ATGC_dict["A%"] + ATGC_dict["T%"]
        ATGC_dict["GC#"] = ATGC_dict["G#"] + ATGC_dict["C#"]
        ATGC_dict["GC%"] = ATGC_dict["G%"] + ATGC_dict["C%"]
        ATGC_dict["not_ATGC#"] = GSL * (1 - (ATGC_dict["AT%"] + ATGC_dict["GC%"]))
    return ATGC_dict


def codon_dict(GS="", choice=""):
    """
    ### Convert a coding sequence (CDS) into a dictionary with items (location : codon)
    """
    codon_dict = {}
    if choice == "pair":
        optional_length = 3
    elif type(choice) == int:
        optional_length = choice
    else:
        optional_length = 0
    GSL = len(GS)
    GSLtest0 = GSL / 3
    GSLtest = GSL / 3
    position = 1
    if GSLtest - int(GSLtest) == 0:
        while GSLtest > 0:
            codon = GS[(position - 1) * 3 : position * 3 + optional_length]
            codon_dict[position] = codon
            GSLtest = GSLtest - 1
            position = position + 1
    else:
        print("Cannot divide by 3.")
    if choice == "pair":
        del codon_dict[len(codon_dict)]
    return codon_dict


def codon_list(GS=""):
    """
    ### Convert a coding sequence (CDS) into a list of codons
    """
    codon_d = codon_dict(GS)
    codon_list = list(codon_d.values())
    return codon_list


def codon_usage(seq=""):
    """
    ### Calculate codon usage (or codon percentages) from a coding sequence
    """
    codon_table_list = db.codon_table_list[:]
    if len(seq) == 0:
        codon_usage_dict = db.empty_codon_usage.copy()
    else:
        codon_usage_dict = {}
        codon_l = codon_list(seq)
        total_codon = len(seq) / 3
        for i in codon_table_list:
            codon_usage_dict[i] = codon_l.count(i) / total_codon
    return codon_usage_dict


def aa_usage(seq=""):
    """
    ### Calculate amino acid usage (or amino acid percentages) from a coding sequence
    """
    aa_table_list = db.aa_table_list[:]
    translation_aa_codon_dict = db.translation_aa_codon_dict.copy()
    aa_usage_dict = {}
    aa = 0
    codon_usage_dict = codon_usage(seq)
    for i in aa_table_list:
        for ii in translation_aa_codon_dict[i]:
            aa = aa + codon_usage_dict[ii]
        aa_usage_dict[i] = aa
        aa = 0
    return aa_usage_dict


def codon_usage_to_aa_usage(codon_usage={}):
    """
    ### Generate a dictionary of AA usage (Amino acid percentages) from a dictionary of CU (codon usage)
    """
    dict_output = {}
    codon_list = []
    aa_table_list = db.aa_table_list[:]
    translation_aa_codon_dict = db.translation_aa_codon_dict.copy()
    for i in aa_table_list:
        codon_list = translation_aa_codon_dict[i]
        codon_p = 0
        for ii in codon_list:
            codon_p = codon_p + codon_usage[ii]
        dict_output[i] = codon_p
    return dict_output


def codon_usage_to_RSCU_usage(codon_usage_dict={}):
    """
    ### Generate RSCU dictionary from a dictionary of CU (codon usage)
    """
    dict_output = {}
    translation_codon_aa_dict = db.translation_codon_aa_dict.copy()
    translation_aa_codon_dict = db.translation_aa_codon_dict.copy()

    aa_usage_dict = codon_usage_to_aa_usage(codon_usage_dict)
    for i in codon_usage_dict:
        aa = translation_codon_aa_dict[i]
        if codon_usage_dict[i] == 0.0:
            rscu = 0.0
        else:
            rscu = (codon_usage_dict[i] / aa_usage_dict[aa]) * len(
                translation_aa_codon_dict[aa]
            )
        dict_output[i] = rscu
    return dict_output


########################################################################################################################
### Taxonomy related functions ###
########################################################################################################################


def taxonomy_search(keyword=""):
    """
    ### Search NCBI taxonomy database with keyword to find taxonomy_id
    """
    handle = Entrez.esearch(db="taxonomy", term=keyword, api_key=NCBI_API_KEY)
    record = Entrez.read(handle)
    id_list = record["IdList"]
    return id_list


def taxonomy_lineage_download(tax_id=''):
    '''
    ### Download and organise lineage information with taxonomy_id
    '''
    handle = Entrez.efetch(db='taxonomy', id=tax_id, api_key=NCBI_API_KEY)
    record = handle.read()
    run_list = record.decode().strip().split("\n")
    current_section = None
    current_sub_section = None
    dict_taxon = {}
    for line in run_list:
        if '<LineageEx>' in (line[0:16]):
            current_section = 'Lineage'
        elif '</LineageEx>' in (line[0:16]):
            current_section = None
        if current_section == None:
            if '<ScientificName>' in line:
                virus_name = del_str(line,['    <ScientificName>','</ScientificName>'])
        elif current_section == 'Lineage':
            if '<Taxon>' in line[0:16]:
                current_sub_section = 'Taxon'
            elif '</Taxon>' in line[0:16]:
                current_sub_section = None
            if current_sub_section == 'Taxon':
                if '<ScientificName>' in line:
                    value = del_str(line,['            <ScientificName>','</ScientificName>'])
                elif '<Rank>' in line:
                    key = del_str(line,['            <Rank>','</Rank>'])
                    while key in dict_taxon:
                        key += '_'
                    dict_taxon[key] = value
    return dict_taxon

    
def taxonomy_encoder(DF_In, Unknown_Label="Unclassified"):
    """
    ### One-hot encodes a pandas DataFrame where each column represents a categorical variable with a hierarchical taxonomy.
        True = 1, False = -1, Unknown = 0
    """
    if not isinstance(DF_In, pd.DataFrame):
        raise ValueError("Input must be pd.DataFrame")

    DF_Output = pd.DataFrame()

    for col in DF_In.columns.tolist():

        data_list = DF_In[col].tolist()
        class_list = list(DF_In[col].unique())
        if Unknown_Label in class_list:
            class_list.remove(Unknown_Label)

        for cla in class_list:

            negative_cla = class_list[:]
            negative_cla.remove(cla)

            encoded = one_hot_encode_with_unknown(
                data_list, cla, negative_cla, Unknown_Label
            )

            DF_Output[f"{col}_{cla}"] = encoded
    return DF_Output


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
