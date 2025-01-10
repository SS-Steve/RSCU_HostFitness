### System requirement ###

import numpy as np
import pandas as pd
from scipy import stats
import sklearn
import time
import math
import os
import multiprocessing
import tqdm
from tqdm.notebook import tqdm as tb

import Utils.ss_utils as ss
import Utils.ss_database as db



### -------------------------------------------------------------------------------- ###
### Define X
### -------------------------------------------------------------------------------- ###

# Load Ready-To-Go datasets
RSCU = ss.load_RTG_Data("RSCU")
CU = ss.load_RTG_Data("CU")
AAU = ss.load_RTG_Data("AAU")
CDS_LENGTH = ss.load_RTG_Data("CDS_LENGTH")
ATGC = ss.load_RTG_Data("ATGC")
HOST = ss.load_RTG_Data("HOST")
INFO = ss.load_RTG_Data("INFO")
PARTITE = ss.load_RTG_Data("PARTITE")
TAXONOMY = ss.load_RTG_Data("TAXONOMY")
START_STOP = ss.load_RTG_Data("START_STOP")
HS_CORR = ss.load_RTG_Data("HS_CORR")
HS_CORR_AA = ss.load_RTG_Data("HS_CORR_AA")

# Generate training data X with selected datasets
Dataset_Note = "RSC, CDS_LENGTH, TAXONOMY"
print (f'Dataset: {Dataset_Note}')

X = pd.concat([RSCU, CDS_LENGTH, TAXONOMY], axis=1)



### -------------------------------------------------------------------------------- ###
### Define y
### -------------------------------------------------------------------------------- ###

Host_Selection_List = ["human", "vertebrates", "land plants", "invertebrates", "bacteria"]

POS_LABEL = "human"

y = HOST[POS_LABEL]


### -------------------------------------------------------------------------------- ###
### Define model training parameters
### -------------------------------------------------------------------------------- ###

CPU_USED = os.cpu_count() - 3

CV = 10
TEST_SIZE = 0.05
RANDOM_STATE = 8

NOTING = f"[{POS_LABEL}]"

MODEL_CHOOSE = "RF"
UPSAMPLE_METHOD = "SMOTE"


OUTPUT_ADD_NOTE = (
    f"{Dataset_Note}_{TEST_SIZE}_RS{RANDOM_STATE}_{UPSAMPLE_METHOD}_{NOTING}"
)


# ==================================================== #
# Determine use of self-define parameter ranges (False if want to run default grid)
MANNUAL_INPUT_PARAM_GRID = True
# ==================================================== #
PARAM_GRID = {
    #"n_estimators": list(range(30, 37, 1)),
    "criterion": ["gini"],  # 'gini'
    'splitter': ['random'],
    "max_depth": [None], # , # [None] list(range(40, 501, 10))+ 
    "min_samples_split": list(range(2, 51, 1)),
    "min_samples_leaf": list(range(1, 22, 1)),
    "max_features": ['sqrt'],
    "class_weight": ["balanced"],
}
# ==================================================== #




X_raw = X.copy()
y_raw = y.copy()
X = np.array(X)
y = np.array(y)


### Choose the model and pre-set parameter gird ###
if MANNUAL_INPUT_PARAM_GRID is True:
    MODEL, _, MODELNAME = ss.classifier_and_parameter_grid_select(MODEL_CHOOSE)
else:
    MODEL, PARAM_GRID, MODELNAME = ss.classifier_and_parameter_grid_select(MODEL_CHOOSE)

### Generate parameters grid ###
run_params_list = ss.cross_dict_value_list(PARAM_GRID)


### Train-test split ###
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=TEST_SIZE, shuffle=True, random_state=RANDOM_STATE
)
print("Finish Train-test split.")


if UPSAMPLE_METHOD in ["SMOTE"]:
    ### SMOTE upsampling ###
    from imblearn.over_sampling import SMOTE
    sm = SMOTE(random_state=RANDOM_STATE)
    X_train, y_train = sm.fit_resample(X_train, y_train)
    
elif UPSAMPLE_METHOD in ["COPY"]:
    ### Copy-upsampling ###
    X_train = pd.DataFrame(X_train)
    y_train = pd.Series(y_train)
    X_train, y_train = ss.balance_binary_sample_by_copying(X_train, y_train)
#elif UPSAMPLE_METHOD in ["COPY"]:
    
    
    
print("Finish resampling.")

print("Optimising starts.")
print(f"Total iterations = {len(run_params_list)}")

### Define a process for multiprocessing ###
def PROCESS(params_dict):
    OUTPUT = params_dict.copy()
    model = sklearn.base.clone(MODEL)
    YTRAIN_COUNT = {f"train_{k}": v for k, v in ss.list_element_count(y_train).items()}
    YTEST_COUNT = {f"test_{k}": v for k, v in ss.list_element_count(y_test).items()}
    ### Main program of model fitting and evaluation ###
    try:
        
        # Add parameters
        model.set_params(**params_dict)
        if "random_state" in model.get_params():
            model.set_params(random_state=RANDOM_STATE)
        if "n_jobs" in model.get_params():
            model.set_params(n_jobs=None)
        
        # Train models
        train_rt = ss.runtime()
        model.fit(X_train, y_train)
        TRAIN_RT = ss.runtime(train_rt, True)
        
        # Calculate different scores
        cv_score_list = list(
            ss.cross_val_score(
                model, X_train, y_train, cv=CV, scoring="balanced_accuracy"
            )
        )
        CV_MEAN = np.mean(cv_score_list)
        CV_STD = np.std(cv_score_list)
        for cvi, single_cv in enumerate(cv_score_list):
            OUTPUT["cv_" + str(cvi + 1)] = single_cv

        AS = ss.accuracy_score(model, X_test, y_test)
        BAS = ss.balanced_accuracy_score(model, X_test, y_test)
        PS = ss.precision_score(model, X_test, y_test, POS_LABEL)
        RS = ss.recall_score(model, X_test, y_test, POS_LABEL)
        F1S = ss.f1_score(model, X_test, y_test, POS_LABEL)

        predict_rt = ss.runtime()
        model.predict(X_raw)
        PREDICT_RT = ss.runtime(predict_rt, True)

        NOTE = np.nan

    except Exception as e:

        AS, BAS, CV_MEAN, CV_STD, PS, RS, F1S, TRAIN_RT, PREDICT_RT = [np.nan] * 9
        NOTE = str(e)

    # Organise into output dictionary
    OUTPUT["note"] = NOTE
    OUTPUT.update(YTRAIN_COUNT)
    OUTPUT.update(YTEST_COUNT)
    OUTPUT["train_time"] = TRAIN_RT
    OUTPUT["predict_time"] = PREDICT_RT
    OUTPUT["precision_score"] = PS
    OUTPUT["recall_score"] = RS
    OUTPUT["f1_score"] = F1S
    OUTPUT["cv_mean"] = CV_MEAN
    OUTPUT["cv_std"] = CV_STD
    OUTPUT["accuracy"] = AS
    OUTPUT["balanced_accuracy"] = BAS

    return OUTPUT


### Run multiprocessing ###
POOL = multiprocessing.Pool(processes=CPU_USED)
c = 0
FINAL_RESULT = {}

for _RESULT_ in tb(
    POOL.imap_unordered(PROCESS, run_params_list), total=len(run_params_list)
):
    FINAL_RESULT[c] = _RESULT_
    c += 1
    
    # Save data every 50 runs
    if c % 50 == 0:

        print(f"Process: {c} / {len(run_params_list)}")

        ### Turn into dataframe
        DF_OUTPUT = pd.DataFrame.from_dict(FINAL_RESULT, "index")
        DF_OUTPUT["model_index"] = DF_OUTPUT.index
        DF_OUTPUT = ss.df_col_relocate(DF_OUTPUT, "model_index", 0)
        DF_OUTPUT = DF_OUTPUT.sort_values(
            "f1_score", ascending=False, ignore_index=True
        )

        ### Saving the data
        FILENAME_OUTPUT = f"[{MODELNAME}]_{OUTPUT_ADD_NOTE}"
        ss.create_folder(FILENAME_OUTPUT, "Output")
        ss.database_save(FINAL_RESULT, FILENAME_OUTPUT + "//FINAL_RESULT", "Output")
        ss.save_file(DF_OUTPUT, FILENAME_OUTPUT + "//DF_OUTPUT")
        ss.dict_save_as_txt(PARAM_GRID, FILENAME_OUTPUT + "//PARAM_GRID")

        noting_string = f"Done: {c} / {len(run_params_list)}"
        ss.save_file(noting_string, FILENAME_OUTPUT + "//Note", "txt")


### Turn into dataframe
DF_OUTPUT = pd.DataFrame.from_dict(FINAL_RESULT, "index")
DF_OUTPUT["model_index"] = DF_OUTPUT.index
DF_OUTPUT = ss.df_col_relocate(DF_OUTPUT, "model_index", 0)
DF_OUTPUT = DF_OUTPUT.sort_values("f1_score", ascending=False, ignore_index=True)

### Saving the data
FILENAME_OUTPUT = f"[{MODELNAME}]_{OUTPUT_ADD_NOTE}"
ss.create_folder(FILENAME_OUTPUT, "Output")
ss.database_save(FINAL_RESULT, FILENAME_OUTPUT + "//FINAL_RESULT", "Output")
ss.save_file(DF_OUTPUT, FILENAME_OUTPUT + "//DF_OUTPUT")
ss.dict_save_as_txt(PARAM_GRID, FILENAME_OUTPUT + "//PARAM_GRID")

noting_string = f"Done: {c} / {len(run_params_list)}"
ss.save_file(noting_string, FILENAME_OUTPUT + "//Note", "txt")
