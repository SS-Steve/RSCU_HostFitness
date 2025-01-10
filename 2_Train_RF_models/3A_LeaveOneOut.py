
### This is for Figure 3A.

### System requirement ###

import numpy as np
import pandas as pd
from scipy import stats
import sklearn
import os
from Utils import ss_utils as ss 
from Utils import ss_database as db
 

import optuna
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate


# -------------------------------------------------------------------------------- #
# Define X

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

X_Dict_Feature = {
    #'RSCU':RSCU,
    'CU':CU,
    'AAU':AAU,
    'ATGC':ATGC,
    'CDS_LENGTH':CDS_LENGTH,
    'PARTITE':PARTITE,
    'TAXONOMY':TAXONOMY,
    'START_STOP':START_STOP,
    'HS_CORR':HS_CORR,
    'HS_CORR_AA':HS_CORR_AA,
    }

# Choose dataset here
# Dataset_Note = "RSCU"
# X = pd.concat([RSCU], axis=1)

# Dataset_Note = "RSCUpTAXONOMY"
# X = pd.concat([RSCU, TAXONOMY], axis=1)

Dataset_Note = "RSCUpTAXONOMYpCDS_LENGTH"
X = pd.concat([RSCU, TAXONOMY, CDS_LENGTH], axis=1)

print (f'Dataset: {Dataset_Note}')
    
# -------------------------------------------------------------------------------- #
# Define y
animal_list = ["human", "vertebrates", "land plants", "invertebrates", "bacteria"]

POS_LABEL = "human"
y = HOST[POS_LABEL]
# -------------------------------------------------------------------------------- #

Run_List_Index_Comb = (0, 1000)  ### For running a part of the samples

Study_Title = f"{Dataset_Note}_LOO_[{Run_List_Index_Comb}]"

CPU_USED = os.cpu_count() - 3
#CPU_USED = int(os.cpu_count() / 2) - 5
#CPU_USED     = 5

CV = 20
TEST_SIZE = 'LOO'
RANDOM_STATE = 8
Optuna_Trial_N = 5

MODELNAME = "RF"
UPSAMPLE_METHOD = "SMOTE"

Optuna_CV_Score = 'recall' # Choose 'balanced_accuracy' or 'recall'
CV_Scoring_List = ['average_precision', 'roc_auc', 'precision', 'recall', 'f1', 
                   'accuracy', 'balanced_accuracy']

# Custom output note on output folder
score_savename_ref = {'balanced_accuracy':'BA',
                      'f1':'F1', 
                      'recall':'Recall'}
NOTING = f"[{POS_LABEL}]_for{score_savename_ref[Optuna_CV_Score]}"
OUTPUT_ADD_NOTE = (
    f"{Study_Title}_{TEST_SIZE}_RS{RANDOM_STATE}_{UPSAMPLE_METHOD}_{NOTING}"
)


# Numerise y to 0 and 1
POS_LABEL_Raw = POS_LABEL
NUMERIC_REF = {POS_LABEL: 1, f"not_{POS_LABEL}": 0}
y = ss.df_col_numerise(y, NUMERIC_REF)
POS_LABEL = NUMERIC_REF[POS_LABEL]

X_raw = X.copy()
X = np.array(X)
y_raw = y.copy()
y = np.array(y)

feautre_list = X_raw.columns.tolist()

# Create LOO split list
from sklearn.model_selection import LeaveOneOut
loo = LeaveOneOut()
train_test_split_ind_combination_list = []
for train_ind, test_ind in loo.split(X):
    train_test_split_ind_combination_list.append((train_ind, test_ind))
print (f'{Dataset_Note} - Done LOO split.')


# Create Output folder
FILENAME_OUTPUT = f"[{MODELNAME}]_{OUTPUT_ADD_NOTE}"
ss.create_folder(FILENAME_OUTPUT, "Output")
ss.create_folder(f"{FILENAME_OUTPUT}//study", "Output")

print (f"{Dataset_Note} - START.")
FINAL_RESULT = {}
Run_List_Total = list(range(len(train_test_split_ind_combination_list)))

Run_List = Run_List_Total[Run_List_Index_Comb[0]: Run_List_Index_Comb[1]] # Get a part of the total list


for train_test_split_iter_n in Run_List:
    
    # Get train, test index
    train_ind, test_ind = train_test_split_ind_combination_list[train_test_split_iter_n]
    
    # Train-test split
    X_train = X[train_ind]
    y_train = y[train_ind]
    X_test  = X[test_ind]
    y_test  = y[test_ind]
    
    n_raw_Train = len(y_train)
    n_raw_Test  = len(y_test)
    
    ### SMOTE upsampling ###
    from imblearn.over_sampling import SMOTE
    sm = SMOTE(random_state=RANDOM_STATE)
    X_train, y_train = sm.fit_resample(X_train, y_train)
    n_resampled_Train = len(y_train)
    n_resampled_Test  = len(y_test)
        
    print(f"{Dataset_Note} - Finish resampling.")
    print(f"{Dataset_Note} - Optimising starts.")

    def objective(trial):

        # Set hyperparameter search space
        n_estimators      = trial.suggest_int('n_estimators', 2, 300)
        criterion         = trial.suggest_categorical('criterion', ['entropy', 'gini'])
        min_samples_split = trial.suggest_int('min_samples_split', 2, 20)
        min_samples_leaf  = trial.suggest_int('min_samples_leaf', 1, 10)
        max_features      = trial.suggest_categorical('max_features', ['sqrt', 'log2'])
        class_weight      = trial.suggest_categorical('class_weight', ['balanced'])

        # Initialize and train the model with hyperparameters suggested by Optuna
        clf = RandomForestClassifier(
            n_estimators=n_estimators, 
            criterion=criterion,
            min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            max_features=max_features,
            class_weight=class_weight,
        )
        
        scores = cross_val_score(clf, X_train, y_train, cv=CV, scoring=Optuna_CV_Score, n_jobs=CPU_USED)

        return scores.mean()

    study = optuna.study.create_study(study_name=Dataset_Note, direction='maximize')
    study.optimize(objective, n_trials=Optuna_Trial_N, n_jobs=CPU_USED)


    print (f"{Dataset_Note} - Done hyper-parameters tuning.")

    # Create output dictionary
    OUTPUT = {'Dataset':Dataset_Note,
              'Sample_Index':train_test_split_iter_n,
              'POS_LABEL':POS_LABEL_Raw,
              'CV':CV,
              'Test_Ratio':TEST_SIZE,
              'Random_State':RANDOM_STATE,
              'Optuna_Trial_N':Optuna_Trial_N,
              'Model':MODELNAME,
              'Resample_Nethod':UPSAMPLE_METHOD,
              'Feature_Number':len(feautre_list),
              'Feature_List':feautre_list}


    # Get a dictionary of best hyper-parameters
    params_dict = study.best_params

    OUTPUT.update(params_dict)
    YTRAIN_COUNT = {f"train_{k}": v for k, v in ss.list_element_count(y_train).items()}
    YTEST_COUNT = {f"test_{k}": v for k, v in ss.list_element_count(y_test).items()}


    ### Main program of model fitting and evaluation ###

    # Set parameters
    model = RandomForestClassifier()
    model.set_params(**params_dict)
    if "random_state" in model.get_params():
        model.set_params(random_state=RANDOM_STATE)
    if "n_jobs" in model.get_params():
        model.set_params(n_jobs=CPU_USED)

    # Cross-validations
    cv_res = cross_validate(model, X_train, y_train, scoring=CV_Scoring_List, cv=CV, n_jobs=CPU_USED)

    # Organise CV results into Output
    d_t = {}
    for scoring in CV_Scoring_List:
        data_list = cv_res[f"test_{scoring}"]
        d_t[f"[cv]_{scoring}_mean"] = data_list.mean()
        d_t[f"[cv]_{scoring}_std"] = data_list.std()
        for ind, element in enumerate(data_list):
            OUTPUT[f"[cv]_{scoring}_{ind}"] = element
    OUTPUT.update(d_t)

    # Train model
    model.fit(X_train, y_train)

    # Calculate model performance and organise them to dictionary
    Test_Scores_Dict1 = ss.classifier_true_false_positive_negative(model, X_test, y_test, POS_LABEL, out=dict)
    for key, value in Test_Scores_Dict1.items():
        OUTPUT[f"[Test]_{key}"] = value

    FINAL_RESULT[train_test_split_iter_n] = OUTPUT
    DF_OUTPUT = pd.DataFrame.from_dict(FINAL_RESULT, "index")

    ss.database_save(study, f"{FILENAME_OUTPUT}//study//study_index_{train_test_split_iter_n}", 'Output')
    study._storage.delete_study(study._study_id)
    ss.save_file(DF_OUTPUT, f"{FILENAME_OUTPUT}//DF_OUTPUT")

    print (f"{Dataset_Note} - Results saved.")
    
    

print ('All Optimisation Done!')




