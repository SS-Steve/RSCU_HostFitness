
### This is for Supplemental Figure 2B.

### System requirement ###

import numpy as np
import pandas as pd
from scipy import stats
import sklearn
from tqdm.notebook import tqdm as tb
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
    #'TAXONOMY':TAXONOMY,
    'START_STOP':START_STOP,
    'HS_CORR':HS_CORR,
    'HS_CORR_AA':HS_CORR_AA,
    }

#Dataset_Note = "RSCU"
Base_X = pd.concat([RSCU, TAXONOMY], axis=1)

#print (f'Dataset: {Dataset_Note}')

    
# -------------------------------------------------------------------------------- #
# Define y
animal_list = ["human", "vertebrates", "land plants", "invertebrates", "bacteria"]

POS_LABEL_Raw = "human"
# -------------------------------------------------------------------------------- #

Study_Title = 'ExtraFeatureWithTaxonomy'

CPU_USED = os.cpu_count() - 5
#CPU_USED = int(os.cpu_count() / 2)
#CPU_USED     = 5

CV = 20
#TEST_SIZE = 0.7
RANDOM_STATE = 8
Optuna_Trial_N = 50

MODELNAME = "RF"
UPSAMPLE_METHOD = "SMOTE"

#Optuna_CV_Score = 'recall'
CV_Scoring_List = ['average_precision', 'roc_auc', 'precision', 'recall', 'f1', 
                   'accuracy', 'balanced_accuracy']

#Run_List_Test_Size = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
Run_List_Test_Size = [0.1, 0.2, 0.3]

for TEST_SIZE in Run_List_Test_Size:

    Dict_Optuna_CV_Score_OutputName = {'balanced_accuracy':'BA',
                                       'f1':'F1',
                                       'recall':'Recall'}

    Run_List_Optuna_CV_Score = ['balanced_accuracy', 'f1', 'recall']

    for Optuna_CV_Score in Run_List_Optuna_CV_Score:
        
        print (f"###################################################")
        print (f"Test_Size : {TEST_SIZE} - Score : {Optuna_CV_Score}")
        print (f"###################################################")

        # Custom output note on output folder
        NOTING = f"[{POS_LABEL_Raw}]_for{Dict_Optuna_CV_Score_OutputName[Optuna_CV_Score]}"
        OUTPUT_ADD_NOTE = (
            f"{Study_Title}_{TEST_SIZE}_RS{RANDOM_STATE}_{UPSAMPLE_METHOD}_{NOTING}"
        )

        y = HOST[POS_LABEL_Raw]

        # Numerise y to 0 and 1
        POS_LABEL = POS_LABEL_Raw
        NUMERIC_REF = {POS_LABEL: 1, f"not_{POS_LABEL}": 0}
        y = ss.df_col_numerise(y, NUMERIC_REF)
        POS_LABEL = NUMERIC_REF[POS_LABEL]

        y_raw = y.copy()
        y = np.array(y)

        # Create Output folder
        FILENAME_OUTPUT = f"[{MODELNAME}]_{OUTPUT_ADD_NOTE}"
        ss.create_folder(FILENAME_OUTPUT, "Output")

        FINAL_RESULT = {}
        
        Run_List_Dataset = [None] + list(X_Dict_Feature.keys())

        for Index, ExtraFeature in enumerate(Run_List_Dataset):
            
            if ExtraFeature is None:
                Dataset_Note = f"RSCUpTAXONOMY"
                print (f'Dataset: {Dataset_Note}')
                X = pd.concat([Base_X], axis=1)
            else:
                Dataset_Note = f"RSCUpTAXONOMYp{ExtraFeature}"
                print (f'Dataset: {Dataset_Note}')
                X = pd.concat([Base_X, X_Dict_Feature[ExtraFeature]], axis=1)

            X_raw = X.copy()
            X = np.array(X)
            
            feautre_list = X_raw.columns.tolist()

            ### Train-test split ###
            from sklearn.model_selection import train_test_split
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=TEST_SIZE, shuffle=True, random_state=RANDOM_STATE
            )
                
            print(f"{Dataset_Note} - Finish Train-test split.")

            if UPSAMPLE_METHOD in ["SMOTE"]:
                ### SMOTE upsampling ###
                from imblearn.over_sampling import SMOTE
                sm = SMOTE(random_state=RANDOM_STATE)
                X_train, y_train = sm.fit_resample(X_train, y_train)
                
                
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
                      'Extra_Feature':ExtraFeature,
                      'POS_LABEL':POS_LABEL_Raw,
                      'CV':CV,
                      'Test_Ratio':TEST_SIZE,
                      'Train_Ratio':1-TEST_SIZE,
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
            Test_Scores_Dict2 = ss.classifier_test_scores(model, X_test, y_test, POS_LABEL, out=dict)
            for key, value in Test_Scores_Dict2.items():
                OUTPUT[f"[Test]_{key}"] = value

            FINAL_RESULT[Index] = OUTPUT
            DF_OUTPUT = pd.DataFrame.from_dict(FINAL_RESULT, "index")

            ss.database_save(study, f"{FILENAME_OUTPUT}//study_{Dataset_Note}", 'Output')
            study._storage.delete_study(study._study_id)
            ss.save_file(DF_OUTPUT, f"{FILENAME_OUTPUT}//DF_OUTPUT")

            print (f"{Dataset_Note} - Results saved.")
            
        

print ('All Optimisation Done!')




