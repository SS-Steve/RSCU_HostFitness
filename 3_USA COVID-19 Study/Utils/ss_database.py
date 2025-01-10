### Import basic modules

import pandas as pd
import numpy as np


########################################################################################################################
### NCBI Download and Codon Bias ###
########################################################################################################################

codon_table_list = ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA',
                    'TAG','TGT','TGC','TGA','TGG','CTT','CTC','CTA','CTG','CCT','CCC',
                    'CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG','ATT',
                    'ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG',
                    'AGT','AGC','AGA','AGG','GTT','GTC','GTA','GTG','GCT','GCC','GCA',
                    'GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']

aa_table_list = ['S','R','L','A','G','P','T','V','I','F',
                 'N','D','H','C','Y','Q','K','E','M','W',
                 '*']

translation_aa_codon_dict = {'S':['TCT','TCC','TCA','TCG','AGT','AGC'], 
                             'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
                             'L':['CTT','CTC','CTA','CTG','TTA','TTG'], 
                             'A':['GCT','GCC','GCA','GCG'],
                             'G':['GGT','GGC','GGA','GGG'], 
                             'P':['CCT','CCC','CCA','CCG'], 
                             'T':['ACT','ACC','ACA','ACG'],
                             'V':['GTT','GTC','GTA','GTG'], 
                             'I':['ATT','ATC','ATA'], 
                             'F':['TTT','TTC'],
                             'N':['AAT','AAC'], 
                             'D':['GAT','GAC'], 
                             'H':['CAT','CAC'], 
                             'C':['TGT','TGC'],
                             'Y':['TAT','TAC'], 
                             'Q':['CAA','CAG'], 
                             'K':['AAA','AAG'], 
                             'E':['GAA','GAG'],
                             'M':['ATG'], 
                             'W':['TGG'], 
                             '*':['TAA','TAG','TGA']}

translation_codon_aa_dict = {'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                             'AGA':'R','AGG':'R','CTT':'L','CTC':'L','CTA':'L','CTG':'L','TTA':'L','TTG':'L','GCT':'A','GCC':'A',
                             'GCA':'A','GCG':'A','GGT':'G','GGC':'G','GGA':'G','GGG':'G','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                             'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GTT':'V','GTC':'V','GTA':'V','GTG':'V','ATT':'I','ATC':'I',
                             'ATA':'I','TTT':'F','TTC':'F','AAT':'N','AAC':'N','GAT':'D','GAC':'D','CAT':'H','CAC':'H','TGT':'C',
                             'TGC':'C','TAT':'Y','TAC':'Y','CAA':'Q','CAG':'Q','AAA':'K','AAG':'K','GAA':'E','GAG':'E','ATG':'M',
                             'TGG':'W','TAA':'*','TAG':'*','TGA':'*'}

empty_codon_usage = {'TTT': 0.0, 'TTC': 0.0, 'TTA': 0.0, 'TTG': 0.0, 'TCT': 0.0, 
                     'TCC': 0.0, 'TCA': 0.0, 'TCG': 0.0, 'TAT': 0.0, 'TAC': 0.0, 
                     'TAA': 0.0, 'TAG': 0.0, 'TGT': 0.0, 'TGC': 0.0, 'TGA': 0.0, 
                     'TGG': 0.0, 'CTT': 0.0, 'CTC': 0.0, 'CTA': 0.0, 'CTG': 0.0, 
                     'CCT': 0.0, 'CCC': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CAT': 0.0, 
                     'CAC': 0.0, 'CAA': 0.0, 'CAG': 0.0, 'CGT': 0.0, 'CGC': 0.0, 
                     'CGA': 0.0, 'CGG': 0.0, 'ATT': 0.0, 'ATC': 0.0, 'ATA': 0.0, 
                     'ATG': 0.0, 'ACT': 0.0, 'ACC': 0.0, 'ACA': 0.0, 'ACG': 0.0, 
                     'AAT': 0.0, 'AAC': 0.0, 'AAA': 0.0, 'AAG': 0.0, 'AGT': 0.0, 
                     'AGC': 0.0, 'AGA': 0.0, 'AGG': 0.0, 'GTT': 0.0, 'GTC': 0.0, 
                     'GTA': 0.0, 'GTG': 0.0, 'GCT': 0.0, 'GCC': 0.0, 'GCA': 0.0, 
                     'GCG': 0.0, 'GAT': 0.0, 'GAC': 0.0, 'GAA': 0.0, 'GAG': 0.0, 
                     'GGT': 0.0, 'GGC': 0.0, 'GGA': 0.0, 'GGG': 0.0}

multibox_codon_table_list = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 
                             'TGT', 'TGC', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 
                             'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 
                             'ATA', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 
                             'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 
                             'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']

########################################################################################################################
### Model training functions ###
########################################################################################################################

classifier_shortname_fullname_dict = {
    'AB'       : 'AdaBoost',
    'DT'       : 'Decision_Tree',
    'ET'       : 'Extra_Trees',
    'Gaussian' : 'Gaussian_Process',
    'KNN'      : 'K_Neighbors',
    'LDA'      : 'Linear_Discriminant_Analysis',
    'LGBM'     : 'LightGBM',
    'LogReg'   : 'LogisticRegression',
    'LSVC'     : 'LinearSVM',
    'MNB'      : 'Multinomial_Naive_Bayes',
    'NB'       : 'Gaussian_Naive_Bayes',
    'QDA'      : 'Quadratic_Discriminant_Analysis',
    'RF'       : 'Random_Forest',
    'Ridge'    : 'Ridge',
    'SGD'      : 'Stochastic_Gradient_Descent',
    'SVC'      : 'SVM',
    'XGB'      : 'XGBoost',
    'MLC'      : 'Neural_Network',
    }

Logistic_Regression_parameter_grid = {
    'penalty'       : ('l1', 'l2', 'elasticnet', 'none'),
    'dual'          : (True, False),
    'tol'           : list(np.logspace(-8, 2, 11)),
    'C'             : [0.1] + list(np.arange(0.5, 3.1, 0.5)),
    'fit_intercept' : (True, False),
    'solver'        : ('newton-cg','lbfgs', 'liblinear','sag', 'saga')
    }

Ridge_parameter_grid = {
    'alpha'         : list(np.arange(0.1, 0.21, 0.01)),
    'fit_intercept' : (True, False),
    'normalize'     : (True, False),
    'tol'           : list(np.logspace(-6,2, 9)),
    'solver'        : ('auto', 'svd','cholesky','lsqr','sparse_cg','sag','saga','lbfgs'),
    'positive'      : (True, False),
    }

Decision_Tree_parameter_grid = {
    'criterion'         : ('gini', 'entropy'),
    'splitter'          : ('best','random'),
    'max_depth'         : list(range(5, 105, 10)) + [None],
    'min_samples_split' : list(range(2, 53, 3)),
    'min_samples_leaf'  : list(range(1, 51, 3)),
    'max_features'      : ('auto', 'sqrt', 'log2'),
    'class_weight'      : ['balanced']
    }

Random_Forest_parameter_grid = {
    'n_estimators'      : [1, 5] + list(range(10, 230, 30)),
    'criterion'         : ('gini', 'entropy'),
    'max_depth'         : list(range(5, 41, 10)) + [None],
    'min_samples_split' : list(range(2, 18, 3)),
    'min_samples_leaf'  : list(range(1, 18, 3)),
    'max_features'      : ('sqrt','log2'),
    'class_weight'      : ['balanced'],
    }

Extra_Tree_parameter_grid = {
    'criterion'         : ('gini','entropy','log_loss'),
    'splitter'          : ('random', 'best'),
    'min_samples_split' : list(np.arange(0.5, 5.5, 0.5)),
    'min_samples_leaf'  : list(range(1,11,1)),
    'max_features'      : ('auto','sqrt','log2'),
    }

KNearestNeighbors_parameter_grid = {
    'n_neighbors' : [1] + list(range(5, 105, 5)),
    'weights'     : ('uniform','distance'),
    'algorithm'   : ('auto', 'ball_tree', 'kd_tree', 'brute'),
    'leaf_size'   : list(range(10, 110, 10)),
    'p'           : [1, 2],
    }

Gaussian_parameter_grid = {
    'n_restarts_optimizer' : list(range(0, 11, 1)),
    'max_iter_predict'     : list(range(0, 201, 50)),
    }

Gaussian_NaiveBayes_parameter_grid = {
    'var_smoothing' : list(np.logspace(100,-9, num=100))
    }

Multinomial_NaiveBayes_parameter_grid = {
    'alpha'     : list(np.arange(0, 5.1, 0.1)),
    'fit_prior' : (True, False),
    }

LinearSVC_parameter_grid = {
    'penalty'      : ('l1', 'l2'),
    'loss'         : ('hinge', 'squared_hinge'),
    'dual'         : (True, False),
    'tol'          : list(np.logspace(-1, -10, 10)),
    'C'            : list(np.arange(0.2, 2.2, 0.2)),
    'multi_class'  : ('ovr', 'crammer_singer'),
    'max_iter'     : [10000],
    }

SVC_parameter_grid = {
    'C'                       : list(range(1,21,2)),
    'kernel'                  : ('linear','poly','rbf','sigmoid'),
    'gamma'                   : ('scale','auto'),
    'tol'                     : [0.0001],
    'decision_function_shape' : ('ovo','ovr')
    }

LDA_parameter_grid = {
    'solver'    : ('svd', 'lsqr','eigen'),
    'shrinkage' : (None, 'auto'),
    'tol'       : list(np.logspace(-10, 1, 100))
    }

QDA_parameter_grid = {
    'tol' : list(np.logspace(-10, 1, 100)),
    }

Stochastic_Gradient_Descent_parameter_grid = {
    'loss'       : ('hinge','log_loss','log','modified_huber','squared_hinge',
                    'perceptron','squared_error','huber','epsilon_insensitive','squared_epsilon_insensitive'),
    'penalty'    : ('l2', 'l1','elasticnet'),
    'alpha'      : list(np.logspace(-10, 0, 11)),
    'tol'        : list(np.logspace(-8, 2, 11)),
    }

AdaBoost_parameter_grid = {
    'n_estimators'  : list(range(10, 210, 10)),
    'learning_rate' : list(np.arange(0.5, 3.5, 0.5)),
    'algorithm'     : ('SAMME', 'SAMME.R'),
    }

LGBM_parameter_grid = {
    'boosting_type' : ('gbdt','dart','log','goss','rf',
                       'perceptron','squared_error','huber',
                       'epsilon_insensitive','squared_epsilon_insensitive'),
    'num_leaves'    : list(range(5, 55, 5)),
    'max_depth'     : list(range(5, 41, 10)) + [-1],
    'n_estimators'  : [1, 5] + list(range(10, 230, 30)),
    }

XGB_parameter_grid = {
    'booster'      : ('gbtree','gblinear','dart'),
    'num_leaves'   : list(range(5, 55, 5)),
    'max_depth'    : list(range(5, 41, 10)) + [-1],
    'n_estimators' : list(range(200, 530, 30)),
}
