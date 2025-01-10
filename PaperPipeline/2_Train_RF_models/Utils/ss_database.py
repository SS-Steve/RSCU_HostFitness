### Import basic modules

import pandas as pd
import numpy as np




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
