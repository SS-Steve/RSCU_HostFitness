Name

2. Train_RF_models


Description

The folder '2. Train_RF_models' is mainly for training Random Forest models with RTG_Data from '1. Download NCBI accession'. The scripts
are .py file instead of .ipynb files because of using High-Performance-Computing.


Script files:

2A_TestRatioGradient.py

    This code file is for training and optimising Random Forest models with Optuna Hyper-parameters tunning. And the models are trained 
    with RSCU data to predict Host data including human, vertebrates, invertebrates, land plants and bacteria, and to test the impacts 
    of different train-test-split ratio (Train data ratio) in trained model performances. The results are shown in Figure 2A.

SF2A_ExtraFeature.py

    This code file is for training and optimising Random Forest models with Optuna Hyper-parameters tunning. Extra features or extra datasets
    are adding separately into RSCU dataset, or Dataset-R, to predict Host data of human. This tests the impacts of different train-test-split 
    ratio (Train data ratio) in trained model performances. The results are shown in Supplemental Figure 2A.

SF2B_ExtraFeature_withTaxonomy.py

    This code file is for training and optimising Random Forest models with Optuna Hyper-parameters tunning. Extra features or extra datasets
    are adding separately into combination of RSCU dataset and Taxonomy dataset, or Dataset-RT, to predict Host data of human. This tests the 
    impacts of different train-test-split ratio (Train data ratio) in trained model performances. The results are shown in Supplemental Figure 2B.

2B_TestRatio01_AllHost.py

    This code file is for training and optimising Random Forest models with Optuna Hyper-parameters tunning. And the models are trained
    with different combined datasets including Dataset-R, Dataset-RT, Dataset-RTC to predict Host data including human, vertebrates, 
    invertebrates, land plants and bacteria, and to test the models performances. The results are shown in Figure 2B.

3A_LeaveOneOut.py

    This code file is for training and optimising Random Forest models with Optuna Hyper-parameters tunning. And the models are trained
    with different combined datasets including Dataset-R, Dataset-RT, Dataset-RTC to predict Host data of human. Leave-One-Out train-test-split
    method is used in train-test-split to test performances in predicting one viruses with data from other viruses. The results are shown in Figure 3A.

SF3A_LOO_False_Samples_ReRun.py

    This code file is for training and optimising Random Forest models with Optuna Hyper-parameters tunning. This re-do the training in code file 
    '3A_LeaveOneOut.py' for false-predicted samples, and the Optuna tunning trails number is increased into 50. The results are shown in Supplemental Figure 3A.

4A_Final_RF_Model_Train.py

    This code file is for training the final Random Forest models with Optuna Hyper-parameters tunning. The models are trained with all samples
    without train-test-split, and the models are used to predict unknown viruses, which are used in the USA COVID-19 HIP study and the evolution
    path simulation study. The results are not shown, but the models are used for generating results in Figure 4 and Figure 5.

    
    
Authors

	Shuquan Su

		Email:
			Shuquan.Su@student.uts.edu.au

		Address:
			Data Science Institute (DSI), 
			School of Computer Science, 
			Faculty of Engineering and Information Technology (FEIT), 
			University of Technology Sydney (UTS)
			Sydney, New South Wales, Australia

