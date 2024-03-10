#################################################
# CONSTANTS
#################################################
colon_data = "/projectnb/mlresearch/test/ani/cancer/data/binary/TCGA-COAD.htseq_fpkm.tsv"
colon_phenotype = "/projectnb/mlresearch/test/ani/cancer/data/binary/TCGA-COAD.GDC_phenotype.tsv"
colon_enhanced_genes = "/projectnb/mlresearch/test/ani/cancer/diff_exp_csvs/Colon_DE.csv"

lung_data = '/projectnb/mlresearch/data/binary/Lung_Adenocarcinoma_LUAD/TCGA-LUAD.htseq_fpkm.tsv'
lung_phenotype = '/projectnb/mlresearch/data/binary/Lung_Adenocarcinoma_LUAD/TCGA-LUAD.GDC_phenotype.tsv'

prostate_data = '/projectnb/mlresearch/data/binary/Prostate_Cancer_PRAD/TCGA-PRAD.htseq_fpkm.tsv'
prostate_phenotype = '/projectnb/mlresearch/data/binary/Prostate_Cancer_PRAD/TCGA-PRAD.GDC_phenotype.tsv'

#################################################
# IMPORTS
#################################################

import pandas as pd
import gzip
import requests
import numpy as np 
import tensorflow as tf
from tensorflow.keras import models,layers,losses,metrics,optimizers,Input
from data_augmentation.tm_function import *
from machines.Machines import BinaryModels
from sklearn.model_selection import train_test_split
from sklearn.metrics import  accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, auc

#################################################
# PIPELINE FUNCTIONS
#################################################

# Get data into a dataframe (Done)
# Augment data using TM function / linear interpolation (Do after getting machines set up)
# Input into Binary machines
# Return metrics

def binary_pipeline(data, phenotype, model_type, TM_use = False, value = None, enhanced_genes = None, lin_int = False, 
                    noise = False, noise_type = None, noise_level = None):

    #TM Function if necessary (use values from before)
    if TM_use == True:
        if enhanced_genes == None or value == None:
              raise ValueError("enhanced_genes and value must not be None")
        data = TM(data, enhanced_genes, value)
    
    else:
        data = sequencing_dataframe(data)

    #Linear interpolation if necessary 
    if lin_int == True:
        data = linear_interpolation(data)
    
    # Noise injection if necessary
    if noise == True:
        if noise_type == None:
            raise ValueError("noise_type and noise_level must not be None")

        data = noise_injection(noise_type, data, noise_level)

    #First seperate the cancer from healthy, add indicators, and shuffle them back together
    cancer, healthy = separate_cancer(phenotype, data)
    cancer = cancer.sample(n=len(healthy), random_state=42)
    combined_samples = pd.concat([cancer, healthy], ignore_index=True)
    combined_samples = combined_samples.sample(frac=1).reset_index(drop=True)
    combined_samples = combined_samples.iloc[:, 1:]
    X = combined_samples.drop(columns=['sample_type.samples','id', 'is_cancer']).values.astype(np.float32)
    y = combined_samples['is_cancer']
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
    

    #Pipe resulting dataframes into svm using split and stuff
    #Use string given in function
    binary_model = BinaryModels(model_type, X_train=X_train, y_train=y_train)
    model = binary_model.get_model()

    y_pred = model.predict(X_val)

    # Compute evaluation metrics
    accuracy = accuracy_score(y_val, y_pred)
    precision = precision_score(y_val, y_pred, average='binary', pos_label=1)
    recall = recall_score(y_val, y_pred, average='binary', pos_label=1)
    f1 = f1_score(y_val, y_pred, average='binary', pos_label=1)
    metrics = {'accuracy': accuracy, 'precision': precision, 'recall': recall, 'f1': f1}

    #compile parameters
    params = {'TM': TM_use, 'value': value, 'lin_int': lin_int, 'noise': noise, 'noise_type': noise_type}
    return params, metrics

BINARY_PARAM_GRID = {
        'value': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        'noise_type': ["uniform", "salt and pepper",  "poisson", "exponential"],
    }

import random 
def select_random_binary_params():
    value = random.choice(BINARY_PARAM_GRID['value'])
    noise_type = random.choice(BINARY_PARAM_GRID['noise_type'])
    return value, noise_type
#################################################
# FULL PIPELINE
#################################################
# print("Colon")
# params, metrics = binary_pipeline(colon_data, colon_phenotype, 'svm')
# print("Params: ")
# print(params)
# print("Metrics: ")
# print(metrics)

# print("Colon TM")
# params, metrics = binary_pipeline(colon_data, colon_phenotype, 'svm', TM_use = True, value = 4, enhanced_genes = colon_enhanced_genes)
# print("Params: ")
# print(params)
# print("Metrics: ")
# print(metrics)

# print("Colon lin int")
# params, metrics = binary_pipeline(lung_data, lung_phenotype, 'svm', lin_int = True)
# print("Params: ")
# print(params)
# print("Metrics: ")
# print(metrics)

# print("Colon noise")
# params, metrics = binary_pipeline(lung_data, lung_phenotype, 'svm', noise = True, noise_type = "poisson")
# print("Params: ")
# print(params)
# print("Metrics: ")
# print(metrics)

# print("Colon all")
# params, metrics = binary_pipeline(lung_data, lung_phenotype, 'svm',  TM_use = True, value = 4, enhanced_genes = colon_enhanced_genes, lin_int = True,  noise = True, noise_type = "poisson")
# print("Params: ")
# print(params)
# print("Metrics: ")
# print(metrics)

# print("Colon all")
# params, metrics = binary_pipeline(lung_data, lung_phenotype, 'svm',  TM_use = True, value = 4, enhanced_genes = colon_enhanced_genes, noise = True, noise_type = "poisson")
# print("Params: ")
# print(params)
# print("Metrics: ")
# print(metrics)
