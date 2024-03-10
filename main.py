#################################################
# CONSTANTS
#################################################

colon_data = "/projectnb/mlresearch/data/binary/Colon_Cancer_COAD/TCGA-COAD-Copy1.htseq_fpkm.tsv"
colon_phenotype = "/projectnb/mlresearch/data/binary/Colon_Cancer_COAD/TCGA-COAD-Copy1.GDC_phenotype.tsv"
colon_enhanced_genes = "/projectnb/mlresearch/test/ani/cancer/diff_exp_csvs/Colon_DE.csv"

lung_data = "/projectnb/mlresearch/data/binary/Lung_Adenocarcinoma_LUAD/TCGA-LUAD.htseq_fpkm.tsv"
lung_phenotype = "/projectnb/mlresearch/data/binary/Lung_Adenocarcinoma_LUAD/TCGA-LUAD.GDC_phenotype.tsv"
lung_enhanced_genes = "/projectnb/mlresearch/test/ani/cancer/diff_exp_csvs/Lung_DE.csv"

prostate_data = "/projectnb/mlresearch/test/ani/cancer/data/binary/TCGA-COAD.htseq_fpkm.tsv"
prostate_phenotype = "/projectnb/mlresearch/test/ani/cancer/data/binary/TCGA-COAD.GDC_phenotype.tsv"
prostate_enhanced_genes = "/projectnb/mlresearch/test/ani/cancer/diff_exp_csvs/Prostate_DE.csv"

oral_normal = "/projectnb/mlresearch/test/ani/cancer/data/Image/oral_normal"
oral_scc = "/projectnb/mlresearch/test/ani/cancer/data/Image/oral_scc"

all_normal= "/projectnb/mlresearch/data/image/ALL/ALL_Benign"
all_scc= "/projectnb/mlresearch/data/image/ALL/ALL_Pre_B"

breast_normal= "/projectnb/mlresearch/data/image/Breast/Breast_Benign"
breast_scc= "/projectnb/mlresearch/data/image/Breast/Breast_Malignant"

colon_normal= "/projectnb/mlresearch/data/image/Colon/Colon_Benign_Tissue"
colon_scc= "/projectnb/mlresearch/data/image/Colon/Colon_Adenocarcinoma"

lung_normal= "/projectnb/mlresearch/data/image/Lung/Lung_Benign_Tissue"
lung_scc= "/projectnb/mlresearch/data/image/Lung/Lung_Adenocarcinoma"


#################################################
# IMPORTS
#################################################


import sys
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
from binary_pipeline_functions import *
import gc
import cv2
import itertools
import random
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, roc_curve, auc
import pandas as pd
import numpy as np
import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import glob
from PIL import Image
from sklearn.metrics import accuracy_score
from machines.Xception import Xception, BasicConv, ResidualUnit, MiddleFlowUnit, ExitFlow
import tensorflow as tf
from tensorflow.keras import models,layers,losses,metrics,optimizers,Input
gpus = tf.config.experimental.list_physical_devices('GPU')
print(gpus)
if gpus:
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
from machines.Xception import Xception
from machines.ResNet import ResNet_152
from machines.DenseNet import DenseNet_264
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.metrics import accuracy_score
from data_augmentation.image_processor_all import ImageProcessor,ImageProcessor2,augment_images_with_proportions_2,augment_images_with_proportions,combined_augment_function
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score
from Image_pipeline_functions import *
from features.image_features import *
# main file
# Our main file will be the control mechanism and part of creating more strategies according to the accuracy of the model.
# we will test the whole pillars using this main function.

#################################################
# MAIN FUNCTIONS
#################################################

#Here cancer is a list of which cancers you've selected
def main(cancer_list, binary_boolean, image_boolean, model_type, TM_use=False,
         lin_int=False, noise=False, augmentation=False, features=False):
    print("MAIN")
    params = None
    metrics = None
    cancer_list = cancer_list.split(',')
    
    if binary_boolean == "True":
        print("BINARY")
        for cancer in cancer_list:
            if cancer.lower() == "colon":
                data = colon_data
                phenotype = colon_phenotype
                enhanced_genes = colon_enhanced_genes
            elif cancer.lower() == "lung":
                data = lung_data
                phenotype = lung_phenotype
                enhanced_genes = lung_enhanced_genes
            elif cancer.lower() == "prostate":
                data = prostate_data
                phenotype = prostate_phenotype
                enhanced_genes = prostate_enhanced_genes
            else:
                raise ValueError(f"Unsupported cancer type: {cancer}")
                
            if data is None or phenotype is None:
                raise ValueError("You must enter phenotype and RNA sequencing data.")
            else:
                value, noise_type = select_random_binary_params()
                params, metrics = binary_pipeline(data, phenotype, model_type, TM_use, value, enhanced_genes, lin_int,
                                                       noise, noise_type)
                # TODO: This only returns the last cancer entered right now, have to fix that

    if image_boolean == "True":
        print("IMAGE")
        for cancer in cancer_list:
            if cancer.lower() == "oral":
                folder1 = oral_normal
                folder2 = oral_scc
            elif cancer.lower() == "all":
                folder1 = all_normal
                folder2 = all_scc
            elif cancer.lower() == "breast":
                folder1 = breast_normal
                folder2 = breast_scc
            elif cancer.lower() == "colon":
                folder1 = colon_normal
                folder2 = colon_scc
            elif cancer.lower() == "lung":
                folder1 = lung_normal
                folder2 = lung_scc
            else:
                raise ValueError(f"Unsupported cancer type: {cancer}")

            if folder1 is None or folder2 is None:
                raise ValueError("You must enter Image Data")
            
            participation_factors = generate_random_participation_factors()
            random_features = generate_random_features(participation_factors)
            normal_df, scc_df, params = pipeline_image_df(folder1, folder2, augmentation, num_images=5)
            
            if features:
                if participation_factors is None:
                    raise ValueError("You must enter participation factors")
                normal_df, scc_df, params = pipeline_features(normal_df, scc_df, participation_factors, params=random_features)
            
            metrics = pipeline_machine(normal_df, scc_df, model_type, epochs=1)

    return params, metrics

# Function Use Examples (in command line):

# Colon Binary on svm machine
#python main.py "Colon" True False "svm"

# Colon Binary, Lung Binary, Prostate Binary on svm with linear interpolation
# python main.py "Colon,Lung,Prostate" True False "svm" lin_int

# Colon 
# python main.py "Colon" True False "knn" TM_use noise

# Oral image on Xception with augmentation and features
# python main.py "Oral" False True Xception augmentation features
        

#################################################
# CALL MAIN
#################################################

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python main.py cancer_list binary image ...")
        sys.exit(1)
    
    cancer_list = sys.argv[1]
    binary = sys.argv[2]
    image = sys.argv[3]
    model_type = sys.argv[4]

    TM_use = False
    lin_int = False
    noise = False
    augmentation = False
    features = False


    for arg in sys.argv[5:]:
        if arg.lower() == "tm_use":
            TM_use = True
        elif arg.lower() == "lin_int":
            lin_int = True
        elif arg.lower() == "noise":
            noise = True
        elif arg.lower() == "augmentation":
            augmentation = True
        elif arg.lower() == "features":
            features = True

    print(cancer_list, binary, image, model_type, TM_use, lin_int, noise, augmentation, features)
    print(main(cancer_list, binary, image, model_type, TM_use, lin_int, noise, augmentation, features))
