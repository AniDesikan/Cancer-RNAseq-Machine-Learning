#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import gzip
import requests
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report

#Initialize all the file paths for the csvs containing cancer RNA sequencing information
#The TCGA datasets were taken from the Xena UCSC data portal
colon_phenotype = "TCGA-COAD.GDC_phenotype.tsv.gz"
lung_phenotype = "TCGA-LUAD.GDC_phenotype.tsv.gz"
colon_genes = "colon_tm.csv"
lung_genes = "lung_tm.csv"
colon_sequencing = "TCGA-COAD.htseq_fpkm.tsv.gz"
lung_sequencing = "TCGA-LUAD.htseq_fpkm.tsv.gz"

#Open the phenotype csvs for sorting the samples into primary tumor vs not cancer later on
with gzip.open(colon_phenotype, 'rt') as file:
    colon_phenotype = pd.read_csv(file, sep='\t')

with gzip.open(lung_phenotype, 'rt') as file:
    lung_phenotype = pd.read_csv(file, sep='\t')

#These are the dataframes only containing enhanced genes in the colon and lung cancer from the tm_function
with open(colon_genes, 'rt') as file:
    colon_tm = pd.read_csv(file) 
    
with open(lung_genes, 'rt') as file:
    lung_tm = pd.read_csv(file)

#Full rna sequencing datasets to see how well it commpares to the shortened ones
with gzip.open(colon_sequencing, 'rt') as file:
    og_colon = pd.read_csv(file, sep='\t')
    
with gzip.open(lung_sequencing, 'rt') as file:
    og_lung = pd.read_csv(file, sep='\t')

def arrange_df(phenotype, genes):
    gene_names = genes.iloc[:,0]
    columns = phenotype['submitter_id.samples'].tolist()
    selected_columns = genes.loc[:, genes.columns.isin(columns)]
    selected_columns_t = selected_columns.transpose().reset_index()
    selected_columns_t = selected_columns_t.rename(columns={'index': 'SampleID'})
    selected_rows = phenotype[phenotype['submitter_id.samples'].isin(selected_columns.columns)]
    selected_rows = selected_rows[['submitter_id.samples', 'sample_type.samples']]
    merged_df = pd.merge(selected_columns_t, selected_rows, how='inner', left_on='SampleID', right_on='submitter_id.samples')
    merged_df_t = merged_df.T
    return(merged_df_t.iloc[:-2].append(merged_df_t.iloc[-1]).transpose())

def machine(df):
    X = df.drop(columns=['SampleID', 'sample_type.samples'])
    y = df['sample_type.samples']    
    test_size=0.5
        # Split the data into training and testing sets (e.g., 80% train, 20% test)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size, random_state=42)
    svm_model = SVC(kernel='linear', C=1.0)
    # Train the SVM model
    svm_model.fit(X_train, y_train)
    # Make predictions on the testing set
    y_pred = svm_model.predict(X_test)
    # Evaluate the model
    accuracy = accuracy_score(y_test, y_pred)
    classification_report_result = classification_report(y_test, y_pred)
    # Display results
    print('Results of testing on ' + df + 'with test size' + test_size)
    print(f'Accuracy: {accuracy}')
    print('Classification Report:')
    print(classification_report_result)
    
#    
colon_merged_df= arrange_df(colon_phenotype, colon_tm)
lung_merged_df= arrange_df(lung_phenotype, lung_tm)
total_colon_merged_df = arrange_df(colon_phenotype, og_colon)
total_lung_merged_df = arrange_df(lung_phenotype, og_lung)

# Merged 
colon_merged_df_2 = colon_merged_df[colon_merged_df["sample_type.samples"] == "Primary Tumor"]
colon_merged_df_2["sample_type.samples"] = "colon"
lung_merged_df_2 = lung_merged_df[lung_merged_df["sample_type.samples"] == "Primary Tumor"]
lung_merged_df_2["sample_type.samples"] = "lung"
merged = pd.concat([lung_merged_df_2, colon_merged_df_2], ignore_index=True)


machine(colon_merged_df)
machine(lung_merged_df)
machine(total_colon_merged_df)
machine(total_lung_merged_df)
machine(merged)


