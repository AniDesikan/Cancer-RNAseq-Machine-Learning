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
colon_phenotype = "data/TCGA-COAD.GDC_phenotype.tsv.gz"
lung_phenotype = "data/TCGA-LUAD.GDC_phenotype.tsv.gz"
lymphoma_phenotype = "data/TCGA-DLBC.GDC_phenotype.tsv.gz"
leukemia_phenotype = "data/TCGA-LAML.GDC_phenotype.tsv.gz"
prostate_phenotype = "data/TCGA-PRAD.GDC_phenotype.tsv.gz"

colon_genes = "data/colon_tm.csv"
lung_genes = "data/lung_tm.csv"
lymphoma_genes = "data/lymphoma_tm.csv"
leukemia_genes = "data/leukemia_tm.csv"
prostate_genes = "data/prostate_tm.csv"

colon_sequencing = "data/TCGA-COAD.htseq_fpkm.tsv.gz"
lung_sequencing = "data/TCGA-LUAD.htseq_fpkm.tsv.gz"
lymphoma_sequencing = "data/TCGA-DLBC.htseq_fpkm.tsv.gz"
leukemia_sequencing = "data/TCGA-LAML.htseq_fpkm.tsv.gz"
prostate_sequencing = "data/TCGA-PRAD.htseq_fpkm.tsv.gz"

#Open the phenotype csvs for sorting the samples into primary tumor vs not cancer later on
with gzip.open(colon_phenotype, 'rt') as file:
    colon_phenotype = pd.read_csv(file, sep='\t')

with gzip.open(lung_phenotype, 'rt') as file:
    lung_phenotype = pd.read_csv(file, sep='\t')

with gzip.open(lymphoma_phenotype, 'rt') as file:
    lymphoma_phenotype = pd.read_csv(file, sep='\t')

with gzip.open(leukemia_phenotype, 'rt') as file:
    leukemia_phenotype = pd.read_csv(file, sep='\t')

with gzip.open(prostate_phenotype, 'rt') as file:
    prostate_phenotype = pd.read_csv(file, sep='\t')

#These are the dataframes only containing enhanced genes in the colon and lung cancer from the tm_function
# with open(colon_genes, 'rt') as file:
#     colon_tm = pd.read_csv(file) 
    
# with open(lung_genes, 'rt') as file:
#     lung_tm = pd.read_csv(file)
    
# with open(lymphoma_genes, 'rt') as file:
#     lymphoma_tm = pd.read_csv(file)

# with open(leukemia_genes, 'rt') as file:
#     leukemia_tm = pd.read_csv(file)

# with open(prostate_genes, 'rt') as file:
#     prostate_tm = pd.read_csv(file)

#Full rna sequencing datasets to see how well it commpares to the shortened ones
with gzip.open(colon_sequencing, 'rt') as file:
    colon_sequencing = pd.read_csv(file, sep='\t')

with gzip.open(lung_sequencing, 'rt') as file:
    lung_sequencing = pd.read_csv(file, sep='\t')

with gzip.open(lymphoma_sequencing, 'rt') as file:
    lymphoma_sequencing = pd.read_csv(file, sep='\t')

with gzip.open(leukemia_sequencing, 'rt') as file:
    leukemia_sequencing = pd.read_csv(file, sep='\t')

with gzip.open(prostate_sequencing, 'rt') as file:
    prostate_sequencing = pd.read_csv(file, sep='\t')

#This is a function to arrange the rna sequencing dataframe in a way that is easy to put into the machine function
def arrange_df(phenotype, genes):
    #Don't need the gene names for the machine function but store them here in case they're needed later
    gene_names = genes.iloc[:,0]
    #The following lines format the dataframe and inner merge the phenotype samples and sequencing samples so only the samples that have phenotypes are returned
    columns = phenotype['submitter_id.samples'].tolist()
    selected_columns = genes.loc[:, genes.columns.isin(columns)]
    selected_columns_t = selected_columns.transpose().reset_index()
    #Rename the index column to sampleID after resetting the index
    selected_columns_t = selected_columns_t.rename(columns={'index': 'SampleID'})
    selected_rows = phenotype[phenotype['submitter_id.samples'].isin(selected_columns.columns)]
    selected_rows = selected_rows[['submitter_id.samples', 'sample_type.samples']]
    merged_df = pd.merge(selected_columns_t, selected_rows, how='inner', left_on='SampleID', right_on='submitter_id.samples')
    #Since we're selecting the columns, we had to transpose, transpose back
    merged_df_t = merged_df.T
    #Return the dataframe while getting rid of the second to last column since we don't want it anymore, the last one is good enough
    return(merged_df_t.iloc[:-2]._append(merged_df_t.iloc[-1]).transpose())

#Function that splits, trains and test the reaaranged sequencing dataframe on an SVM model and then prints the results
def machine(df):
    y = df['sample_type.samples']  
    X = df.drop(columns=['SampleID', 'sample_type.samples'])
    test_size=0.5
        # Split the data into training and testing sets (e.g., 80% train, 20% test) according to the test_size variable
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    svm_model = SVC(kernel='linear', C=1.0)
    # Train the SVM model
    svm_model.fit(X_train, y_train)
    # Make predictions on the testing set
    y_pred = svm_model.predict(X_test)
    # Evaluate the model
    accuracy = accuracy_score(y_test, y_pred)
    classification_report_result = classification_report(y_test, y_pred)
    # Display results
    print('Results of testing with test size ' + str(test_size))
    print(f'Accuracy: {accuracy}') 
    print('Classification Report:')
    print(classification_report_result)
    
#Add controls and arrange for the machine for leukemia and lymphoma, since there aren't any healthy controls in the datasets to begin with
#TODO: Make this a function
with open("data/healthy_blood_rnaseq.csv", 'rt') as file:
    healthy_blood = pd.read_csv(file)
new_column_names = ['Ensembl_ID', 'control_1', 'control_2','control_3','control_4','control_5','control_6','control_7','control_8','control_9',]
healthy_blood.columns = new_column_names
leukemia_w_controls = pd.merge(leukemia_sequencing, healthy_blood, how='inner', left_on='Ensembl_ID', right_on='Ensembl_ID')
lymphoma_w_controls = pd.merge(lymphoma_sequencing, healthy_blood, how='inner', left_on='Ensembl_ID', right_on='Ensembl_ID')
lymphoma_column_names = lymphoma_w_controls.columns
lymphoma_row_values = lymphoma_column_names.map(lambda x: "healthy" if "control" in x else "tumor")

# Add the new row to the DataFrame
lymphoma_w_controls.loc['sample_type.samples'] = lymphoma_row_values
lymphoma = lymphoma_w_controls.T

leukemia_column_names = leukemia_w_controls.columns
leukemia_row_values = leukemia_column_names.map(lambda x: "healthy" if "control" in x else "tumor")

# Add the new row to the DataFrame
leukemia_w_controls.loc['sample_type.samples'] = leukemia_row_values
leukemia = leukemia_w_controls.T
# leukemia = leukemia.drop('new_row', axis=1)

lymphoma.drop('Ensembl_ID', axis=0, inplace=True)
leukemia.drop('Ensembl_ID', axis=0, inplace=True)
leukemia['SampleID'] = leukemia.index

# Reset the index to default integer index
leukemia.reset_index(drop=True, inplace=True)
lymphoma['SampleID'] = lymphoma.index

# Reset the index to default integer index
lymphoma.reset_index(drop=True, inplace=True)




# Arrange the rna sequencing dataframes that are shortened and the full ones
# colon_merged_df= arrange_df(colon_phenotype, colon_tm)
# lung_merged_df= arrange_df(lung_phenotype, lung_tm)
# lymphoma_merged_df= arrange_df(lymphoma_phenotype, lung_tm)
# leukemia_merged_df= arrange_df(leukemia_phenotype, lung_tm)
# prostate_merged_df= arrange_df(prostate_phenotype, lung_tm)

total_colon_merged_df = arrange_df(colon_phenotype, colon_sequencing)
total_lung_merged_df = arrange_df(lung_phenotype, lung_sequencing)
# For leukemia and lymphoma add the controls and arrange in a seperate function 
total_lymphoma_merged_df= lymphoma
total_leukemia_merged_df= leukemia
total_prostate_merged_df= arrange_df(prostate_phenotype, prostate_sequencing)


# See if the machine is able to tell the difference between the colon and lung datasets
# colon_merged_df_2 = colon_merged_df[colon_merged_df["sample_type.samples"] == "Primary Tumor"]
# colon_merged_df_2["sample_type.samples"] = "colon"
# lung_merged_df_2 = lung_merged_df[lung_merged_df["sample_type.samples"] == "Primary Tumor"]
# lung_merged_df_2["sample_type.samples"] = "lung"
# merged = pd.concat([lung_merged_df_2, colon_merged_df_2], ignore_index=True)

#Run the machine function
# print("colon filtered")
# machine(colon_merged_df)
# print("lung filtered")
# machine(lung_merged_df)
# print("lymphoma filtered")
# machine(lymphoma_merged_df)
# print("leukemia filtered")
# machine(leukemia_merged_df)
# print("prostate filtered")
# machine(prostate_merged_df)

print("colon unfiltered")
machine(total_colon_merged_df)
print("lung unfiltered")
machine(total_lung_merged_df)
print("lymphoma unfiltered")
machine(total_lymphoma_merged_df)
print("leukemia unfiltered")
machine(total_leukemia_merged_df)
print("prostate unfiltered")
machine(total_prostate_merged_df)


