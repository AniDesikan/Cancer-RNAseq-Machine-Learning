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
    df = pd.read_csv(file, sep='\t')

with gzip.open(lung_phenotype, 'rt') as file:
    lung_df = pd.read_csv(file, sep='\t')

#These are the dataframes only containing enhanced genes in the colon and lung cancer from the tm_function
with open(colon_genes, 'rt') as file:
    colon_tm = pd.read_csv(file) 
    
with open(lung_genes, 'rt') as file:
    lung_tm = pd.read_csv(file)

#Full rna sequencing datasets to see how well 
with gzip.open(colon_sequencing, 'rt') as file:
    og_colon = pd.read_csv(file, sep='\t')
    
with gzip.open(lung_sequencing, 'rt') as file:
    og_lung = pd.read_csv(file, sep='\t')


gene_names = colon_tm.iloc[:,0]
columns = df['submitter_id.samples'].tolist()
selected_columns = colon_tm.loc[:, colon_tm.columns.isin(columns)]
selected_columns_t = selected_columns.transpose().reset_index()
selected_columns_t = selected_columns_t.rename(columns={'index': 'SampleID'})
selected_rows = df[df['submitter_id.samples'].isin(selected_columns.columns)]
selected_rows = selected_rows[['submitter_id.samples', 'sample_type.samples']]
merged_df = pd.merge(selected_columns_t, selected_rows, how='inner', left_on='SampleID', right_on='submitter_id.samples')
merged_df_t = merged_df.T
colon_merged_df_2 = merged_df_t.iloc[:-2].append(merged_df_t.iloc[-1]).transpose()


gene_names = lung_tm.iloc[:,0]
columns = lung_df['submitter_id.samples'].tolist()
selected_columns = lung_tm.loc[:, lung_tm.columns.isin(columns)]
selected_columns_t = selected_columns.transpose().reset_index()
selected_columns_t = selected_columns_t.rename(columns={'index': 'SampleID'})
selected_rows = lung_df[lung_df['submitter_id.samples'].isin(selected_columns.columns)]
selected_rows = selected_rows[['submitter_id.samples', 'sample_type.samples']]
merged_df = pd.merge(selected_columns_t, selected_rows, how='inner', left_on='SampleID', right_on='submitter_id.samples')
merged_df_t = merged_df.T
lung_merged_df_2 = merged_df_t.iloc[:-2].append(merged_df_t.iloc[-1]).transpose()



colon_merged_df_2 = colon_merged_df_2[colon_merged_df_2["sample_type.samples"] == "Primary Tumor"]
colon_merged_df_2["sample_type.samples"] = "colon"
lung_merged_df_2 = lung_merged_df_2[lung_merged_df_2["sample_type.samples"] == "Primary Tumor"]
lung_merged_df_2["sample_type.samples"] = "lung"
merged = pd.concat([lung_merged_df_2, colon_merged_df_2], ignore_index=True)


colon_samples = og.columns[1:].tolist()
lung_samples = og_lung.columns[1:].tolist()
merged_og = pd.merge(og, og_lung, how='inner', left_on='Ensembl_ID', right_on='Ensembl_ID')
merged_og["sample_type.samples"] = np.where(merged_og.isin(colon_samples), "colon", np.where(merged_og.columns,isin(lung_samples), "lung", ""))



merged_og



#From the model with selected genes
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report

# Assuming df is your Pandas DataFrame
# Replace 'TargetColumnName' with the actual name of your target column
X = colon_merged_df_2.drop(columns=['SampleID', 'sample_type.samples'])
y = colon_merged_df_2['sample_type.samples']

# Split the data into training and testing sets (e.g., 80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.8, random_state=42)

svm_model = SVC(kernel='linear', C=1.0)

# Train the SVM model
svm_model.fit(X_train, y_train)

# Make predictions on the testing set
y_pred = svm_model.predict(X_test)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
classification_report_result = classification_report(y_test, y_pred)

# Display results
print(f'Accuracy: {accuracy}')
print('Classification Report:')
print(classification_report_result)




#From the model with selected genes
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report

# Assuming df is your Pandas DataFrame
# Replace 'TargetColumnName' with the actual name of your target column
X = lung_merged_df_2.drop(columns=['SampleID', 'sample_type.samples'])
y = lung_merged_df_2['sample_type.samples']

# Split the data into training and testing sets (e.g., 80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.8, random_state=42)

svm_model = SVC(kernel='linear', C=1.0)

# Train the SVM model
svm_model.fit(X_train, y_train)

# Make predictions on the testing set
y_pred = svm_model.predict(X_test)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
classification_report_result = classification_report(y_test, y_pred)

# Display results
print(f'Accuracy: {accuracy}')
print('Classification Report:')
print(classification_report_result)







og_selected_columns = og.loc[:, og.columns.isin(columns)]
og_selected_columns_t = og_selected_columns.transpose().reset_index()
og_selected_columns_t = og_selected_columns_t.rename(columns={'index': 'SampleID'})
og_merged_df = pd.merge(og_selected_columns_t, selected_rows, how='inner', left_on='SampleID', right_on='submitter_id.samples')
og_merged_df_t = og_merged_df.T
# Exclude the second-to-last row
og_merged_df_2 = og_merged_df_t.iloc[:-2].append(og_merged_df_t.iloc[-1]).transpose()


#Dataset with all genes
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report

# Assuming df is your Pandas DataFrame
# Replace 'TargetColumnName' with the actual name of your target column
X = og_merged_df_2.drop(columns=['SampleID', 'sample_type.samples'])
y = og_merged_df_2['sample_type.samples']

# Split the data into training and testing sets (e.g., 80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.8, random_state=42)

svm_model = SVC(kernel='linear', C=1.0)

# Train the SVM model
svm_model.fit(X_train, y_train)

# Make predictions on the testing set
y_pred = svm_model.predict(X_test)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
classification_report_result = classification_report(y_test, y_pred)

# Display results
print(f'Accuracy: {accuracy}')
print('Classification Report:')
print(classification_report_result)


# In[10]:


#Dataset with all genes


# Assuming df is your Pandas DataFrame
# Replace 'TargetColumnName' with the actual name of your target column
X = merged.drop(columns=['SampleID', 'sample_type.samples'])
y = merged['sample_type.samples']

# Split the data into training and testing sets (e.g., 80% train, 20% test)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42)

svm_model = SVC(kernel='linear', C=1.0)

# Train the SVM model
svm_model.fit(X_train, y_train)

# Make predictions on the testing set
y_pred = svm_model.predict(X_test)

# Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
classification_report_result = classification_report(y_test, y_pred)

# Display results
print(f'Accuracy: {accuracy}')
print('Classification Report:')
print(classification_report_result)

