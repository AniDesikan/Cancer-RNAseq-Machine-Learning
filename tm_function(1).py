#!/usr/bin/env python
# coding: utf-8

# In[12]:


# The goal of this function will be to take in normalized RNA seq datasets and then return a dataset with only genes that 
# are most important for the tissues that we are looking at.

import pandas as pd
import gzip
import requests


# In[32]:


sequencing = "TCGA-COAD.htseq_fpkm.tsv.gz"
lung_sequencing = "TCGA-LUAD.htseq_fpkm.tsv.gz"
colon_enhanced_genes = "tissue_specificity_rna_colon_Group.tab"
lung_enhanced_genes = "tissue_category_rna_lung_Tissue.tsv"


# In[33]:


with gzip.open(lung_sequencing, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
    lung = pd.read_csv(file, sep='\t')
#Make a list of file_paths for each tissue


# In[15]:


# Code that creates the lists of genes in cancer and the list of enriched genes in each tissue to extract for testing
# Make a list of file paths for each sequencing data
# Specify the file path

with gzip.open(sequencing, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(file, sep='\t')
#Make a list of file_paths for each tissue

# Open the gzip-compressed file and read it with pandas
with open(colon_enhanced_genes, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
    colon_genes = pd.read_csv(file, sep='\t')
colon_gene_names = colon_genes['Gene']
colon_gene_list = colon_gene_names.tolist()

with open(lung_enhanced_genes, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
    lung_genes = pd.read_csv(file, sep='\t')
lung_gene_names = lung_genes['Gene']
lung_gene_list = lung_gene_names.tolist()


# In[26]:


# Function to convert gene names in enriched genes to ensembl
def convert_gene_names_to_ensembl(gene_names):
    ensembl_ids = []
    # Ensembl REST API endpoint for mapping gene names to Ensembl IDs
    api_endpoint = "https://rest.ensembl.org/xrefs/symbol/homo_sapiens/"
    for gene_name in gene_names:
        # Make a request to the Ensembl API
        response = requests.get(f"{api_endpoint}{gene_name}?content-type=application/json")
        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Parse the JSON response and extract the Ensembl ID
            data = response.json()
            if data:
                ensembl_id = data[0]["id"]
                ensembl_ids.append(ensembl_id)
            else:
                ensembl_ids.append(None)
        else:
            # If the request was not successful, append None to the list
            ensembl_ids.append(None)
    return ensembl_ids
ensembl_ids = convert_gene_names_to_ensembl(gene_names)


# In[24]:


# Convert gene names to Ensembl IDs
# colon_ensembl_ids = convert_gene_names_to_ensembl(colon_gene_names)
# lung_ensembl_ids = convert_gene_names_to_ensembl(lung_gene_names)
gene_names = pd.concat([colon_gene_names, lung_gene_names], ignore_index=True)
gene_names


# In[27]:


#TM function that returns the dataframe 
def TM(df, ensembl_ids):
    genes_in_df = df["Ensembl_ID"].tolist()
    genes_in_df_no_decimal = []
    for gene in genes_in_df:
        genes_in_df_no_decimal.append(gene.split('.')[0])
    genes_in_df_no_decimal
    df["genes"] = genes_in_df_no_decimal
    df = df[df["genes"].isin(ensembl_ids)]
    return df


# In[29]:


df = TM(df, ensembl_ids)
df


# In[34]:


lung = TM(lung, ensembl_ids)
lung


# In[30]:


csv_file_path = "/home/ani/ML_Project_2024/tm_function/colon_tm.csv"
df.to_csv(csv_file_path, index=False)


# In[35]:


csv_file_path = "/home/ani/ML_Project_2024/tm_function/lung_tm.csv"
lung.to_csv(csv_file_path, index=False)

