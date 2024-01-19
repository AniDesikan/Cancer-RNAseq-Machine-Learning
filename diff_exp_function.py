#!/usr/bin/env python
# coding: utf-8

# In[1]:


# The goal of this function will be to take in normalized RNA seq datasets and then return a dataset with only genes that 
# are most important for the tissues that we are looking at.

import pandas as pd
import gzip
import requests


# In[ ]:


sequencing = "TCGA-COAD.htseq_fpkm.tsv.gz"
enhanced_genes = "tissue_specificity_rna_colon_Group.tab"


# In[2]:


# Code that creates the lists of genes in cancer and the list of enriched genes in each tissue to extract for testing
# Make a list of file paths for each sequencing data
# Specify the file path

with gzip.open(sequencing, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(file, sep='\t')
#Make a list of file_paths for each tissue

# Open the gzip-compressed file and read it with pandas
with open(enhanced_genes, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
    prostate_genes = pd.read_csv(file, sep='\t')
gene_names = prostate_genes['Gene']
gene_list = gene_names.tolist()


# In[ ]:


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
# Convert gene names to Ensembl IDs
ensembl_ids = convert_gene_names_to_ensembl(gene_names)


# In[ ]:


#TM function that returns the dataframe 
def diff_exp(df, ensembl_ids):
    genes_in_df = df["Ensembl_ID"].tolist()
    genes_in_df_no_decimal = []
    for gene in genes_in_df:
        genes_in_df_no_decimal.append(gene.split('.')[0])
    genes_in_df_no_decimal
    df["genes"] = genes_in_df_no_decimal
    df = df[df["genes"].isin(ensembl_ids)]
    return df


# In[ ]:


diff_exp(df, ensembl_ids)

