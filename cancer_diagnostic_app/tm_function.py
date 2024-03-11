#!/usr/bin/env python
# coding: utf-8

# The goal of this function will be to take in normalized RNA seq datasets and then return a dataset with only genes that 
# are most important for the tissues that we are looking at.
#################################################
# IMPORTS
#################################################
import pandas as pd
import gzip
import requests
import numpy as np
#################################################
# CONSTANTS
#################################################



#################################################
# HELPER FUNCTIONS
#################################################

# Code that creates the lists of genes in cancer and the list of enriched genes in each tissue to extract for testing
# Make a list of file paths for each sequencing data
# Specify the file path

def sequencing_dataframe(sequencing):
    with open(sequencing, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
        df = pd.read_csv(file, sep='\t')
    return df

def enhanced_gene_dataframe(genes):
    with open(genes, 'rt') as file:
    # Read the TSV file into a pandas DataFrame
        df= pd.read_csv(file)
    df.columns.values[0] = 'Gene'
    df_gene_names = df['Gene']
    df_gene_list = df_gene_names.tolist()
    return df, df_gene_names, df_gene_list

def filter_genes(df, value):
    df = df.sort_values(by='adj.P.Val', ascending=True)
    number = 2 ** value
    df = df.iloc[:number, :]
    df_gene_names = df['Gene']
    return df_gene_names
    
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

def separate_cancer(phenotypes, df):
    if type(df) == str:
        data = sequencing_dataframe(df)
    else:
        data = df
    phenotypes = sequencing_dataframe(phenotypes)
    df = data.T
    df['id'] = df.index
    ph = phenotypes[['submitter_id.samples', 'sample_type.samples']]
    merged_df = ph.merge(df, left_on='submitter_id.samples', right_on='id', how='inner')
    primary_tumor_ids = merged_df.loc[merged_df['sample_type.samples'] == 'Primary Tumor', 'submitter_id.samples']
    cancer = merged_df[merged_df['submitter_id.samples'].isin(primary_tumor_ids)]
    healthy = merged_df[~merged_df['submitter_id.samples'].isin(primary_tumor_ids)]
    cancer['is_cancer'] = 1
    healthy['is_cancer'] = 0
    return cancer, healthy

#################################################
# TM FUNCTION
#################################################

#TM function that returns the dataframe 
# Inputs are the path to the dataframe, the enhanced genes, and the value of how many genes to return
def TM(df_path, enhanced_genes, value):

    #new stuff
    data = sequencing_dataframe(df_path)
    genes, gene_names, gene_list = enhanced_gene_dataframe(enhanced_genes)
    filtered = filter_genes(genes, value)
    ensembl_ids = convert_gene_names_to_ensembl(filtered)

    genes_in_df = data["Ensembl_ID"].tolist()
    genes_in_df_no_decimal = []
    for gene in genes_in_df:
        genes_in_df_no_decimal.append(gene.split('.')[0])
    genes_in_df_no_decimal
    data["Ensembl_ID"] = genes_in_df_no_decimal
    data = data[data["Ensembl_ID"].isin(ensembl_ids)]
    return data

#################################################
# NOISE INJECTION FUNCTION
#################################################

# This gives the random addition of Gaussian noise with variance of 0.1 on normalized data (unless another noise_level number is given)
# Found it on this paper on GAN Data Augmentation
# https://academic.oup.com/bioinformatics/article/39/Supplement_1/i111/7210506

def add_gaussian_noise(df, noise_level=0.1):
    noisy_df = df.copy()
    # Select numeric columns to add noise
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    # Add Gaussian noise to numeric columns
    noisy_df[numeric_cols] += np.random.normal(scale=noise_level, size=noisy_df[numeric_cols].shape)
    return noisy_df

# Different type of noise models suggestion from Hakan:

def add_uniform_noise(df, noise_level=0.1):
    noisy_df = df.copy()
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    noisy_df[numeric_cols] += np.random.uniform(-noise_level, noise_level, size=noisy_df[numeric_cols].shape)
    return noisy_df
    
def add_salt_and_pepper_noise(df, noise_level=0.1):
    noisy_df = df.copy()
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    n_salt = int(noise_level * noisy_df[numeric_cols].size * 0.5)
    for col in numeric_cols:
        indices = np.random.randint(0, len(noisy_df), n_salt)
        noisy_df[col].iloc[indices[:n_salt//2]] = noisy_df[col].max()
        noisy_df[col].iloc[indices[n_salt//2:]] = noisy_df[col].min()
    return noisy_df
    
def add_poisson_noise(df, noise_level=1):
    noisy_df = df.copy()
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    lambda_ = noise_level
    noisy_df[numeric_cols] += np.random.poisson(lambda_, size=noisy_df[numeric_cols].shape)
    return noisy_df

def add_exponential_noise(df, noise_level=1):
    noisy_df = df.copy()
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    noisy_df[numeric_cols] += np.random.exponential(scale=noise_level, size=noisy_df[numeric_cols].shape)
    return noisy_df

def noise_injection(type, df, noise_level):
    type = type.lower()
    if noise_level == None and type == "uniform" or type == "salt and pepper":
        noise_level = 0.1
    elif noise_level == None and type == "poisson" or type == "exponential":
        noise_level = 1
    if type == "uniform":
        noisy_df = add_uniform_noise(df, noise_level)
    if type == "salt and pepper":
        noisy_df = add_salt_and_pepper_noise(df, noise_level)
    if type == "poisson":
        noisy_df = add_poisson_noise(df, noise_level)
    if type == "exponential":
        noisy_df = add_exponential_noise(df, noise_level)
    return noisy_df


#################################################
# INTERPOLATION FUNCTION(s)
#################################################

# There are 2 kinds of interpolation: The more complex "spline interpolation"
# and the simpler, less computationally expensive "linear interpolation"

# Linear interpolation is what the GAN paper above used, I can't find any papers that use interpolation on spline interpolation

# The paper also did a really weird kind of interpolation that I don't really understand, so I've also got another function that does the same
# thing but with a more traditional linear interpolation formula

# Link to code (data_augmentation function): https://github.com/KBRI-Neuroinformatics/WGAN-for-RNASeq-analysis/blob/master/WGAN-for-RNASeq-analysis/preprocess.py

##################################################

def linear_interpolation(data):
    interpolated_samples = []
    data = data.reset_index(drop=True)
    for i in range(len(data.index)):
        try:
            A = data.iloc[i-1,1:]
            B = data.iloc[i,1:]
            interpolated_row = []
            for gene1, gene2 in zip(A, B):
                interpolated_gene = (gene1 + gene2) / 2
                interpolated_row.append(interpolated_gene)
            interpolated_samples.append(interpolated_row)
        except:
            print("failed at index "+ str(i))
    
    # Add index to interpolated samples
    interpolated_df = pd.DataFrame(interpolated_samples, columns=data.columns[1:])
    data = data.iloc[:,1:]
    
    # Concatenate interpolated_df with the original data
    concatenated_df = pd.concat([data, interpolated_df])

    return concatenated_df

