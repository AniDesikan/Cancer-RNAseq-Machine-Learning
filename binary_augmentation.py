#!/usr/bin/env python
# coding: utf-8

#This code is for functions that help in bootstrapping RNA sequencing data in order to train the machines
# At max the datasets we found off of Xena have about 500 samples, so we want to inject noise into them and bootstrap them.

import pandas as pd
import gzip
import requests
import numpy as np

#Read the csv to augment
with open("colon_tm.csv", 'rt') as file:
    colon_tm = pd.read_csv(file)

#This is a function that injects noise into an entered pandas dataframe, with noise_level being an argument that produces a random seeded number
def inject_feature_noise(data, noise_level):
    if not (0 <= noise_level <= 1):
        raise ValueError("Noise level should be between 0 and 1.")
    #select the second columbn to second to last column and put gene names in a new variable for safekeeping
    new_data = data.iloc[:, 1:]
    new_data = new_data.iloc[:, :-1]
    gene_names = data.iloc[:, 0]
        
#     if not np.issubdtype(new_data.dtype, np.number):
#         raise ValueError("Input data should be of numeric type.")

    # Check if input is a NumPy array or Pandas DataFrame
    if isinstance(new_data, np.ndarray):
        injected_data = new_data.copy()
    elif "pandas" in str(type(data)).lower():
        injected_data = new_data.copy().values
    else:
        raise ValueError("Unsupported data type. Use NumPy array or Pandas DataFrame.")
    # Determine the number of features (columns) in the data
    # Possibly getting rid of this in the future, who knows
    num_features = injected_data.shape[1]
    # Generate random noise for each feature
    # Calculate the row-wise average
    row_avg = np.mean(injected_data, axis=1)
    # Set the noise level to be 5% of the row-wise average(this is the noise_level argument)
    noise_percent = noise_level * row_avg
    # Generate noise with the updated noise level
    noise = np.random.normal(loc=0, scale=noise_percent[:, np.newaxis], size=(len(injected_data), num_features))
    # Inject noise into the data
    injected_data += noise
    injected_data = np.column_stack((gene_names, injected_data))
    return pd.DataFrame(injected_data)

inject_feature_noise(colon_tm, 0.05)

#TODO: Add bootstrapping function and teston colon_tm
