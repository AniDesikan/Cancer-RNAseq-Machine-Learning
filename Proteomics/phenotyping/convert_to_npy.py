import os
import scipy.io
import numpy as np

def convert_mat_to_npy(directory):
    mass_directory = "/Users/ani/Downloads/LungData/LUAD_IMC_Segmentation"
    print(directory)
    # Iterate over all files in the given directory
    for file_name in os.listdir(directory):
        # Check if the file is a .mat file
        if file_name.endswith('.mat'):
            mat_file_path = os.path.join(directory, file_name)
            npy_file_path = os.path.join(directory, file_name.replace('.mat', '.npy'))
            
            # Load the .mat file
            mat_data = scipy.io.loadmat(mat_file_path)
            
            # Save the loaded data as a .npy file
            np.save(npy_file_path, mat_data)
            print(f'Converted {mat_file_path} to {npy_file_path}')

# Specify the directory containing the .mat files
mass_directory = "/Users/ani/Downloads/LungData/LUAD_IMC_Segmentation"
for directory in os.listdir(mass_directory):
    if directory != ".DS_Store":
        directory = "/Users/ani/Downloads/LungData/LUAD_IMC_Segmentation/" + directory
        print(directory)
        convert_mat_to_npy(directory)