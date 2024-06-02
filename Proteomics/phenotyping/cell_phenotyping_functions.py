############
# IMPORTS

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import ndimage
from skimage import io, filters, morphology, segmentation
from skimage.util import img_as_float64
import cv2 
import glob 
import pathlib as Path
from scipy.ndimage import label
############
# Main Code
# Goal is to phenotype cells based on Image mass cytometry information

#######################
# Function: apply_mask_segmentation_files_weights
# PURPOSE: Makes the mask for cell segmentation with the given .h5 weights
# INPUTS:
        # filePath, which I think has a part that has images and text? 
        # metalJMFtablePath is file path to an array of the metals and presences
        # outcomeJMFPath a path that gets made into a directory
        # additional_dilation_table is a optional table with metals and dilations to dilate the final mask image
# OUTPUTS:
        # Makes a directory with the name of outcomeJMFPath variable
        # writes the result and finalMask images to the new directory 

def apply_mask_segmentation_files_weights(filePath, metalJMFtablePath, outcomeJMFPath, additional_dilation_table=None):
    # Read the metal list table as a pandas data frame
    metalListTable = pd.read_csv(metalJMFtablePath)
    lookUpTable = metalListTable["Metal"]
    topKS = metalListTable["Presence"]
    
    # Read the image file
    image = io.imread(filePath)
    
    # Get the file name from the file path
    fileName = filePath.split("/")[-1]
    
    # Create the output directory
    outcomeDir = outcomeJMFPath + "/" + fileName
    os.makedirs(outcomeDir, exist_ok=True)
    
    # Get the text data from the image
    textData = image.textdata
    
    # Loop through each metal in the lookup table
    for lt in range(len(lookUpTable)):
        metal_query = lookUpTable.iloc[lt]
        ind_itr = np.where(textData == metal_query)[0]
        
        # Read the image channel
        result = image[:, :, ind_itr]
        
        # Convert the image to double
        resultIm = img_as_float64(result)
        
        # Compute the mask with weights
        finalMask = compute_mask_with_weights(resultIm, topKS.iloc[lt])
        
        # Apply additional dilation if specified
        if additional_dilation_table is not None:
            metal_dil = additional_dilation_table["Metal"]
            dilation_level = additional_dilation_table["Dilation"]
            if metal_query in metal_dil:
                finalMask = morphology.dilate(finalMask, morphology.disk(dilation_level[metal_dil == metal_query]))
        
        # Save the result and final mask images
        io.imsave(outcomeDir + "/" + metal_query + "_1.png", result)
        io.imsave(outcomeDir + "/" + metal_query + "_2.png", finalMask)
    
##################
# function compute_mask_with_weights
# PURPOSE: Compute the mask with weights using k-means segmentation and morphological operations.
# INPUTS:
#     result: The input image
#     topK: The number of top clusters to consider
# OUTPUTS:
#     finalMask: The computed mask with weights
def compute_mask_with_weights(image, top_clusters):

    filter_size = 3
    filtered_image = filters.median(image, np.ones((filter_size, filter_size)))
    
    initial_clusters = 5
    
    while True:
        cluster_labels, _ = segmentation.kmeans(filtered_image, initial_clusters)
        cluster_intensities = np.zeros(initial_clusters)
        
        for i in range(initial_clusters):
            cluster_intensities[i] = np.mean(filtered_image[cluster_labels == i])
        
        sorted_cluster_indices = np.argsort(cluster_intensities)[::-1]
        cluster_mask = np.zeros_like(cluster_labels)
        for i in range(top_clusters):
            cluster_mask[cluster_labels == sorted_cluster_indices[i]] = 1
        
        cluster_mask = morphology.remove_small_objects(cluster_mask, 20)
        # cluster_mask = morphology.dilate(cluster_mask, morphology.disk(2))  # Uncomment if needed
        mean_mask_value = np.mean(cluster_mask)
        
        if np.abs(1 - mean_mask_value) > 0.1:
            break
        else:
            if initial_clusters >= 2 and top_clusters >= 1:
                initial_clusters -= 1
                top_clusters -= 1
            else:
                break
    
    final_mask = cluster_mask
    return final_mask


################
# function compute_area_of_a_cell_from_indices
# PURPOSE: Compute the area of a cell from its indices.
# INPUTS:
#     curB: A list of indices representing the cell boundary.
#     sizeImage: The size of the image.
# OUTPUTS:
#     regInd: The indices of the cell area.
    
def compute_area_of_a_cell_from_indices(curB, sizeImage):

    # Unravel the indices into separate arrays for X and Y coordinates
    X, Y = np.unravel_index(curB, sizeImage)
    
    # Find the minimum and maximum X and Y coordinates
    minX, minY = np.min(X), np.min(Y)
    maxX, maxY = np.max(X), np.max(Y)
    
    # Shift the X and Y coordinates so that the minimum is at (1,1)
    X = X - minX + 1
    Y = Y - minY + 1
    
    # Create a 2D array of zeros with the size of the bounding box
    markedBoundary = np.zeros((maxX - minX + 1, maxY - minY + 1))
    
    # Shift the X and Y coordinates by one element to the right (wraps around to the start)
    X2 = np.roll(X, -1)
    Y2 = np.roll(Y, -1)
    
    # Interpolate new points between the existing points to create a smoother boundary
    X = np.concatenate((X, (X + 4 * X2) / 5, (2 * X + 3 * X2) / 5, (3 * X + 2 * X2) / 5, (4 * X + X2) / 5, (X + X2) / 2))
    Y = np.concatenate((Y, (Y + 4 * Y2) / 5, (2 * Y + 3 * Y2) / 5, (3 * Y + 2 * Y2) / 5, (4 * Y + Y2) / 5, (Y + Y2) / 2))
    
    # Ensure the coordinates are within the bounds of the image
    X = np.maximum(np.minimum(np.floor(X).astype(int), markedBoundary.shape[0] - 1), 0)
    Y = np.maximum(np.minimum(np.floor(Y).astype(int), markedBoundary.shape[1] - 1), 0)
    
    # Convert the coordinates to a linear index
    cur_nuc_inds = np.ravel_multi_index((X, Y), markedBoundary.shape)
    
    # Set the points in the markedBoundary array to 1
    markedBoundary.flat[cur_nuc_inds] = 1
    
    # Fill in any holes in the markedBoundary array
    markedBoundary = ndimage.binary_fill_holes(markedBoundary)
    
    # Find the coordinates of the points in the markedBoundary array
    regX, regY = np.where(markedBoundary)
    
    # Shift the coordinates back to their original position
    regX = np.maximum(np.minimum(regX + minX, sizeImage[0] - 1), 0)
    regY = np.maximum(np.minimum(regY + minY, sizeImage[1] - 1), 0)
    
    # Convert the coordinates back to a linear index
    regInd = np.ravel_multi_index((regX, regY), sizeImage)
    
    return regInd

###################
# function compute_presence_matrix
# PURPOSE: Compute the presence matrix.
# INPUTS:
#     outcomeJMFPath: The path to the outcome JMF folder.
#     metalList: A list of metal names.
#     scanName: The name of the scan.
#     regIndPath: The path to the region indices file.
# OUTPUTS:
#     presence_matrix: A matrix where each row represents a cell and each column represents a metal, with values indicating the presence of the metal in the cell.

def compute_presence_matrix(outcomeJMFPath, metalList, scanName, regIndPath):

    # Construct the path to the JMF folder
    JMFPath = os.path.join(outcomeJMFPath, scanName)
    
    # Get the number of indices and the number of cells
    nOfInds = len(metalList)
    allRegInds = np.load(regIndPath)
    nOfCells = len(allRegInds)
    
    # Initialize the presence matrix with zeros
    presence_matrix = np.zeros((nOfCells, nOfInds))
    
    # Set the alpha value
    alpha = 0.9
    
    # Loop over each metal
    for j in range(nOfInds):
        # Read the mask image
        curMask = io.imread(os.path.join(JMFPath, f'{metalList[j]}_2.png'))
        
        # If the metal is Yb176Di, dilate the mask
        if metalList[j] == 'Yb176Di':
            curMask = morphology.dilation(curMask, morphology.disk(2))
        
        # Read the image
        curImage = io.imread(os.path.join(JMFPath, f'{metalList[j]}_1.png'))
        
        # Convert the image to double precision and normalize it
        curImage = curImage.astype(np.double)
        if np.any(curMask == 1):
            curImage = curImage / np.mean(curImage[curMask == 1])
        
        # Loop over each cell
        for i in range(nOfCells):
            # Calculate the presence of the metal in the cell
            presence_matrix[i, j] = np.sum(curMask[allRegInds[i]]) / len(allRegInds[i])
    
    return presence_matrix

####################
# function compute_presence_matrix_with_overlay
# PURPOSE: Compute the presence matrix with overlay.
# INPUTS:
#     outcomeJMFPath: The path to the outcome JMF folder.
#     metalList: A list of metal names.
#     scanName: The name of the scan.
#     regIndPath: The path to the region indices file.
#     nucleiPath: The path to the nuclei data file.
# OUTPUTS:
#     Figures of the mask image and raw image with cel boundaries

def compute_presence_matrix_with_overlay(outcomeJMFPath, metalList, scanName, regIndPath, nucleiPath):

    # Construct the path to the JMF folder
    jmf_folder_path = os.path.join(outcomeJMFPath, scanName)
    
    # Get the number of indices and the number of cells
    number_of_metals = len(metalList)
    all_region_indices = np.load(regIndPath)
    nuclei_data = np.load(nucleiPath)
    cell_boundaries = nuclei_data['Boundaries']
    image_size = nuclei_data['nucleiImage'].shape
    number_of_cells = len(all_region_indices)
    
    # Loop over each metal
    for metal_index in range(number_of_metals):
        # Read the mask image
        current_mask = io.imread(os.path.join(jmf_folder_path, f'{metalList[metal_index]}_2.png'))
        
        # Read the image
        current_image = io.imread(os.path.join(jmf_folder_path, f'{metalList[metal_index]}_1.png'))
        
        # Construct the paths to the cell mask and image files
        cell_mask_file_path = os.path.join(jmf_folder_path, f'{metalList[metal_index]}_cells_2.png')
        cell_image_file_path = os.path.join(jmf_folder_path, f'{metalList[metal_index]}_cells_1.png')
        
        # Convert the image to double precision and normalize it
        current_image = current_image.astype(np.double)
        if np.any(current_mask == 1):
            current_image = current_image / np.mean(current_image[current_mask == 1])
        
        # Display the mask image with the cell boundaries
        plt.imshow(current_mask)
        plt.hold(True)
        for cell_index in range(number_of_cells):
            current_boundary = cell_boundaries[cell_index]
            X, Y = np.unravel_index(current_boundary, image_size)
            plt.plot(Y, X, 'r.')
        plt.savefig(cell_mask_file_path)
        plt.close()
        
        # Display the image with the cell boundaries
        plt.imshow(current_image)
        plt.hold(True)
        for cell_index in range(number_of_cells):
            current_boundary = cell_boundaries[cell_index]
            X, Y = np.unravel_index(current_boundary, image_size)
            plt.plot(Y, X, 'r.')
        plt.savefig(cell_image_file_path)

###################
# function create_cell_intensity_indicator
# PURPOSE: Create cell intensity indicators.
# INPUTS:
#     fileNamesForThisGroup: A list of file names for this group.
#     indexFile: The index of the file.
#     pano_files: A list of panorama files.
#     regIndices_files: A list of region indices files.
#     nuc_csv_files: A list of nucleus CSV files.
#     ring_csv_files: A list of ring CSV files.
#     nuc_csv_source_folder: The path to the nucleus CSV source folder.
#     ring_csv_source_folder: The path to the ring CSV source folder.
#     IndicatorsList: A list of indicators.
#     ring_rad: The ring radius.
# OUTPUTS:
#     nuc_pathAddressArray: An array of nucleus path addresses.
#     reg_pathAddressArray: An array of region path addresses.

def create_cell_intensity_indicator(file_names_for_this_group, index_file, pano_files, reg_indices_files, nuc_csv_files, ring_csv_files, nuc_csv_source_folder, ring_csv_source_folder, indicators_list, ring_rad):
    nuc_path_address_array = []
    reg_path_address_array = []

    if not isinstance(file_names_for_this_group[index_file], str):
        for j in range(len(pano_files)):
            full_file_path = pano_files[j]
            scan_name = os.path.splitext(os.path.basename(full_file_path))[0]
            if scan_name == file_names_for_this_group[index_file]:
                pano_file_path = full_file_path
    else:
        cur_cell = file_names_for_this_group[index_file]
        pano_file_path = []
        for tt in range(len(cur_cell)):
            query_str = cur_cell[tt]
            for j in range(len(pano_files)):
                full_file_path = pano_files[j]
                scan_name = os.path.splitext(os.path.basename(full_file_path))[0]
                if scan_name == query_str:
                    pano_file_path.append(full_file_path)

    # find all region indices file addresses
    if not isinstance(file_names_for_this_group[index_file], str):
        for j in range(len(reg_indices_files)):
            full_file_path = reg_indices_files[j]
            folder_path = os.path.dirname(full_file_path)
            endout = re.split('/', folder_path)
            scan_name = endout[-1]
            if scan_name == file_names_for_this_group[index_file]:
                all_reg_inds_path = full_file_path
    else:
        cur_cell = file_names_for_this_group[index_file]
        final_scan_name = ''
        all_reg_inds_path = []
        for tt in range(len(cur_cell)):
            query_str = cur_cell[tt]
            for j in range(len(reg_indices_files)):
                full_file_path = reg_indices_files[j]
                folder_path = os.path.dirname(full_file_path)
                endout = re.split('/', folder_path)
                scan_name = endout[-1]
                if scan_name == query_str:
                    all_reg_inds_path.append(full_file_path)
            cur_scan_name = os.path.splitext(os.path.basename(os.path.dirname(all_reg_inds_path[tt])))[0]
            if tt > 0:
                if len(cur_scan_name) < len(final_scan_name):
                    final_scan_name = cur_scan_name
            else:
                final_scan_name = cur_scan_name

    # find csv (.csv) file addresses
    combined_flag = 0
    nuc_csv_file_path_source_found = 0
    if not isinstance(file_names_for_this_group[index_file], str):
        for j in range(len(nuc_csv_files)):
            full_file_path = nuc_csv_files[j]
            scan_name = os.path.splitext(os.path.basename(full_file_path))[0]
            if scan_name == file_names_for_this_group[index_file]:
                nuc_csv_file_path_source_found = 1
    else:
        combined_flag = 1
        cur_cell = file_names_for_this_group[index_file]
        for tt in range(len(cur_cell)):
            query_str = cur_cell[tt]
            for j in range(len(nuc_csv_files)):
                full_file_path = nuc_csv_files[j]
                scan_name = os.path.splitext(os.path.basename(full_file_path))[0]
                if scan_name == query_str:
                    nuc_csv_file_path_source_found = 1

    ring_csv_file_path_source_found = 0
    if not isinstance(file_names_for_this_group[index_file], str):
        for j in range(len(ring_csv_files)):
            full_file_path = ring_csv_files[j]
            scan_name = os.path.splitext(os.path.basename(full_file_path))[0]
            if scan_name == file_names_for_this_group[index_file]:
                ring_csv_file_path_source_found = 1
    else:
        cur_cell = file_names_for_this_group[index_file]
        for tt in range(len(cur_cell)):
            query_str = cur_cell[tt]
            for j in range(len(ring_csv_files)):
                full_file_path = ring_csv_files[j]
                scan_name = os.path.splitext(os.path.basename(full_file_path))[0]
                if scan_name == query_str:
                    ring_csv_file_path_source_found = 1

    # Need to figure out what's going on with this function so we aren't saving the address arrays here
    # Can't find the function anywhere in the code, it's supposed to save as an MAT file but very confusing so ignoring it

    # ring_mode = 0
    # nuc_path_address_array = save_csv_intensity_matrix_mat_filepaths(nuc_csv_source_folder, final_scan_name, nuc_csv_file_path_source_found, combined_flag, pano_file_path, all_reg_inds_path, indicators_list, ring_mode, ring_rad)

    # ring_mode = 1
    # reg_path_address_array = save_csv_intensity_matrix_mat_filepaths(ring_csv_source_folder, final_scan_name, ring_csv_file_path_source_found, combined_flag, pano_file_path, all_reg_inds_path, indicators_list, ring_mode, ring_rad)

####################
# function create_figure_here
# PURPOSE: Create a figure with cell boundaries overlaid on the nuclei image.
# INPUTS:
#     cell_type_file_path: The path to the cell type file.
#     input_cell_type_folder: The path to the input cell type folder.
#     input_cell_segmentation_folder: The path to the input cell segmentation folder.
#     output_figure_folder: The path to the output figure folder.
#     all_labels: A list of all cell labels.
#     type_colors: A list of colors for each cell type.
# OUTPUTS:
#     A figure saved to the output_figure_folder showing the nuclei image with cell boundaries overlaid.

def create_figure_here(cell_type_file_path, input_cell_type_folder, input_cell_segmentation_folder, output_figure_folder, all_labels, type_colors):
    path_structure = cell_type_file_path.replace(input_cell_type_folder, '')
    folder_name, file_name, _ = os.path.basename(path_structure), os.path.basename(cell_type_file_path), ''

    k = folder_name.find('MATData')
    if k == -1 or k + 7 >= len(folder_name):
        folder_name = ''
    else:
        folder_name = folder_name[k + 8:]

    file_path = os.path.join(input_cell_segmentation_folder, folder_name, file_name, 'nuclei_multiscale.npy')
    nuc_data = np.load(file_path, allow_pickle=True)
    cell_types = np.load(cell_type_file_path, allow_pickle=True)
    boundaries = nuc_data['Boundaries']
    nuclei_image = nuc_data['nucleiImage']

    pdf_cell_type_to_save_path = os.path.join(output_figure_folder, 'PDF', folder_name, f'{file_name}.pdf')
    if not os.path.exists(os.path.join(output_figure_folder, 'PDF', folder_name)):
        os.makedirs(os.path.join(output_figure_folder, 'PDF', folder_name))

    plt.figure()
    plt.imshow(np.zeros_like(nuclei_image))
    plt.hold(True)
    n_of_cells = len(cell_types)
    for i in range(n_of_cells):
        cur_b = boundaries[i]
        if len(cur_b) > 3:
            x, y = np.unravel_index(cur_b, nuclei_image.shape)
            type_cell = cell_types[i]

            res_ind = np.where(np.array(all_labels) == type_cell)[0]
            if len(res_ind) == 0:
                res_ind = len(type_colors)
            plt.fill(y, x, color=ListedColormap(type_colors)(res_ind))

    plt.savefig(pdf_cell_type_to_save_path)


#################
# function create_masks
# PURPOSE: Create masks for each metal in the metal list.
# INPUTS:
#    file_path (str): The path to the input file.
#    metal_jmf_table_path (str): The path to the metal JMF table.
#    outcome_jmf_path (str): The path to the outcome JMF folder.
#    metal_list_thresh_vals (list): A list of threshold values for each metal.
# OUTPUTS:  
#    original and masked images for each metal in the list
def create_masks(file_path, metal_jmf_table_path, outcome_jmf_path, metal_list_thresh_vals):

    # Read the metal list table
    metal_list_table = np.genfromtxt(metal_jmf_table_path, dtype=str)
    metal_list_metals = metal_list_table[:, 0]

    # Get the file name and read the input file
    file_name = os.path.basename(file_path)
    pano = io.imread(file_path)
    all_headers = pano.textdata

    # Create the output folder if it doesn't exist
    os.makedirs(os.path.join(outcome_jmf_path, file_name), exist_ok=True)

    # Loop through each metal in the list
    for a_i, metal_query in enumerate(metal_list_metals):
        # Find the index of the current metal in the headers
        ind_itr = np.where([metal_query in header for header in all_headers])[0]

        # Read the current channel and apply median filter
        cur_ind_matrix = io.imread(file_path, channel=ind_itr[0])
        cur_ind_matrix = filters.median(cur_ind_matrix, np.ones((3, 3)))

        # Threshold the image and create a binary mask
        thresh = metal_list_thresh_vals[a_i]
        bw = cur_ind_matrix
        bw[cur_ind_matrix < thresh] = 0
        bw[cur_ind_matrix >= thresh] = 1

        # Remove small objects from the mask
        bw = ~morphology.remove_small_objects(~bw, 5)
        bw = morphology.remove_small_objects(bw, 5)

        # Save the original and masked images
        orig_ipath = os.path.join(outcome_jmf_path, file_name, f'{metal_query}_1.png')
        mask_ipath = os.path.join(outcome_jmf_path, file_name, f'{metal_query}_2.png')
        io.imsave(orig_ipath, io.imread(file_path, channel=ind_itr[0]))
        io.imsave(mask_ipath, bw)


#######################
# Function create_masks_for_indicators
# PURPOSE: Create masks for indicators in a panorama image based on metal list indicators.

# INPUTS:
# - panorama_file_path (string): File path to the panorama image.
# - adjusting_indicators (list): List of adjusting indicators.
# - metal_list_indicator (list): List of metal list indicators.

# OUTPUTS:
# - Original and masked images for each indicator saved in the 'JMFNew' directory.

# Description:
# This code reads a panorama image, extracts text data from the blue channel,
# and creates binary masks for each indicator based on the metal list indicators.
# The masks are then saved as images in the 'JMFNew' directory.
#
# Note that in the Matlab code, bwareaopen is used to remove small objects from the binary mask. 
# In this Python version, I've replaced it with a combination of cv2.bitwise_not and cv2.erode to achieve a similar effect.
# However, the exact behavior may differ slightly due to the different implementations of morphological operations in OpenCV and Matlab.

def create_masks_for_indicators(panorama_file_path, adjusting_indicators, metal_list_indicator):
    # Read the panorama image
    panorama_image = cv2.imread(panorama_file_path)
    
    # Extract the text data from the blue channel (assuming it's stored there)
    all_headers = panorama_image[:, :, 2]
    
    # Iterate over the adjusting indicators
    for adjusting_indicator_index in range(len(adjusting_indicators)):
        current_indicator = adjusting_indicators[adjusting_indicator_index]
        
        # Initialize metal query and threshold query
        metal_query = ''
        threshold_query = 50
        
        # Find the corresponding metal query and threshold query in the metal list indicator
        # I have no idea what lists they were referencing here, it doesn't even show up in their code. Commenting this entire section out
        # for metal_list_index in range(len(metal_list_indicator)):
        #     if metal_list_indicator[metal_list_index] == current_indicator:
        #         metal_query = metal_list_metals[metal_list_index]
        #         threshold_query = metal_list_thresh[metal_list_index]
        #         break
        
        # Find the indices of the current indicator in the headers
        indicator_indices = np.where(np.char.find(all_headers, current_indicator + '(') != -1)
        
        # Extract the current indicator matrix from the panorama image
        current_indicator_matrix = panorama_image[indicator_indices[0], :, :]
        
        # Load the metal column vector from file
        metal_column_vector = np.load(os.path.join('MetalColumnVectors', metal_query + '.npy'))
        
        # Calculate the histogram and threshold
        histogram, bins = np.histogram(metal_column_vector, bins=100)
        num_bins = 100
        threshold = np.exp(bins[num_bins - threshold_query])
        
        # Create a binary mask for the current indicator
        binary_mask = current_indicator_matrix.copy()
        binary_mask[binary_mask < threshold] = 0
        binary_mask[binary_mask >= threshold] = 255
        
        # Apply morphological operations (different from Matlab's bwareaopen)
        binary_mask = cv2.bitwise_not(binary_mask)
        binary_mask = cv2.bitwise_not(cv2.erode(binary_mask, np.ones((5, 5), np.uint8)))
        binary_mask = cv2.erode(binary_mask, np.ones((5, 5), np.uint8))
        
        # Save the original and masked images
        file_name = os.path.basename(panorama_file_path)
        os.makedirs(os.path.join('JMFNew', file_name), exist_ok=True)
        original_image_path = os.path.join('JMFNew', file_name, f'{current_indicator}_1.png')
        mask_image_path = os.path.join('JMFNew', file_name, f'{current_indicator}_2.png')
        cv2.imwrite(original_image_path, current_indicator_matrix)
        cv2.imwrite(mask_image_path, binary_mask)



############################
# PURPOSE: Dilate and erode priority masks based on metal list indicators and priority tables.

# INPUTS:
# - filePath (string): File path to the image.
# - metalJMFtablePath (string): File path to the metal JMF table.
# - outcomeJMFPath (string): File path to the outcome JMF directory.
# - additionalDilationTable (table): Table containing additional dilation information.
# - priorityTable (table): Table containing priority information.
# - groupName (string): Name of the group to process.

# OUTPUTS:
# - Dilated and eroded masks saved in the outcome JMF directory.

# Description:
# This code reads the metal JMF table, priority table, and additional dilation table,
# and applies dilation and erosion operations to the priority masks based on the metal list indicators.
# The resulting masks are saved in the outcome JMF directory.


def dilate_erode_priority_masks(filePath, metalJMFtablePath, outcomeJMFPath, additionalDilationTable, priorityTable, groupName):
    # Read the metal JMF table
    metalListTable = np.genfromtxt(metalJMFtablePath, dtype=str)
    metalIndicators = metalListTable[:, 0]
    
    # Read the priority table
    lowPriority = priorityTable['Low_Priority']
    highPriority = priorityTable['High_Priority']
    priorityKind = priorityTable['Kind']
    priorityGroup = priorityTable['Group']
    
    # Get the file name and create the output directory
    fileName = os.path.basename(filePath)
    outputDir = os.path.join(outcomeJMFPath, fileName)
    os.makedirs(outputDir, exist_ok=True)
    
    # Read the additional dilation table
    dilationMetals = additionalDilationTable['Metal']
    dilationLevels = additionalDilationTable['Dilation']
    
    # Iterate over the metal indicators
    for metalIndex in range(len(metalIndicators)):
        metalIndicator = metalIndicators[metalIndex]
        
        # Read the current mask
        currentMask = io.imread(os.path.join(outputDir, f'{metalIndicator}_2.png'), as_gray=True)
        currentMask = np.where(currentMask > 0, 1, 0)
        
        # Check if additional dilation is required
        if len(dilationMetals) > 0:
            if metalIndicator in dilationMetals:
                dilationIndex = np.where(dilationMetals == metalIndicator)[0][0]
                currentMask = morphology.dilate(currentMask, morphology.disk(dilationLevels[dilationIndex]))
                io.imsave(os.path.join(outputDir, f'{metalIndicator}_2.png'), currentMask)
    
    # Iterate over the priority table
    for priorityIndex in range(len(lowPriority)):
        if priorityKind[priorityIndex] == 'Strong' and priorityGroup[priorityIndex] == groupName:
            # Read the low and high priority masks
            lowMask = io.imread(os.path.join(outputDir, f'{lowPriority[priorityIndex]}_2.png'), as_gray=True)
            highMask = io.imread(os.path.join(outputDir, f'{highPriority[priorityIndex]}_2.png'), as_gray=True)
            lowMask = np.where(lowMask > 0, 1, 0)
            highMask = np.where(highMask > 0, 1, 0)
            
            # Apply dilation and erosion operations
            lowMask = lowMask - morphology.dilate(np.logical_and(lowMask, highMask), morphology.disk(3))
            lowMask = np.where(lowMask < 0, 0, lowMask)
            lowMask = morphology.remove_small_objects(lowMask, 40)
            highMask = morphology.dilate(highMask, morphology.disk(1))
            
            # Save the updated masks
            io.imsave(os.path.join(outputDir, f'{lowPriority[priorityIndex]}_2.png'), lowMask)
            io.imsave(os.path.join(outputDir, f'{highPriority[priorityIndex]}_2.png'), highMask)
        elif priorityKind[priorityIndex] == 'Weak' and priorityGroup[priorityIndex] == groupName:
            # Read the low and high priority masks
            lowMask = io.imread(os.path.join(outputDir, f'{lowPriority[priorityIndex]}_2.png'), as_gray=True)
            highMask = io.imread(os.path.join(outputDir, f'{highPriority[priorityIndex]}_2.png'), as_gray=True)
            lowMask = np.where(lowMask > 0, 1, 0)
            highMask = np.where(highMask > 0, 1, 0)
            
            # Apply dilation and erosion operations
            highMask = morphology.remove_small_objects(highMask, 50)
            lowMask = lowMask - np.logical_and(lowMask, morphology.dilate(highMask, morphology.disk(3)))
            lowMask = np.where(lowMask < 0, 0, lowMask)
            io.imsave(os.path.join(outputDir, f'{lowPriority[priorityIndex]}_2.png'), lowMask)
            io.imsave(os.path.join(outputDir, f'{highPriority[priorityIndex]}_2.png'), highMask)
        elif priorityKind[priorityIndex] == 'SuperWeak' and priorityGroup[priorityIndex] == groupName:
            # Read the low and high priority masks
            lowMask = io.imread(os.path.join(outputDir, f'{lowPriority[priorityIndex]}_2.png'), as_gray=True)
            highMask = io.imread(os.path.join(outputDir, f'{highPriority[priorityIndex]}_2.png'), as_gray=True)
            lowMask = np.where(lowMask > 0, 1, 0)
            highMask = np.where(highMask > 0, 1, 0)
            
            # Apply dilation and erosion operations
            lowMask = lowMask - np.logical_and(lowMask, highMask)
            lowMask = np.where(lowMask < 0, 0, lowMask)
            io.imsave(os.path.join(outputDir, f'{lowPriority[priorityIndex]}_2.png'), lowMask)


##########################
# function findUniqueStrCell
# PURPOSE: Find unique strings in a list.

# INPUTS:
# - A (list): List of strings.

# OUTPUTS:
# - result (list): List of unique strings.

# Description:
# This function iterates over the input list A and adds each unique string to the result list.

def findUniqueStrCell(A):
    # Initialize an empty list to store the unique strings
    result = []
    
    # Initialize a variable to keep track of the last seen string
    prev_type = ''
    
    # Iterate over the input list
    for item in A:
        # Check if the current string is different from the last seen string
        if item != prev_type:
            # Add the current string to the result list
            result.append(item)
            # Update the last seen string
            prev_type = item
    
    # Return the list of unique strings
    return result

##############################
# function findWherePathAddressesAreStr=ored
# PURPOSE: Find the index where a path address is stored.

# INPUTS:
# - allPaths (list): List of paths.
# - regPaths (list or string): Registered paths.

# OUTPUTS:
# - result (int): Index where the path adderss is stored.

# Description:
# This function searches for the index in allPaths where the registered path (regPaths) is stored.

def findWherePathAddressesAreStored(allPaths, regPaths):
    # Check if regPaths is a list
    if isinstance(regPaths, list):
        # Get the first element of the list
        query = regPaths[0]
    else:
        # regPaths is a string
        query = regPaths
    
    # Extract the file name from the query path
    query = os.path.basename(os.path.dirname(query))
    
    # Iterate over allPaths
    for i, curPaths in enumerate(allPaths):
        # Check if curPaths is a list
        if isinstance(curPaths, list):
            # Iterate over curPaths
            for curPath in curPaths:
                # Extract the file name from the current path
                fileName = os.path.basename(os.path.dirname(curPath))
                # Check if the file name matches the query
                if fileName == query:
                    # Return the index
                    return i
        else:
            # curPaths is a string
            # Extract the file name from the current path
            fileName = os.path.basename(os.path.dirname(curPaths))
            # Check if the file name matches the query
            if fileName == query:
                # Return the index
                return i
    
    # If no match is found, return None
    return None

#################
# FUNCTION: find_all_files

# PURPOSE: Find all files for a given group.

# INPUTS:
# - indexGroup (int): Index of the group.
# - nuc_csv_source_folder (str): Folder path for nuclide CSV files.
# - ring_csv_source_folder (str): Folder path for ring CSV files.
# - groupPatientIDTableData (list): Table data for patient IDs.
# - imageIDColumn (list): List of image IDs.
# - fileNameColumn (list): List of file names.

# OUTPUTS:
# - fileNamesForThisGroup (list): List of file names for the current group.
# - nuc_csv_files (list): List of nuclide CSV files.
# - ring_csv_files (list): List of ring CSV files.

# Description:
# This function finds all files (nuclide and ring CSV files) for a given group and extracts the corresponding file names.

def find_all_files(indexGroup, nuc_csv_source_folder, ring_csv_source_folder, groupPatientIDTableData, imageIDColumn, fileNameColumn):
    # Initialize lists to store file names and CSV files
    fileNamesForThisGroup = []
    nuc_csv_files = []
    ring_csv_files = []
    
    # Find all nuclide CSV files
    for i in range(100):
        nuc_csv_files.extend(glob.glob(os.path.join(nuc_csv_source_folder, '*') + '*.csv'))
    
    # Find all ring CSV files
    for i in range(100):
        ring_csv_files.extend(glob.glob(os.path.join(ring_csv_source_folder, '*') + '*.csv'))
    
    # Get the current group keys
    cur_group_keys = groupPatientIDTableData[:, indexGroup]
    cur_group_keys = [key for key in cur_group_keys if key != '']
    
    # Initialize a list to store file names for the current group
    fileNamesForThisGroup = [None] * len(cur_group_keys)
    
    # Iterate over the current group keys
    for i, cur_key in enumerate(cur_group_keys):
        # Find the row index for the current key in the image ID column
        rowIndex = [j for j, x in enumerate(imageIDColumn) if cur_key in x]
        
        # If only one row index is found, assign the corresponding file name
        if len(rowIndex) == 1:
            fileNamesForThisGroup[i] = fileNameColumn[rowIndex[0]]
        # If multiple row indices are found, assign a list of file names
        else:
            fileNamesForThisGroup[i] = [fileNameColumn[j] for j in rowIndex]
    
    return fileNamesForThisGroup, nuc_csv_files, ring_csv_files

###################################
# Purpose: Find registration files.

# INPUTS:
# - fileNamesForThisGroup (list): List of file names.
# - indexFile (int): Index of the file.
# - regIndices_files (list): List of registration index files.

# Outputs:
# - allRegIndsPath (list): List of registration file paths.
# - finalScanName (str): Final scan name.

# Description:
# This function finds the registration files corresponding to the file names in fileNamesForThisGroup.
# Registration files contain information about the alignment of multiple images or datasets in imaging mass cytometry. 
# The function returns the paths to the registration files and the final scan name, allowing for the integration of data from multiple images or datasets.

def find_reg_files(fileNamesForThisGroup, indexFile, regIndices_files):
    # Check if fileNamesForThisGroup is a string or a list of strings
    if isinstance(fileNamesForThisGroup[indexFile], str):
        # Iterate over the registration index files
        for reg_file in regIndices_files:
            # Get the full file path and scan name
            fullFilePath = reg_file.name
            folderPath, _, _ = os.path.basename(fullFilePath)
            endout = folderPath.split(os.sep)
            scanName = endout[-1]
            # Check if the scan name matches the file name
            if scanName == fileNamesForThisGroup[indexFile]:
                allRegIndsPath = fullFilePath
                break
        # Extract the final scan name
        _, finalScanName, _ = os.path.basename(os.path.dirname(allRegIndsPath))
    else:
        # fileNamesForThisGroup is a list of strings
        curCell = fileNamesForThisGroup[indexFile]
        finalScanName = ''
        allRegIndsPath = []
        # Iterate over the file names in curCell
        for queryStr in curCell:
            # Iterate over the registration index files
            for reg_file in regIndices_files:
                # Get the full file path and scan name
                fullFilePath = reg_file.name
                folderPath, _, _ = os.path.basename(fullFilePath)
                endout = folderPath.split(os.sep)
                scanName = endout[-1]
                # Check if the scan name matches the query string
                if scanName == queryStr:
                    allRegIndsPath.append(fullFilePath)
                    break
            # Extract the final scan name
            _, cur_scanName, _ = os.path.basename(os.path.dirname(allRegIndsPath[-1]))
            if len(cur_scanName) < len(finalScanName) or not finalScanName:
                finalScanName = cur_scanName
    return allRegIndsPath, finalScanName

##############################
# FUNCTION: fixFillingNuclei

# PURPOSE: Fix filling nuclei.

# INPUTS:
# - inputImage (numpy array): Input image.
# - nucleiImage (numpy array): Nuclei image.

# OUTPUTS:
# - inputImage (numpy array): Output image with filled nuclei.

# Description:
# This function fills nuclei in the input image based on the nuclei image.

def fixFillingNuclei(inputImage, nucleiImage):
    # Label connected regions in the inverse of inputImage
    labeledImage, numRegions = label(~inputImage)
    
    # Iterate over the labeled regions
    for regionIndex in range(1, numRegions+1):
        # Find the pixels in the current region
        regionPixels = np.where(labeledImage == regionIndex)
        
        # Calculate the ratio of covered pixels in the region
        coveredPixelsRatio = np.sum(nucleiImage[regionPixels] != 0) / len(regionPixels[0])
        
        # If the ratio is greater than 0.7, fill the region in inputImage
        if coveredPixelsRatio > 0.7:
            inputImage[regionPixels] = 1
    
    return inputImage

#########################
# FUNCTION: get_the_first_list_indicators

# Purpose: Get the first list indicators.

# Inputs:
# - rule_table_location (str): Location of the rule table file.

# Outputs:
# - firstListIndicators (list): List of first indicators.
# - tableHeaders (list): List of table headers.
# - ruleTable (pandas DataFrame): Rule table data.
# - relationship (list): List of relationships.
# - before_type (list): List of before types.
# - indicators (list): List of indicators.
# - metals (list): List of metals.
# - after_type (list): List of after types.
# - sensity_level (list): List of sensitivity levels.

# Description:
# This function reads the rule table from a file and extracts the first list indicators and other related data.

def get_the_first_list_indicators(rule_table_location):

    pass_no = 1

    # Read the rule table from the file
    ruleTable = pd.read_table(rule_table_location)

    # Get the table headers
    tableHeaders = ruleTable.columns.tolist()

    # Get the first list indicators
    firstListIndicators = ruleTable[ruleTable.columns[ruleTable.columns.str.contains(f'METAL_{pass_no}')]].values.tolist()

    # Remove empty values from the list
    firstListIndicators = [x for x in firstListIndicators if x != '']

    # Get the relationship, before type, indicators, metals, after type, and sensitivity level data
    relationship = ruleTable[ruleTable.columns[ruleTable.columns.str.contains(f'RELATIONSHIP_{pass_no}')]].values.tolist()
    before_type = ruleTable[ruleTable.columns[ruleTable.columns.str.contains(f'TYPE_BEFORE_{pass_no}')]].values.tolist()
    indicators = ruleTable[ruleTable.columns[ruleTable.columns.str.contains(f'IND_{pass_no}')]].values.tolist()
    metals = firstListIndicators
    after_type = ruleTable[ruleTable.columns[ruleTable.columns.str.contains(f'TYPE_AFTER_{pass_no}')]].values.tolist()
    sensity_level = ruleTable[ruleTable.columns[ruleTable.columns.str.contains(f'Sensitivity_Level_{pass_no}')]].values.tolist()
    return firstListIndicators, tableHeaders, ruleTable, relationship, before_type, indicators, metals, after_type, sensity_level

############################################
# FUNCTION: main_cell_type_assignment_v3

# PURPOSE: Assign cell types based on rules.

# INPUTS:
# - rule_table_location (str): Location of the rule table file.
# - outcomeJMFPath (str): Path to the outcome JMF file.
# - scanName (str): Name of the scan.
# - regIndPath (str): Path to the registration indices file.
# - ringIndPath (str): Path to the ring indices file.

# OUTPUTS:
# - allCellTypes (list): List of assigned cell types.

# Description:
# This function assigns cell types based on rules defined in the rule table.

def main_cell_type_assignment_v3(rule_table_location, outcomeJMFPath, scanName, regIndPath, ringIndPath):
    # Get the first list indicators and other data from the rule table
    firstListIndicators, tableHeaders, ruleTable, relationship, before_type, indicators, metals, after_type, sensity_level = get_the_first_list_indicators(rule_table_location)
    
    # Get the number of passes from the table headers
    nOfPasses = len(tableHeaders) // 6
    
    # Import the registration indices
    allRegInds = np.genfromtxt(regIndPath)
    
    # Get the number of cells
    nOfCells = len(allRegInds)
    
    # Initialize the allCellTypes list
    allCellTypes = [[] for _ in range(nOfCells)]
    
    # Iterate over the passes
    for pass_no in range(1, nOfPasses + 1):
        # Run the first pass or subsequent passes
        if pass_no == 1:
            allCellTypes = run_first_pass(outcomeJMFPath, scanName, regIndPath, ringIndPath, firstListIndicators, relationship, indicators, after_type, metals, allCellTypes)
        else:
            allCellTypes = type_assignment_per_pass_v3_priority(tableHeaders, allCellTypes, pass_no, ruleTable, outcomeJMFPath, scanName, regIndPath, ringIndPath)
    
    return allCellTypes

############################################
# FUNCTION: obtain_thresh_vals

# PURPOSE: Obtain threshold values for metals.

# INPUTS:
# - metalJMFtablePath (str): Path to the metal JMF table file.

# OUTPUTS:
# - metalListThreshVals (list): List of threshold values for metals.

# Description:
# This function reads the metal JMF table, extracts the metal names and presence thresholds,
# and calculates the threshold values using the metal column vectors.

def obtain_thresh_vals(metalJMFtablePath):
    # Read the metal JMF table
    metalListTable = pd.read_table(metalJMFtablePath)
    
    # Extract the metal names and presence thresholds
    metalListMetals = metalListTable['Metal']
    metalListThresh = metalListTable['Presence']
    
    # Initialize the list to store threshold values
    metalListThreshVals = [0] * len(metalListMetals)
    
    # Iterate over the metals
    for aI, metal_query in enumerate(metalListMetals):
        # Get the threshold query
        thresh_query = metalListThresh[aI]
        
        # Load the metal column vector
        lR = np.load(os.path.join('MetalColumnVectors', f'{metal_query}.npy'))
        
        # Calculate the histogram bins
        bins = np.histogram(lR, bins=100)[1]
        
        # Calculate the threshold value
        metalListThreshVals[aI] = np.exp(bins[99 - thresh_query])
    
    return metalListThreshVals

############################################
# FUNCTION: prepare_JMFMats

# Purpose: Prepare JMF matrices.

# Inputs:
# - scanName (str): Name of the scan.
# - metalJMFtablePath (str): Path to the metal JMF table file.
# - outcomeJMFPath (str): Path to the outcome JMF file.
# - cIMatOuputFolder (str): Output folder for the cIMat file.
# - regIndPath (str): Path to the registration indices file.
# - ringIndPath (str): Path to the ring indices file.

# Outputs:
# - cIMat (numpy array): Combined presence matrix.

# Description:
# This function prepares the JMF matrices by computing the presence matrices for nuclei and rings,
# and combining them into a single matrix.

def prepare_JMFMats(scanName, metalJMFtablePath, outcomeJMFPath, cIMatOuputFolder, regIndPath, ringIndPath):
    # Read the metal JMF table
    metalJMFtable = pd.read_table(metalJMFtablePath)
    
    # Extract the metal list
    metalList = metalJMFtable['Metal']
    
    # Compute the presence matrices
    nuc_presence_matrix = compute_presence_matrix(outcomeJMFPath, metalList, scanName, regIndPath)
    ring_presence_matrix = compute_presence_matrix(outcomeJMFPath, metalList, scanName, ringIndPath)
    
    # Combine the presence matrices
    cIMat = np.maximum(nuc_presence_matrix, ring_presence_matrix)
    
    # Create the output folder if it doesn't exist
    if not os.path.exists(cIMatOuputFolder):
        os.makedirs(cIMatOuputFolder)
    
    # Save the cIMat file
    np.save(os.path.join(cIMatOuputFolder, f'{scanName}.npy'), cIMat)
    
    return cIMat

############################################
# FUNCTION: prepare_JMFMats_with_overlay

# PURPOSE: Prepare JMF matrices with overlay.

# INPUTS:
# - scanName (str): Name of the scan.
# - metalJMFtablePath (str): Path to the metal JMF table file.
# - outcomeJMFPath (str): Path to the outcome JMF file.
# - regIndPath (str): Path to the registration indices file.
# - ringIndPath (str): Path to the ring indices file.
# - nucleiPath (str): Path to the nuclei file.

# Description:
# This function prepares the JMF matrices with overlay by computing the presence matrices
# with the overlay information.

def prepare_JMFMats_with_overlay(scanName, metalJMFtablePath, outcomeJMFPath, regIndPath, ringIndPath, nucleiPath):
    # Read the metal JMF table
    metalJMFtable = pd.read_table(metalJMFtablePath)
    
    # Extract the metal list
    metalList = metalJMFtable['Metal']
    
    # Compute the presence matrices with overlay
    compute_presence_matrix_with_overlay(outcomeJMFPath, metalList, scanName, regIndPath, ringIndPath, nucleiPath)

############################################
# Function: read_tables_list_files

# Purpose: Reads in various tables and files from a configuration file

# Inputs:
#   config_path: string containing the path to a configuration file

# Outputs:
#   IndicatorsList: list of indicators
#   nuc_csv_source_folder: string containing the path to a folder containing nuclei CSV files
#   ring_csv_source_folder: string containing the path to a folder containing ring CSV files
#   rule_table_location: string containing the path to a file containing a rule table
#   nOfGroups: integer containing the number of groups
#   groupPatientIDTableData: list of group patient ID data
#   headerGroupPatientID: list of header names for the group patient ID table
#   fileNameColumn: list of file names
#   imageIDColumn: list of image IDs
#   regIndices_files: list of region indices files
#   pano_files: list of pano files
#   IndicatorsListPath: string containing the path to a file containing the indicators list table
#   input_pano_folder: string containing the path to a folder containing pano data
#   input_cell_segmentation_folder: string containing the path to a folder containing cell segmentation data
#   output_folder: string containing the path to a folder where output files will be saved

# Description:
#   This function reads in various tables and files from a configuration file and returns the data in a tuple.
# =========================================================

def read_tables_list_files(config_path: str) -> tuple:
    # Compile the configuration file to get the input and output folders
    input_pano_folder, input_cell_segmentation_folder, input_keyPatientID, group_PatientID, output_folder, IndicatorsListPath, nuc_csv_source_folder, ring_csv_source_folder, rule_table_location = compile_config_file(config_path)
    
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Read in the indicators list table
    IndicatorsListTable = pd.read_table(IndicatorsListPath)
    IndicatorsList = IndicatorsListTable.values.tolist()
    
    # Print some messages to the console
    print(f'The code is reading from pano data folder {input_pano_folder}')
    print(f'The code is reading from cell segmentation data folder {input_cell_segmentation_folder}')
    print(f'The code is reading key patient id file from {input_keyPatientID}')
    print(f'The code is reading group patient ids file from {group_PatientID}')
    
    # Read in the group patient ID table
    groupPatientIDTable = pd.read_table(group_PatientID)
    nOfGroups = groupPatientIDTable.shape[1]
    groupPatientIDTableData = groupPatientIDTable.values.tolist()
    headerGroupPatientID = groupPatientIDTable.columns.tolist()
    
    # Read in the key patient ID table
    inputKeyPatientIDTable = pd.read_table(input_keyPatientID)
    headerKeyTable = inputKeyPatientIDTable.columns.tolist()
    columnIndexFileName = [i for i, col in enumerate(headerKeyTable) if 'FileName' in col]
    columnIndexImageID = [i for i, col in enumerate(headerKeyTable) if 'UniqueImageID' in col]
    inputKeyPatientIDData = inputKeyPatientIDTable.values.tolist()
    fileNameColumn = [row[columnIndexFileName[0]] for row in inputKeyPatientIDData]
    imageIDColumn = [row[columnIndexImageID[0]] for row in inputKeyPatientIDData]
    
    # Get the region indices files and pano files
    regIndices_files = []
    for i in range(10):
        regIndices_files.extend(Path(input_cell_segmentation_folder).rglob('allRegionIndices.npy'))
    
    pano_files = []
    for i in range(10):
        pano_files.extend(Path(input_pano_folder).rglob('*.txt'))
    
    # Loop through the pano files and check if a corresponding region indices file exists
    for i in range(len(pano_files)):
        fileNameQ = Path(pano_files[i]).name
        
        foundQuery = 0
        for j in range(len(regIndices_files)):
            fileNameR = Path(regIndices_files[j]).name
            
            if fileNameQ == fileNameR:
                foundQuery = 1
                break
        
        if foundQuery == 0:
            folderPath, fileName, _ = Path(pano_files[i]).parts
            folderStructure = str(Path(folderPath).relative_to(input_pano_folder))
            folderPathNuc = Path(input_cell_segmentation_folder) / folderStructure / fileName
            nucleiFilePath = folderPathNuc / 'nuclei_multiscale.npy'
            data = io.loadmat(nucleiFilePath)
            Boundaries = data['Boundaries']
            # Compute the area of each cell
            I = data['nucleiImage']
            allRegInds = np.empty((len(Boundaries),), dtype=object)
            for k in range(len(Boundaries)):
                curB = Boundaries[k]
                allRegInds[k] = compute_area_of_a_cell_from_indices(curB, I.shape)
            
            # Save the region indices to a file
            io.savemat(folderPathNuc / 'allRegionIndices.npy', {'allRegInds': allRegInds})
            
            # Compute the ring indices
            allRingIndsPath = folderPathNuc / 'allRingIndices.npy'
            nOfCells = len(allRegInds)
            
            allRingInds = np.empty((nOfCells,), dtype=object)
            for iii in range(nOfCells):
                ringImage = np.zeros(I.shape)
                ringImage[allRegInds[iii]] = 1
                ringImageDil = np.binary_dilation(ringImage, np.ones((3, 3)))
                ringImage = ringImageDil - ringImage
                allRingInds[iii] = np.argwhere(ringImage).flatten()
            
            # Save the ring indices to a file
            io.savemat(allRingIndsPath, {'allRingInds': allRingInds})
    
    # Return the outputs
    return IndicatorsList, nuc_csv_source_folder, ring_csv_source_folder, rule_table_location, nOfGroups, groupPatientIDTableData, headerGroupPatientID, fileNameColumn, imageIDColumn, regIndices_files, pano_files, IndicatorsListPath, input_pano_folder, input_cell_segmentation_folder, output_folder

############################################
# FUNCTION: compile_config_file

# PURPOSE: Compile the configuration file into a set of parameters.

# INPUTS:
# - config_file_path (str): Path to the configuration file.

# OUTPUTS:
# - input_pano_folder (str): Input pano folder.
# - input_cell_segmentation_folder (str): Input cell segmentation folder.
# - input_keyPatientID (str): Input key patient ID file.
# - group_PatientID (str): Group patient ID file.
# - output_folder (str): Output folder.
# - IndicatorsListPath (str): Indicators list path.
# - nuc_csv_source_folder (str): Nuclei CSV source folder.
# - ring_csv_source_folder (str): Ring CSV source folder.
# - rule_table_location (list): List of rule table locations.

# Description:
# This function reads the configuration file, extracts the necessary parameters,
# and returns them as a set of outputs.

def compile_config_file(config_file_path):
    # Read the configuration file
    config_table = pd.read_table(config_file_path)
    
    # Extract the parameters and values
    params = config_table['PARAM']
    values = config_table['VALUE']
    
    # Initialize the output variables
    input_pano_folder = ''
    input_cell_segmentation_folder = ''
    input_keyPatientID = ''
    group_PatientID = ''
    output_folder = ''
    IndicatorsListPath = ''
    nuc_csv_source_folder = ''
    ring_csv_source_folder = ''
    rule_table_location = []
    
    # Loop through the parameters and extract the values
    for i, param in enumerate(params):
        if param == 'input_pano_folder':
            input_pano_folder = values[i]
        elif param == 'input_cell_segmentation_folder':
            input_cell_segmentation_folder = values[i]
        elif param == 'input_keyPatientID':
            input_keyPatientID = values[i]
        elif param == 'group_PatientID':
            group_PatientID = values[i]
        elif param == 'output_folder':
            output_folder = values[i]
        elif param == 'IndicatorsListPath':
            IndicatorsListPath = values[i]
        elif param == 'nuc_csv_source_folder':
            nuc_csv_source_folder = values[i]
        elif param == 'ring_csv_source_folder':
            ring_csv_source_folder = values[i]
        elif 'rule_table_location_' in param:
            cohort_name = param.split('_')[1]
            rule_table_location.append((cohort_name, values[i]))
    
    # Return the output variables
    return (input_pano_folder, input_cell_segmentation_folder, input_keyPatientID, 
            group_PatientID, output_folder, IndicatorsListPath, nuc_csv_source_folder, 
            ring_csv_source_folder, rule_table_location)

############################################
# FUNCTION: run_first_pass

# PURPOSE: Run the first pass of the cell type identification algorithm.

# INPUTS:
# - outcomeJMFPath (str): Path to the outcome JMF folder.
# - scanName (str): Name of the scan.
# - regIndPath (str): Path to the registration indices file.
# - ringIndPath (str): Path to the ring indices file.
# - firstListIndicators (list): List of first-level indicators.
# - relationship (list): List of relationships between indicators.
# - indicators (list): List of indicators.
# - after_type (list): List of after-type indicators.
# - metals (list): List of metals.
# - allCellTypes (array): Array of all cell types.

# OUTPUTS:
# - allCellTypes (array): Updated array of all cell types.

# Description:
# This function runs the first pass of the cell type identification algorithm,
# computing the presence matrices and assigning the cell types based on the
# highest indicator values.
# """

def run_first_pass(outcomeJMFPath, scanName, regIndPath, ringIndPath, firstListIndicators, relationship, indicators, after_type, metals, allCellTypes):
    NONE_THRESH = 0.07
    pass_no = 1
    JMFPath = os.path.join(outcomeJMFPath, scanName)
    
    # Load the registration and ring indices
    allRegInds = np.load(regIndPath)
    allRingInds = np.load(ringIndPath)
    
    # Initialize the presence matrices
    nOfCells = len(allRegInds)
    nuc_presence_matrix = np.zeros((nOfCells, len(firstListIndicators)))
    ring_presence_matrix = np.zeros((nOfCells, len(firstListIndicators)))
    
    # Loop through the indicators and compute the presence matrices
    for j, indicator in enumerate(firstListIndicators):
        curMask = io.imread(os.path.join(JMFPath, f'{indicator}_2.png'), as_gray=True)
        for i in range(nOfCells):
            nuc_presence_matrix[i, j] = np.mean(curMask[allRegInds[i]])
            ring_presence_matrix[i, j] = np.mean(curMask[allRingInds[i]])
    
    # Compute the final presence matrix
    presence_matrix = nuc_presence_matrix  # + ring_presence_matrix
    
    # Apply the relationship logic
    rel_exits = 1
    if isinstance(relationship[0], (int, float)):
        rel_exits = 0
    
    arri = []
    for rule_row_itr in range(len(relationship)):
        if len(indicators[rule_row_itr]) > 0:
            if rel_exits and relationship[rule_row_itr] == 'AND':
                rri = rule_row_itr
                cur_presence_matrix = io.imread(os.path.join(JMFPath, f'{firstListIndicators[rri]}_2.png'), as_gray=True)
                arri.append(rri)
                while rri < len(relationship) and relationship[rri] == 'AND':
                    rri += 1
                    if relationship[rri] == 'AND':
                        arri.append(rri)
                        cur_presence_matrix = np.logical_and(cur_presence_matrix, io.imread(os.path.join(JMFPath, f'{firstListIndicators[rri]}_2.png'), as_gray=True))
                io.imsave(os.path.join(JMFPath, f'{after_type[rule_row_itr]}_2.png'), cur_presence_matrix)
                io.imsave(os.path.join(JMFPath, f'{after_type[rule_row_itr]}_1.png'), cur_presence_matrix)
                firstListIndicators.append(after_type[rule_row_itr])
    
    # Remove duplicates from firstListIndicators
    firstListIndicators = list(dict.fromkeys(firstListIndicators))
    
    # Find the highest indicator values and their indices
    highIndicatorVals = np.max(presence_matrix, axis=1)
    highIndicatorIds = np.argmax(presence_matrix, axis=1)
    
    # Assign the cell types based on the highest indicator values
    newCellTypes = []
    for cellItr in range(nOfCells):
        if highIndicatorVals[cellItr] >= NONE_THRESH:
            T = np.where([metal in firstListIndicators[highIndicatorIds[cellItr]] for metal in metals])[0]
            if len(T) > 0:
                newCellTypes.append(after_type[T[0]])
            else:
                newCellTypes.append(firstListIndicators[highIndicatorIds[cellItr]])
        else:
            newCellTypes.append('NONE')
    
    # Update the allCellTypes array
    allCellTypes[:, pass_no] = newCellTypes
    
    return allCellTypes

############################################
# FUNCTION: save_celltypes

# PUPROSE: Save the cell types to a MAT file.

# INPUTS:
# - allCellTypes (array): Array of all cell types.
# - panoFilePath (str): Path to the pano file.
# - input_pano_folder (str): Input pano folder.
# - input_cell_segmentation_folder (str): Input cell segmentation folder. Was part of the old code, but in python removed
# - output_folder (str): Output folder.

# OUTPUTS:
# - cellTypeToSavePath (str): Path to the saved MAT file.

# Description:
# This function saves the cell types to a MAT file in the output folder.
# It extracts the folder name and file name from the pano file path,
# and creates the output folder if it does not exist.

def save_celltypes(allCellTypes, panoFilePath, input_pano_folder, output_folder):
    # Get the last column of allCellTypes
    cellTypes = allCellTypes[:, -1]
    
    # Extract the path structure from the pano file path
    pathStructure = panoFilePath.replace(input_pano_folder, '')
    
    # Split the path structure into folder name, file name, and extension
    folderName, fileName, _ = pathStructure.split('/')
    
    # Create the output file path
    cellTypeToSavePath = os.path.join(output_folder, 'NPYData', folderName, f'{fileName}.npy')
    
    # Create the output folder if it does not exist
    if not os.path.exists(os.path.join(output_folder, 'NPYData', folderName)):
        os.makedirs(os.path.join(output_folder, 'NPYData', folderName))
    
    # Save the cell types to the NPY file
    np.save(cellTypeToSavePath, cellTypes)
    
    # Return the output file path
    return cellTypeToSavePath


############################################
# FUNCTION:type_assignment_per_pass_v3_priority

# PURPOSE: Assign cell types based on rules.

# INPUTS:
# - tableHeaders (list): List of table headers.
# - allCellTypes (array): Array of all cell types.
# - pass_no (int): Pass number.
# - ruleTable (array): Rule table.
# - outcomeJMFPath (str): Path to the outcome JMF folder.
# - scanName (str): Name of the scan.
# - regIndPath (str): Path to the registration indices file.
# - ringIndPath (str): Path to the ring indices file.

# OUTPUTS:
# - allCellTypes (array): Updated array of all cell types.

# Description:
# This function assigns cell types based on rules defined in the rule table.
# It uses the priority system to assign cell types.

def type_assignment_per_pass_v3_priority(tableHeaders, allCellTypes, pass_no, ruleTable, outcomeJMFPath, scanName, regIndPath, ringIndPath):
    BIG_THRESH = 0.5

    # Load the registration and ring indices
    allRegInds = np.load(regIndPath)
    allRingInds = np.load(ringIndPath)

    # Get the JMF path
    JMFPath = os.path.join(outcomeJMFPath, scanName)

    # Initialize the new cell types
    newCellTypes = allCellTypes[:, pass_no - 1]

    # Loop through the rules
    for rule_row_itr in range(len(ruleTable)):
        # Get the indicators, metals, and after type
        before_type = ruleTable[rule_row_itr][next((i for i, x in enumerate(tableHeaders) if f'TYPE_BEFORE_{pass_no}' in x), None)]
        indicators = ruleTable[rule_row_itr][next((i for i, x in enumerate(tableHeaders) if f'IND_{pass_no}' in x), None)]
        metals = ruleTable[rule_row_itr][next((i for i, x in enumerate(tableHeaders) if f'METAL_{pass_no}' in x), None)]
        after_type = ruleTable[rule_row_itr][next((i for i, x in enumerate(tableHeaders) if f'TYPE_AFTER_{pass_no}' in x), None)]
        sensitivity_level = ruleTable[rule_row_itr][next((i for i, x in enumerate(tableHeaders) if f'Sensitivity_Level_{pass_no}' in x), None)]

        # Check if the before type exists
        before_type_exist = 1
        if isinstance(before_type, (int, float)):
            before_type_exist = 0

        # Get the candidate cell IDs
        candidate_cell_ids = []
        if before_type_exist and before_type is not None:
            for i in range(len(newCellTypes)):
                if newCellTypes[i] == before_type:
                    candidate_cell_ids.append(i)
        else:
            for i in range(len(newCellTypes)):
                if newCellTypes[i] == '':
                    candidate_cell_ids.append(i)

        # Get the current rule row indicators
        cur_rule_row_ind = indicators
        sign_pos_neg = cur_rule_row_ind[-1]

        # Check the sign
        if sign_pos_neg == '+':
            sign_pos_neg = 1
        elif sign_pos_neg == '-':
            sign_pos_neg = -1
        else:
            sign_pos_neg = 0

        # Load the metal image
        if metals != '':
            curMask = io.imread(os.path.join(JMFPath, f'{metals}_2.png'), as_gray=True)
        else:
            else_cond = 1

        # Loop through the candidate cell IDs
        for i in range(len(candidate_cell_ids)):
            cI = candidate_cell_ids[i]
            if else_cond == 0:
                score1 = np.mean(curMask[allRegInds[cI]])
                score2 = np.mean(curMask[allRingInds[cI]])
                score_max = score1
                # score_max = max(score1, score2)

            # Check the sign and assign the cell type
            if sign_pos_neg == 1:
                if score_max >= sensitivity_level:
                    newCellTypes[cI] = after_type
            elif sign_pos_neg == -1:
                if score_max < sensitivity_level:
                    newCellTypes[cI] = after_type
            else:
                if pass_no >= 2:
                    newCellTypes[cI] = after_type

    # Assign NONE to empty cell types
    for i in range(len(newCellTypes)):
        if newCellTypes[i] == '':
            newCellTypes[i] = 'NONE'

    # Update the allCellTypes array
    allCellTypes[:, pass_no] = newCellTypes

    return allCellTypes

# Purpose: This function 'uniqueStrCell' performs 'UNIQUE' for a list of string and numeric values.

# Inputs: 
#     inputStrCell - a list of string and numeric values

# Outputs:
#     out - a list of unique string and numeric values (converted to strings)

# Description:
#     This function takes a list of inputs (which can be strings, numbers, or empty/NaN values), 
#     and returns a list of unique strings and numbers (converted to strings). 
#     It excludes empty and NaN values.


def uniqueStrCell(inputStrCell):
    out = []
    
    A = [isinstance(x, str) for x in inputStrCell]
    B = [isinstance(x, (int, float)) for x in inputStrCell]
    C = [x is None or np.isnan(x) for x in inputStrCell]
    D = [x == '' for x in inputStrCell]
    
    numCell = [str(x) for x in inputStrCell if type(x) in (int, float) and not np.isnan(x) and x != '']
    numCell = list(set(numCell))
    
    strCell = list(set([x for x in inputStrCell if isinstance(x, str) and x != '']))
    
    out = strCell + numCell
    
    return out
