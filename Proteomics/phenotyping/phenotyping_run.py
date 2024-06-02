import os
import pandas as pd
from cell_phenotyping_functions import *

def main():
    # Define the paths to various configuration and data files
    config_file_path = 'config_template.xlsx'
    additional_dilation_table_path = 'additional_dilation.xlsx'
    metalJMFtablePath = os.path.join('ConfigFiles', 'metal_indicator_JMF_kmeans.xlsx')
    outcomeJMFPath = 'JMFPath'
    
    # Define output folders
    cIMatOuputFolder = 'presence_matrices'
    output_figure_folder = 'output_folder'
    allCellTypesAddress = os.path.join('ConfigFiles', 'allCellTypesRGB.xlsx')
    ring_rad = 3
    percentile = 0.995
    run_step_JMF = 1
    
    # Read the cell type labels and their RGB values
    allLabelsTable = pd.read_excel(allCellTypesAddress)
    allLabels = allLabelsTable['CellTypes']
    R = allLabelsTable['R']
    G = allLabelsTable['G']
    B = allLabelsTable['B']
    type_colors = np.array([R, G, B]).T / 255
    
    run_celltype_assignment = 1
    
    # Extract the file name from the config file path
    configFileName = Path(config_file_path).stem
    
    # Define the full path to the config file
    config_path = os.path.join('ConfigFiles', config_file_path)
    
    # Read various tables and file lists from the config file
    (IndicatorsList, nuc_csv_source_folder, ring_csv_source_folder, rule_table_location, nOfGroups, 
     groupPatientIDTableData, headerGroupPatientID, fileNameColumn, imageIDColumn, regIndices_files, 
     pano_files, IndicatorsListPath, input_pano_folder, input_cell_segmentation_folder, 
     output_folder) = read_tables_list_files(config_path)
    
    middleString = []
    files = []
    
    # Obtain threshold values if not running the JMF step
    if run_step_JMF == 0:
        metalListThreshVals = obtain_thresh_vals(metalJMFtablePath)
    
    # Create the outcome directory if running the JMF step and it doesn't exist
    if run_step_JMF == 1 and not os.path.exists(outcomeJMFPath):
        os.makedirs(outcomeJMFPath)
    
    # Loop through each file in pano_files (adjust the range as needed)
    for indexFile in range(1):
        panoFilePath = pano_files[indexFile]['name']
        folderPath, scanName = os.path.split(panoFilePath)
        
        # Find the index of the current scan name in the fileNameColumn
        T = [i for i, item in enumerate(fileNameColumn) if scanName in item]
        if T:
            # Extract the folder and file names for the path structure
            pathStructure = panoFilePath[len(input_pano_folder):]
            folderName, fileName = os.path.split(pathStructure)
            cellTypeToSavePath = os.path.join(output_folder, 'MATData', folderName, f'{fileName}.mat')
            
            folderStructure = folderPath[len(input_pano_folder):]
            folderPathReg = os.path.join(input_cell_segmentation_folder, folderStructure, scanName)
            regIndPath = os.path.join(folderPathReg, 'allRegionIndices.mat')
            ringIndPath = os.path.join(folderPathReg, 'allRingIndices.mat')
            nucleiFilePath = os.path.join(folderPathReg, 'nuclei_multiscale.mat')
            
            found_group = ''
            # Find the group for the current image and corresponding rule table location
            for gI in range(nOfGroups):
                curGr = groupPatientIDTableData[:, gI]
                for rowItr in range(len(curGr)):
                    if curGr[rowItr] == imageIDColumn[T]:
                        groupName = headerGroupPatientID[gI][len('Group_'):]
                        for ri in range(len(rule_table_location)):
                            if rule_table_location[ri, 0] == groupName:
                                rule_table_location_for_this_group = rule_table_location[ri, 1]
                                found_group = groupName
            
            # Apply mask segmentation if conditions are met
            if (not os.path.exists(cellTypeToSavePath) and run_step_JMF == 1) or (run_celltype_assignment == 0 and run_step_JMF == 1):
                additional_dilation_table = pd.read_excel(additional_dilation_table_path)
                print('mask segmentation')
                apply_mask_segmentation_files_weights(panoFilePath, metalJMFtablePath, outcomeJMFPath, additional_dilation_table)
            
            # Assign cell types and create figures if required
            if run_celltype_assignment:
                allCellTypes = main_cell_type_assignment_v3(rule_table_location_for_this_group, outcomeJMFPath, scanName, regIndPath, ringIndPath)
                cellTypeToSavePath = save_celltypes(allCellTypes, panoFilePath, input_pano_folder, input_cell_segmentation_folder, output_folder)
                create_figure_here(cellTypeToSavePath, output_folder, input_cell_segmentation_folder, output_figure_folder, allLabels, type_colors)
        else:
            print(f'The program is not producing a celltype for file: [{scanName}]')

if __name__ == "__main__":
    main()