'''
DESCRIPTION:
Script containing useful functions for processing SIGNAL-seq data. 
'''

#  imports
import numpy as np
import pandas as pd
import scanpy as sc
import scprep 

# Perform EMD on all cells as global control
def calculate_emd(emd_data, control_obs, control_condition):

    # Generate condition list
    condition_list = list(emd_data[control_obs].unique())


    # Generate the marker list
    # Pop the sample_id col at index 0
    marker_list = list(emd_data.columns)
    marker_list.pop(0)

    # Initialise empty datframe to contain all signed EMD values
    adt_signed_emds = pd.DataFrame()

    # Callums code to run EMD
    for condition in condition_list:

        # Initialise and empty array to hold EMD clas for each condition
        # These will ultimately be columns in the dataframe
        emd_array = []

        # Isolate query data based on condition
        query_filter_df = emd_data[emd_data[control_obs].str.contains(condition)]

        # Isolate reference data, here it is based on the control_condition
        reference_filter_df = emd_data[emd_data[control_obs].str.contains(control_condition)]

        for marker in marker_list:

                        # Calculate the sign for the EMD representing the direction distribution by media
                        sign = np.sign(query_filter_df[marker].median() - reference_filter_df[marker].median())

                        # If median is 0 use the mean
                        if sign == 0:
                            sign = np.sign(query_filter_df[marker].mean() - reference_filter_df[marker].mean())
                        
                        # runt the EMD and create sign
                        signed_emd = sign*scprep.stats.EMD(
                            query_filter_df[marker], reference_filter_df[marker])

                        #Append the array values    
                        emd_array.append(signed_emd)

        # Append to the EMD dataframe
        adt_signed_emds[condition]=pd.Series(emd_array)
        
    # Add row names
    adt_signed_emds.index = list(marker_list)

    # Check for any NAs
    assert not adt_signed_emds.isna().values.any()

    return(adt_signed_emds)


#### ADT barcode merging functions for SPLiT-seq

# RT barcode merding funtion based on Kuchina et al. SPLiT-seq setup
def merge_rt_barcodes(adata, cell_barcode_1, cell_barcode_2):

  #Get index for both cell barcodes
  index_1 = adata.obs.index == cell_barcode_1
  index_2 = adata.obs.index == cell_barcode_2

  # merge cell barcodes data
  adata.X[index_1,:] += adata.X[index_2,:]

  # Remove cell_barcode_2 from the adata object
  adata = adata[~index_2,:].copy()
    
  return(adata)

# RT barcode matching funtion based on Kuchina et al. SPLiT-seq setup
def pairwise_rt_comparison(ids):
    
    # Iterate throught the list of Ligation barcoiodes
    # Find RT matches with 48 difference
    matched = []
    for i in range(len(ids)):
        for j in range(i, len(ids)):
            if ids[i] - ids[j] == -48:
                matched.append([i, j])
    
    return(matched)

