#!/usr/bin/env python3

import sys
import itertools
import pandas as pd

def extract_sequence_4df(pddf, curr_species1):
    curr_species1 = str(curr_species1)
    sequence_curr_species1 = pddf.loc[pddf['species'] == curr_species1] # Extracting sequence for species
    sequence_curr_species1 = str(sequence_curr_species1['sequence'].values)
    sequence_curr_species1 = sequence_curr_species1.replace('\'','').replace('[','').replace(']','')
    sequence_curr_species1 = sequence_curr_species1.upper()
    return sequence_curr_species1

def counter_dict_conv_values(all_species_counter1):
    for keys1, values1 in all_species_counter1.items():
        keys1 = str(keys1)
        for keys2, values2 in values1.items():
            keys2 = str(keys2)
            for keys3, values3 in values2.items():
                keys3 = str(keys3)
                values3_int = int(str(values3).replace('[','').replace(']',''))
                all_species_counter1[keys1][keys2][keys3] = values3_int
    return all_species_counter1

def pairspecies_counter_dict_conv_values(all_species_counter1):
#    for keys1, values1 in all_species_counter1.items():
#        keys1 = str(keys1)
#    for keys2, values2 in values1.items():
    for keys2, values2 in all_species_counter1.items():
        keys2 = str(keys2)
        for keys3, values3 in values2.items():
            keys3 = str(keys3)
            values3_int = int(str(values3).replace('[','').replace(']',''))
#            all_species_counter1[keys1][keys2][keys3] = values3_int
            all_species_counter1[keys2][keys3] = values3_int
    return all_species_counter1

def counter_dict_2_df(nested_dict_counter):
    ids_listX = []
    df1 = []

    for ids_list1, currX1 in nested_dict_counter.items():
        ids_listX.append(ids_list1)
        df1.append(pd.DataFrame.from_dict(currX1, orient='index'))

    final_df = pd.concat(df1, axis=0, keys=ids_listX, sort=False) # horizontally concatenated DF
    cols_to_sum = final_df.columns[0 : final_df.shape[1]] # getting all columns

    final_df['Z_total'] = final_df[cols_to_sum].sum(axis=1) # creating summation column
    final_df.reset_index(inplace=True) # Converting indices to columns
    final_df = final_df.rename(columns = {'level_0':'Species', 'level_1':'Pattern'})

    final_df = final_df.melt(id_vars=["Species","Pattern"], var_name="Transition", value_name="Count") # Reshaping the DF
    final_df = final_df.sort_values(["Species","Pattern","Transition"], ascending = (True, True, True)) # Sorting the DF
    final_df = final_df.reset_index(drop=True) # Changing the index after sorting
    return final_df

def pairwise_counter_dict_2_df(pairwise_nested_dict_counter):
    ids_listX = []
    df1 = []

    for ids_list1, currX1 in pairwise_nested_dict_counter.items():
        ids_listX.append(ids_list1)
        df1.append(pd.DataFrame.from_dict(currX1, orient='index'))

    final_df = pd.concat(df1, axis=0, keys=ids_listX, sort=False) # horizontally concatenated DF

    final_df.reset_index(inplace=True) # Converting indices to columns
    final_df = final_df.rename(columns = {'level_0':'Speciespairs', 'level_1':'Patterndetails'})

    final_df = final_df.melt(id_vars=["Speciespairs","Patterndetails"], value_name="Count") # Reshaping the DF
    final_df.drop('variable', axis=1, inplace=True) # deleting unnecessary column
    final_df = final_df.sort_values(["Speciespairs","Patterndetails","Count"], ascending = (True, True, True)) # Sorting the DF
    final_df = final_df.reset_index(drop=True) # Changing the index after sorting
    return final_df

