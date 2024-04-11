# Importing the libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 

# Configuration
cluster = 'Gulliver_27'
coords = [146.0801, -54.1154]

# cluster = 'OC-0541'
# coords = [155.133016, -59.700489]

os.makedirs(f'data_saved/{cluster}', exist_ok=True)
os.makedirs(f'plots/{cluster}', exist_ok=True)
os.makedirs(f'data_output/{cluster}', exist_ok=True)

# Importing the datasets
gaia = pd.read_csv(f'data_input/{cluster}/star_table_{cluster}.csv')
uv = pd.read_csv(f'data_input/{cluster}/uv.csv')

members_m = pd.read_csv('data_input/members_m_final.csv')
members_p = pd.read_csv('data_input/members_p_final.csv')

main_sequence_1 = pd.read_csv('data_input/lines_UBPRPG_1.txt', delimiter = '\s')
main_sequence_2 = pd.read_csv('data_input/lines_UBPRPG_2.txt', delimiter = '\s')

gaia['id'] = range(1, len(gaia) + 1)   
uv['id'] = range(1, len(uv) + 1)   

#########################################
gaia_matched_list = pd.DataFrame(columns=['id', 'x', 'y', 'source_id', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'ra', 'dec'])
uv_matched_list = pd.DataFrame(columns=['id', 'x', 'y', 'u', 'source_id', 'u-bp'])
gaia_filtered_list = pd.DataFrame(columns=['id', 'x', 'y', 'source_id', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'ra', 'dec'])
uv_filtered_list = pd.DataFrame(columns=['id', 'x', 'y', 'u', 'source_id', 'u-bp'])

gaia_mags = ['phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag'] # choosing filter for Gaia
for gaia_mag in gaia_mags:
    filter = gaia_mag.replace('phot_', '').replace('_mean_mag', '')

    gaia_selected = gaia[['id', 'ra', 'dec', gaia_mag]]
    gaia_selected.to_csv(f'data_output/{cluster}/gaia_sky_{filter}.dat', sep='\t', index=False, header=False)

    uv_selected = uv[['id', 'ra', 'dec', 'u']]
    uv_selected.to_csv(f'data_output/{cluster}/uv_sky_{filter}.dat', sep='\t', index=False, header=False)

    # Changing coordinates from RA-DEC to x-y
    os.system(f'project_coords data_output/{cluster}/gaia_sky_{filter}.dat 1 2 {coords[0]} {coords[1]} asec outfile=data_output/{cluster}/gaia_cart_{filter}.dat')
    os.system(f'project_coords data_output/{cluster}/uv_sky_{filter}.dat 1 2 {coords[0]} {coords[1]} asec outfile=data_output/{cluster}/uv_cart_{filter}.dat')

    # Importing datasets with new coordinates
    gaia_cart = pd.read_csv(f'data_output/{cluster}/gaia_cart_{filter}.dat', sep='\s+', header=None, names=['id', 'x', 'y', gaia_mag])
    uv_cart = pd.read_csv(f'data_output/{cluster}/uv_cart_{filter}.dat', sep='\s+', header=None, names=['id', 'x', 'y', 'u'])

    # Matching the datasets
    os.system(f'match data_output/{cluster}/gaia_cart_{filter}.dat 1 2 3 data_output/{cluster}/uv_cart_{filter}.dat 1 2 3 outfile=data_output/{cluster}/matched_{filter} id1=0 id2=0 matchrad=3.0 trirad=0.001 nobj=40 recalc')

    # Importing matched datasets
    gaia_matched = pd.read_csv(f'data_output/{cluster}/matched_{filter}.mtA', sep='\s+', header=None, names=['id', 'x', 'y', gaia_mag])
    uv_matched = pd.read_csv(f'data_output/{cluster}/matched_{filter}.mtB', sep='\s+', header=None, names=['id', 'x', 'y', 'u'])

    # List of columns to add
    columns_to_add = ['source_id', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'ra', 'dec']

    # Remove the column represented by the variable gaia_mag
    gaia_matched = gaia_matched.drop(columns=gaia_mag)

    # Merge the dataframes on the 'id' column
    gaia_matched = pd.merge(gaia_matched, gaia[['id'] + columns_to_add], on='id')

    uv_matched = pd.concat([uv_matched, gaia_matched['source_id']], axis=1)
    uv_matched['u-bp'] = uv_matched['u'] - gaia_matched['phot_bp_mean_mag'] # calculating u-bp color

    members_comb = pd.concat([members_m, members_p], ignore_index=True)
    members_comb_filtered = members_comb[members_comb['probability'] >= 0.7]

    mask = gaia_matched['source_id'].isin(members_comb_filtered['source_id'])
    gaia_filtered = gaia_matched[mask]

    mask = uv_matched['source_id'].isin(members_comb_filtered['source_id'])
    uv_filtered = uv_matched[mask]
    uv_filtered['u-bp'] = uv_filtered['u'] - gaia_matched['phot_bp_mean_mag'] # calculating u-bp color

    gaia_matched_list = pd.concat([gaia_matched_list, gaia_matched], ignore_index=True)
    uv_matched_list = pd.concat([uv_matched_list, uv_matched], ignore_index=True)
    gaia_filtered_list = pd.concat([gaia_filtered_list, gaia_filtered], ignore_index=True)
    uv_filtered_list = pd.concat([uv_filtered_list, uv_filtered], ignore_index=True)

gaia_matched_list.drop_duplicates(subset='source_id', keep='first', inplace=True)
uv_matched_list.drop_duplicates(subset='source_id', keep='first', inplace=True)
gaia_filtered_list.drop_duplicates(subset='source_id', keep='first', inplace=True)
uv_filtered_list.drop_duplicates(subset='source_id', keep='first', inplace=True)

print(len(gaia_matched_list), len(uv_matched_list), len(gaia_filtered_list), len(uv_filtered_list))