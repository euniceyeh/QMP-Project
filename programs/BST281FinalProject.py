# BST 281 Final Project
#
# Notes:
# Produces transposed_data_subset.tsv 4/18/18 - MC

import pandas as pd

data = pd.read_csv("taxonomic_profiles.tsv", sep="\t", index_col=0)

control_list = ['HSM67VDR_P', 'HSM6XRUL', 'HSM6XRUN', 'HSM6XRUR', 'HSM6XRQ8', 'HSM7CZ16', 'HSM7CZ18', 'HSM7CZ1A', 'HSM7CZ1C', 'HSM7CZ1E', 'HSM7CZ1G', 'HSM7J4HA', 'HSM7J4HC', 'HSM7J4HE', 'HSM7J4HG', 'HSM7J4HI', 'HSM7J4HK', 'HSM7J4KC', 'HSM7J4KG', 'HSM7J4KI', 'HSM7J4KK', 'HSM7J4KM']
cases_list = ['HSM5MD5X_P', 'HSM5MD62', 'HSM5MD6Y', 'HSM5MD71', 'HSM5MD73', 'HSM5MD75', 'HSM6XRS4', 'HSM6XRS6', 'HSM6XRS8', 'HSM6XRSE', 'HSM7CYZ5', 'HSM7CYZ7', 'HSM7CYZ9', 'HSM7CYZB', 'HSM7CYZD', 'HSM7CYZF', 'HSM7J4QB', 'HSM7J4QD', 'HSM7J4QF', 'HSM7J4QH', 'HSM7J4QJ', 'HSM7J4QL']
column_names = control_list + cases_list

def get_only_until_genus(name):
    if 's__' not in str(name) and 'g__' in str(name) and 'k__Bacteria' in str(name):
        return name

all_row_names = pd.Series(data.index)
row_names = list(all_row_names.apply(get_only_until_genus))
row_names = list(filter(None, row_names))

new_data = data.loc[row_names, column_names]

new_data_transpose = new_data.transpose()

new_data_transpose.to_csv("transposed_data_subset.tsv", sep="\t")





