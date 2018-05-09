"""
BST 281 Final Project

Produces new data files: 
transposed_data_subset_genus_or_species_level.tsv
transposed_data_subset_species_level.tsv
hmp2_metadata_metagenomics_CD_subset.tsv
hmp2_metadata_metagenomics_nonIBD_subset.tsv
hmp2_metadata_metagenomics_CD_subset_transposed.tsv
hmp2_metadata_metagenomics_nonIBD_subset_transposed.tsv
external_IDs_and_sequencing_reads_count_case_H4015.txt
external_IDs_and_sequencing_reads_count_control_H4023.txt

4/28/18 - MC
"""
import pandas as pd
import gzip
import subprocess

def get_only_genus_or_species_level(name):

    if ('s__' not in str(name) and 'g__' in str(name) and 'k__Bacteria' in str(name)) or ('t__' not in str(name) and 's__' in str(name) and 'k__Bacteria' in str(name)):
        return name

def get_only_species_level(name):

    if ('t__' not in str(name) and 's__' in str(name) and 'k__Bacteria' in str(name)):
        return name

def transposed_data_subset_helper(data, column_names, row_names):
    
    row_names = list(filter(None, row_names))

    new_data = data.loc[row_names, column_names]

    new_data_transpose = new_data.transpose()

    return new_data_transpose

def transposed_data_subset_genus_or_species_level(data, column_names):
    """
    Creates transposed_data_subset_genus_or_species_level.tsv
    """

    all_row_names = pd.Series(data.index)
    row_names = list(all_row_names.apply(get_only_genus_or_species_level))
    
    new_data_transpose = transposed_data_subset_helper(data, column_names, row_names)
    new_data_transpose.to_csv("transposed_data_subset_genus_or_species_level.tsv", sep="\t")

def transposed_data_subset_species_level(data, column_names):
    """
    Creates transposed_data_subset_species_level.tsv
    """

    all_row_names = pd.Series(data.index)
    row_names = list(all_row_names.apply(get_only_species_level))
    
    new_data_transpose = transposed_data_subset_helper(data, column_names, row_names)
    new_data_transpose.to_csv("transposed_data_subset_species_level.tsv", sep="\t")






def metadata_subset(data, metadata):

    external_ID_list = list(metadata["External ID"])
    week_num_list = list(metadata["week_num"])
    participant_ID_list = list(metadata["Participant ID"])
    sampleID_list = list(data)
    external_IDs = []
    week_nums = []
    participant_IDs = []
    for i in range(len(external_ID_list)):
        external_ID = external_ID_list[i]
        week_num = week_num_list[i]
        participant_ID = participant_ID_list[i]
        if external_ID in sampleID_list:
            external_IDs.append(external_ID)
            week_nums.append(week_num)
            participant_IDs.append(participant_ID)


    external_ID_week_num_df = pd.DataFrame(dict(zip(external_IDs, week_nums)), index=['week_num'])
    external_ID_participant_ID_df = pd.DataFrame(dict(zip(external_IDs, participant_IDs)), index=['Participant ID'])

    external_ID_week_num_participant_ID_df = external_ID_week_num_df.append(external_ID_participant_ID_df)

    taxonomic_profiles_subset = data.loc[:,external_IDs]
    
    taxonomic_profiles_subset_week_num_pariticipant_id = external_ID_week_num_participant_ID_df.append(taxonomic_profiles_subset)
    
    return taxonomic_profiles_subset_week_num_pariticipant_id


def create_metadata_subset_data(data):
    """
    Creates: 
    taxonomic_profiles_metagenomics_CD_subset_week_num_participant_id.tsv
    taxonomic_profiles_metagenomics_nonIBD_subset_week_num_participant_id.tsv
    """

    metadata_CD = pd.read_csv("hmp2_metadata_metagenomics_CD_subset.csv", sep=",", low_memory=False)
    metadata_nonIBD = pd.read_csv("hmp2_metadata_metagenomics_nonIBD_subset.csv", sep=",", low_memory=False)

    taxonomic_profiles_CD_subset_week_num_participant_id = metadata_subset(data, metadata_CD)
    taxonomic_profiles_CD_subset_week_num_participant_id.to_csv("taxonomic_profiles_metagenomics_CD_subset_week_num_participant_id.tsv", sep="\t")

    taxonomic_profiles_nonIBD_subset_week_num_participant_id = metadata_subset(data, metadata_nonIBD)
    taxonomic_profiles_nonIBD_subset_week_num_participant_id.to_csv("taxonomic_profiles_metagenomics_nonIBD_subset_week_num_participant_id.tsv", sep="\t")

    return (taxonomic_profiles_CD_subset_week_num_participant_id, taxonomic_profiles_nonIBD_subset_week_num_participant_id)


def create_metadata_subset_data_transposed(data):
    """
    Creates: 
    taxonomic_profiles_metagenomics_CD_subset_week_num_participant_id_transposed.tsv
    taxonomic_profiles_metagenomics_nonIBD_subset_week_num_participant_id_transposed.tsv
    """

    (taxonomic_profiles_CD_subset_week_num_participant_id, taxonomic_profiles_nonIBD_subset_week_num_participant_id) = create_metadata_subset_data(data)
    
    taxonomic_profiles_CD_subset_week_num_participant_id_transposed = taxonomic_profiles_CD_subset_week_num_participant_id.transpose()
    taxonomic_profiles_nonIBD_subset_week_num_participant_id_transposed = taxonomic_profiles_nonIBD_subset_week_num_participant_id.transpose()

    taxonomic_profiles_CD_subset_week_num_participant_id_transposed.to_csv("taxonomic_profiles_metagenomics_CD_subset_week_num_participant_id_transposed.tsv", sep="\t")
    taxonomic_profiles_nonIBD_subset_week_num_participant_id_transposed.to_csv("taxonomic_profiles_metagenomics_nonIBD_subset_week_num_participant_id_transposed.tsv", sep="\t")






def sequencing_reads_counts(external_IDs):
    """
    Returns a pandas DataFrame of external IDs and sequencing reads counts
    """

    sequencing_reads_counts = []
    
    for external_ID in external_IDs:

        file_name = "raw_data/" + external_ID + ".fastq.gz"
        num_lines_CompletedProcess = subprocess.run("gzip -dc " + file_name + "| wc -l", shell=True, stdout = subprocess.PIPE, universal_newlines = True)
        num_lines_int = int(num_lines_CompletedProcess.stdout.split()[0])
        sequencing_reads_count = int(num_lines_int / 4)

        sequencing_reads_counts.append(sequencing_reads_count)

    external_IDs_and_sequencing_reads_counts = pd.DataFrame(sequencing_reads_counts, external_IDs)
    external_IDs_and_sequencing_reads_counts.columns = ['count']

    return external_IDs_and_sequencing_reads_counts




def main():

    data = pd.read_csv("taxonomic_profiles.tsv", sep="\t", index_col=0)

    control_list = ['HSM67VDR_P', 'HSM6XRUL', 'HSM6XRUN', 'HSM6XRUR', 'HSM6XRQ8', 'HSM7CZ16', 'HSM7CZ18', 'HSM7CZ1A', 'HSM7CZ1C', 'HSM7CZ1E', 'HSM7CZ1G', 'HSM7J4HA', 'HSM7J4HC', 'HSM7J4HE', 'HSM7J4HG', 'HSM7J4HI', 'HSM7J4HK', 'HSM7J4KC', 'HSM7J4KG', 'HSM7J4KI', 'HSM7J4KK', 'HSM7J4KM']
    cases_list = ['HSM5MD5X_P', 'HSM5MD62', 'HSM5MD6Y', 'HSM5MD71', 'HSM5MD73', 'HSM5MD75', 'HSM6XRS4', 'HSM6XRS6', 'HSM6XRS8', 'HSM6XRSE', 'HSM7CYZ5', 'HSM7CYZ7', 'HSM7CYZ9', 'HSM7CYZB', 'HSM7CYZD', 'HSM7CYZF', 'HSM7J4QB', 'HSM7J4QD', 'HSM7J4QF', 'HSM7J4QH', 'HSM7J4QJ', 'HSM7J4QL']
    cases_and_control_list = control_list + cases_list

    # Creates: 
    # transposed_data_subset_genus_or_species_level.tsv
    transposed_data_subset_genus_or_species_level(data, cases_and_control_list)
    # Creates:
    # transposed_data_subset_species_level.tsv
    transposed_data_subset_species_level(data, cases_and_control_list)

    # Creates: 
    # taxonomic_profiles_metagenomics_CD_subset_week_num_participant_id.tsv
    # taxonomic_profiles_metagenomics_nonIBD_subset_week_num_participant_id.tsv
    create_metadata_subset_data(data)
    
    # Creates: 
    # taxonomic_profiles_metagenomics_CD_subset_week_num_participant_id_transposed.tsv
    # taxonomic_profiles_metagenomics_nonIBD_subset_week_num_participant_id_transposed.tsv
    create_metadata_subset_data_transposed(data)
    
    # Creates:
    # external_IDs_and_sequencing_reads_count_case_H4015.txt
    external_IDs_and_sequencing_reads_count = sequencing_reads_counts(cases_list)
    external_IDs_and_sequencing_reads_count.to_csv("external_IDs_and_sequencing_reads_count_case_H4015.csv", sep=",")
    # Creates:
    # external_IDs_and_sequencing_reads_count_control_H4023.txt
    external_IDs_and_sequencing_reads_count = sequencing_reads_counts(control_list)
    external_IDs_and_sequencing_reads_count.to_csv("external_IDs_and_sequencing_reads_count_control_H4023.csv", sep=",")


if __name__ == "__main__":
    main()