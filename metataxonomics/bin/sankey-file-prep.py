#!/usr/bin/python3

# Author: Harm Laurense
# Last edited: 13-11-2021
# Function: This script is used to create the necessary .csv files to create sankey diagram(s) using R.
# Bugs:
# TODO: This script potentially needs altering for it to be implemented into nextflow.
# TODO: Bugtesting; test what happens when input of files is altered for example.

"""
Parameter 1: Integer/Float for filtering
Parameter 2: String (boolean principle*) for creating a .csv file for each individual sample
* Checks if input is equal to "true". Booleans are a bit tricky / quick to error in command lines flags.
Parameter 3: String input for the metadata file
Parameter 4: String (booleon principle*) for creating a .csv file for each unique rankstat combination of samples
Parameter 5: String input for the biotaviz file.

Typical run::
python3 .\sankey-file-prep.py 0.01 false Metadata.tsv false relative-table_kopie_test.biotaviz_test.txt
"""

import csv
import sys
import traceback
import itertools

# Dictionary used to determine the numbers for linking nodes (based on taxonomic rank)
taxonomic_ranks_dict = {
    "empty": -1,
    "phylum": 1,
    "class": 2,
    "order": 3,
    "family": 4,
    "genus": 5,
    "species": 6}
# Dictionary used to determine the numbers for linking nodes (based on taxonomic rank / last used taxonomic rank)
taxonomic_ranks_last_linked_rank_dict = {
    "phylum": 0,
    "class": 0,
    "order": 0,
    "family": 0,
    "genus": 0,
    "species": 0}

def main(tax_filter, sample_repeat, mappingfile, combine_rankstat):
    """
    Determines what functions are needed to be called based on command line input.
    :param mappingfile: File containing the metadata. Determines which samples are averaged.
    :param tax_filter: Parameter for filtering (low) relative abundance.
    :param sample_repeat: Parameter which determines if files are created for every individual sample.
    :return: .csv files according to user input, to be used in R script for creating the sankey diagrams
    """
    # Create a .csv file for every sample (needed for generating sankey diagram in R script)
    if sample_repeat.lower() == "true":
        total_samples = determine_sample_total()
        for sample in range(total_samples):
            hierarchy_counts(tax_filter, sample, average_all_samples, filename_combination)

    # Create .csv file iterative over all possible rankstat combinations
    if combine_rankstat.lower() == "true":
        # Preparation for creating .csv files for averaged samples (based on the Metadata file)
        all_sets = get_sets(mappingfile)
        sample_index = determine_sample_index()
        combinations_rankstatheaders, unique_sample_combinations = determine_all_rankstat_combinations(all_sets)
        indexed_combinations = combination_to_index(sample_index, unique_sample_combinations)

        # Create sample average file for each unique sample combination
        for index, combination in enumerate(indexed_combinations):
            sample_average(combination, combinations_rankstatheaders[index])

    # Create sample average file for every individuel rankstat column
    all_sets = get_sets(mappingfile)
    sample_index = determine_sample_index()
    filename_rankstatheaders, rankstat_samples = determine_rankstat_samples(all_sets)
    indexed_combinations = combination_to_index(sample_index, rankstat_samples)
    for index, combination in enumerate(indexed_combinations):
        sample_average(combination, filename_rankstatheaders[index])

    # Create sample average file over all samples (includes samples without metadata values)
    sample_average_all()


def determine_sample_total():
    """
    Determine the total amount of samples by counting the columns. The first 2 columns aren't samples and thus skipped.
    :return:Number of samples (total)
    """
    try:
        with open(biotavizfile, "r") as file:
            line = file.readline()
            total_samples = len(line.rstrip().split('\t')[2:])
        return total_samples
    except IndexError:
        print("IndexError; check if the correct file is given as input: ", traceback.print_exc())


def hierarchy_counts(tax_filter, sample, average_samples, filename_combination):
    """
    Generate the files necessary for creating a sankey diagram from the biotaviz file.
    :param tax_filter: Parameter for filtering (low) relative abundance.
    :param sample: Integer (standard 0) used as index. This only changes if sample_repeat is set to true.
    :param average_samples: List of average values (currently from all rankstat sample combinations).
    :param filename_combination: Specific string correlating to the unique combination, needed for unique filenames.
    :return: Variables (link1, link2, label) are determined and finally given to the write_new_biotaviz_file() function.
    Which in return will write the necessary .csv files.
    """
    link1 = ["link1", "link1"]
    link2 = ["link2", "link2"]
    label = [["label", "value"]]
    count = 0
    last_rank = "empty"
    # Currently used only for testing purposes (tax_filter control)
    removed_entries = []
    try:
        with open(biotavizfile, "r") as file:
            for _ in range(1):  # skip column headers + root
                next(file)
            taxonomic_rank = []
            for index, line in enumerate(file.readlines()):
                if line.strip():
                    line = line.rstrip().split('\t')
                    if len(average_samples) > 0:
                        tax_value = float(average_samples[index])
                    else:
                        # Skip first 2 columns
                        tax_value = float(line[sample + 2])
                    # Values of 0 (relative abundance) or below taxonomic filter (standard 1%) aren't used
                    if tax_value >= tax_filter and tax_value > 0:
                        tax_rank = line[1].split('-')[0].rstrip()
                        tax_specific = line[1].split('-')[1].rstrip()
                        taxonomic_rank.append(tax_rank)
                        tax_with_score = tax_specific + ":" + str(round(float(tax_value) * 100, 2)) + "%"
                        # This var replacement was used during testing for various value sizes based on taxonomic rank
                        # tax_value = taxonomic_ranks_valuesize_dict.get(tax_rank)*tax_value
                        label.append([tax_with_score, tax_value])
                    else:
                        removed_entries.append(line)
            # The following is the logic to determine the number combinations for connecting the nodes
            for rank in taxonomic_rank:
                if rank == "domain":
                    continue
                elif taxonomic_ranks_dict.get(rank):
                    count += 1
                    if rank == last_rank:
                        link1.append(taxonomic_ranks_last_linked_rank_dict.get(rank))
                    if taxonomic_ranks_dict.get(last_rank) < taxonomic_ranks_dict.get(rank):
                        if link2[-1] != "link2":
                            link1.append(link2[-1])
                        else:
                            link1.append(0)
                    elif taxonomic_ranks_dict.get(last_rank) > taxonomic_ranks_dict.get(rank):
                        link1.append(taxonomic_ranks_last_linked_rank_dict.get(rank))
                    link2.append(count)
                    last_rank = rank
                    taxonomic_ranks_last_linked_rank_dict.update({rank: link1[-1]})
        if len(label) > 1:
            # print(removed_entries)
            write_new_biotaviz_file(link1, link2, label, sample, average_samples, filename_combination)
        else:
            sys.exit(print(
                "No matches with current criteria found, try lowering the given filter for "
                "relative abundance"))
    except IndexError:
        print("IndexError; check if the correct file is given as input: ", traceback.print_exc())


def write_new_biotaviz_file(link1, link2, label, sample, average_samples, filename_combination):
    """
    Write a .csv file containing the necessary information for creating a sankey diagram (used in: Sankey R module)
    :param link1: List of numbers which represents the node being connected from.
    :param link2: List of (second) numbers which represents the node connected to.
    :param label: List of lists containing the labels (taxonomic rank : % abundance) and value (relative abundance)
    :param sample: Integer (index) used to create an unique filename
    :param average_samples: Needed to determine which filename is to be used.
    :param filename_combination: Specific string correlating to the unique combination, needed for unique filenames.
    :return: .csv file (4 columns)
    """
    try:
        if not average_samples:
            filename = f"biotaviz_sankey_prepfile-{sample}.csv"
        else:
            filename = f"biotaviz_sankey_prepfile-{filename_combination}.csv"
        with open(filename, "w", newline="") as f:
            wr = csv.writer(f)
            # column headers
            for index, number in enumerate(link1):
                wr.writerow([number, link2[index], label[index][0], label[index][1]])
    except IndexError:
        print("IndexError; can't write to biotaviz_sankey_prepfile-{sample}.csv: ", traceback.print_exc())


def get_sets(filename):
    """
    Determine what subsets should be compared (rankstat)
    :param filename: Metadata file (should contain rankstat columns)
    :return: Dictionary containing each rankstat column, containing all unique values with corresponding samples
    """
    lines = load_txt(filename).strip().split('\n')
    headers = {}
    datasets = {}
    for i, element in enumerate(lines[0].split('\t')):
        if element.split('_')[0].lower() == 'rankstat':
            dataset = element.split("_")[1]
            headers[i] = dataset
            datasets[dataset] = {}
    for line in lines[1:]:
        lineg = line.split('\t')
        sample_id = lineg[0]
        for i, element in enumerate(lineg):
            if i in list(headers.keys()) and element != '':
                sample_class = element
                dataset = headers[i]
                if sample_class not in datasets[dataset]:
                    datasets[dataset][sample_class] = []
                datasets[dataset][sample_class].append(sample_id)
    return datasets


def load_txt(filename):
    """
    Opens files and returns their input/content
    :param filename: File
    :return: File content
    """
    fileinput = open(filename)
    text = fileinput.read()
    fileinput.close()
    return text


def determine_sample_index():
    """
    Determines the index for each rankstat column, necessary in following functions to determine the sample average
    :return: Dictionary containing the rankstat columns and their corresponding index.
    """
    try:
        sample_index = {}
        with open(biotavizfile, "r") as file:
            sample_headers = file.readline()  # skip column headers + root
            sample_headers = sample_headers.rstrip().split("\t")
            for index, header in enumerate(sample_headers[2:]):
                sample_index.update({header: index + 2})
        return sample_index
    except IndexError:
        print("IndexError; check if the correct file is given as input: ", traceback.print_exc())


def determine_rankstat_samples(all_sets):
    """
    Create two lists necessary for unique filenames and to determine index + average of samples in following functions
    :param all_sets: Dictionary containing each rankstat column, containing all unique values with corresponding samples
    :return: A list containing the filenames (rankstat+unique_value) and a list with their corresponding samples
    """
    rankstat_samples = []
    filename_rankstatheaders = []
    for keys, values in all_sets.items():
        for item in values.items():
            column_value = str(keys) + "-" + str(item[0])
            filename_rankstatheaders.append(column_value)
            rankstat_samples.append(item[1])

    return filename_rankstatheaders, rankstat_samples


def determine_all_rankstat_combinations(all_sets):
    """
    Create two lists necessary for unique filenames and to determine index + average of samples in following functions
    These are based on all possible combinations between rankstat columns, rather than individual rankstat columns.
    :param all_sets: Dictionary containing each rankstat column, containing all unique values with corresponding samples
    :return: A list containing the filenames (rankstat+unique_value) and a nested list containing their corresponding
    samples (unique sample combination, no duplicates)
    """
    unique_sample_combinations = []
    all_names = sorted(all_sets)
    combinations_rankstatheaders = list(itertools.product(*(all_sets[Name] for Name in all_names)))
    combinations_samples = list(itertools.product(*(all_sets[Name].values() for Name in all_names)))
    # Make a new list containing only unique samples for each combination
    for index, combination in enumerate(combinations_samples):
        unique_sample_combinations.append([])
        for allsamples in combination:
            for sample in allsamples:
                if sample not in unique_sample_combinations[index]:
                    unique_sample_combinations[index].append(sample)
    return combinations_rankstatheaders, unique_sample_combinations


def combination_to_index(sample_index, unique_sample_combinations):
    """
    Create a nested list with the index for each sample of each combination needed for determining the avarage in
    following functions
    :param sample_index: Dictionary containing the rankstat columns and their corresponding index.
    :param unique_sample_combinations: Nested list containing their corresponding (unique/no duplicate) samples
    :return: Nested list containing the index for each sample (for each unique combination)
    """
    indexed_combinations = []
    for index, combination in enumerate(unique_sample_combinations):
        indexed_combinations.append([])
        for unique_sample in combination:
            indexed_combinations[index].append(sample_index.get(unique_sample))
        indexed_combinations[index] = sorted(indexed_combinations[index])
    return indexed_combinations


def sample_average(combination, combination_header):
    """
    Calculates the average values of samples for each combination of rankstat column values.
    :param combination: Unique combination of rankstat column values
    :param combination_header: Unique filename for each combination
    :return: .csv file (used to create sankey for sample average)
    """
    try:
        average_samples = []
        rankstat_values = []
        with open(biotavizfile, "r") as file:
            for _ in range(1):  # skip column headers + root
                next(file)
            for line in file:
                if line.strip():
                    line = line.rstrip().split('\t')
                    for rankstat_column_index in combination:
                        rankstat_values.append(float(line[rankstat_column_index]))
                    average_samples.append(sum(rankstat_values) / len(rankstat_values))
                    rankstat_values = []
        hierarchy_counts(tax_filter, sample, average_samples, combination_header)
    except IndexError:
        print("IndexError; check if the correct file is given as input: ", traceback.print_exc())


def sample_average_all():
    """
    Calculates the average values of all samples, to be used in generating a .csv file. (even without value in rankstat)
    :return: .csv file (used to create sankey for sample average)
    """
    try:
        average_samples_all = []
        with open(biotavizfile, "r") as file:
            for _ in range(1):  # skip column headers + root
                next(file)
            for line in file:
                if line.strip():
                    line = line.rstrip().split('\t')
                    all_samples = []
                    for value in line[2:]:
                        all_samples.append(float(value))
                    average_samples_all.append(sum(all_samples) / len(all_samples))
        filename_part = "AverageAllSamples"
        hierarchy_counts(tax_filter, sample, average_samples_all, filename_part)
    except IndexError:
        print("IndexError; check if the correct file is given as input: ", traceback.print_exc())


if __name__ == "__main__":
    # Parameter for filtering on relative abundance
    tax_filter = 0.01
    # If true, generate a file for each sample column from input file
    sample_repeat = "false"
    # If true, generate a file for each unique combination of the rankstat
    combine_rankstat = "false"
    # Variables which need to be set beforehand
    average_all_samples = []
    sample = 0
    biotavizfile = "relative-table.biotaviz_relative_abundance.txt"
    mappingfile = "Metadata.tsv"
    filename_combination = ""
    try:
        if len(sys.argv) > 1:
            tax_filter = float(sys.argv[1])
            sample_repeat = str(sys.argv[2])
            mappingfile = str(sys.argv[3])
            combine_rankstat = str(sys.argv[4])
            biotavizfile = str(sys.argv[5])
            if tax_filter < 0 or tax_filter > 1:
                sys.exit(print("Use a number between 0 and 1 as parameter for filtering relative abundance"))
        main(tax_filter, sample_repeat, mappingfile, combine_rankstat)
    except ValueError:
        print("Parameter given was not a valid numeric value: ", traceback.print_exc())
        print("If the input is a decimal number, use a decimal point instead of comma (eg 0.01 instead of 0,01)")
    except IndexError:
        print("Not enough parameters were given: ", traceback.print_exc())
