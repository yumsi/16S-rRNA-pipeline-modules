#!/usr/bin/python3
"""
JOS_mapping_rankstat
---------------
.. module:: JOS_mapping_rankstat
  :synopsis: do rankstats from mapping and biotiviz files
.. moduleauthor:: Jos Boekhorst

Perform rankstats from mapping and biotaviz-style input files

Typical run::

JOS_mapping_rankstat -m mapping.txt -b biotaviz.txt

Run the script with '-h' for a list of options.
Low-effort and more transparent alternative for LEfSe

"""

import sys
import argparse
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from math import log


def load_txt(filename):
    fileinput = open(filename)
    text = fileinput.read()
    fileinput.close()
    return text


def get_sets(filename):
    """parse mapping file: what subsets should be compared"""
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


def get_counts(counts_file, comment_columns=(0,), taxon_column=1):
    """return dictionary counts[taxon][sample] = count.
    defaults are for biotaviz style, set comment_columns to [] and taxon_column to 0 for simple table"""
    if comment_columns is None:
        comment_columns = [0]
    lines = load_txt(counts_file).rstrip('\n').split('\n')
    header = lines[0].split('\t')
    samples = []
    counts = {}
    for i, element in enumerate(header):
        if i not in comment_columns and i != taxon_column:
            samples.append(element)
    taxa = []
    postfixes = {}
    for line in lines[1:]:
        lineg = line.split('\t')
        taxon = lineg[taxon_column]
        if taxon in taxa:
            sys.stderr.write(f"WARNING: duplicate taxon {taxon}, adding postfix\n")
            if not taxon in postfixes:
                postfixes[taxon] = 0
            postfixes[taxon] += 1
            taxon = taxon + '_' + str(postfixes[taxon])
        taxa.append(taxon)
        if taxon not in counts:
            counts[taxon] = {}
        relevant_i = 0  # count index
        for i, element in enumerate(lineg):
            if i not in comment_columns and i != taxon_column:
                sample = samples[relevant_i]
                counts[taxon][sample] = float(element)
                relevant_i += 1
    return counts, taxa


def do_kruskal(dataset, taxon, counts):
    values = []
    for subset in dataset.values():
        subvalues = []
        for sample in subset:
            subvalues.append(counts[taxon][sample])
        values.append(subvalues)
    try:
        kruskal_results = kruskal(*values)
    except ValueError:  # happens when all values are identical
        kruskal_results = ["NA", "NA"]
    return kruskal_results[0], kruskal_results[1]


def do_mannwhitneyu(values_A, values_B):
    try:
        MWU_results = mannwhitneyu(values_A, values_B)
    except ValueError:  # happens when all values are identical
        MWU_results = ["NA", "NA"]
    return MWU_results[0], MWU_results[1]


############
# settings #
############


description = "Do rank statistics from mapping and biotaviz-style input files"
max_factorial = 1

# main program
if __name__ == '__main__':
    # argument parsing
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument('-b', dest='biotavizfile', help='biotaviz input file', required=True)
    parser.add_argument('-m', dest='mappingfile', help='mapping input file', required=True)
    parser.add_argument('-o', dest='outfile', help='output file', required=True)
    options = vars(parser.parse_args())
    mappingfile = options['mappingfile']
    countsfile = options['biotavizfile']
    outfile = options['outfile']

    # get comparisons and samples per class from mapping file
    all_sets = get_sets(mappingfile)
    counts, taxa = get_counts(countsfile)

    header = ['comparison', 'taxon', 'classA', 'classB', 'kruskal_H', 'kruskal_p', 'MWU_U', 'MWU_p', 'A_average',
              'B_average',
              'fold_increase']  # , 'A_samples', 'B_samples', 'A_values', ' B_values'] # comment: add actual values

    output = open(outfile, 'w')
    output.write("\t".join(header) + '\n')

    for setname in list(all_sets.keys()):  # loop over datasets ("columns") from mapping data
        for taxon in taxa:
            # noinspection PyDictCreation
            data = {}
            data['comparison'] = setname
            data['taxon'] = taxon
            data['kruskal_H'], data['kruskal_p'] = do_kruskal(all_sets[setname], taxon, counts)
            # pairwise comparisons
            subsets = list(all_sets[setname].keys())
            for iA in range(len(subsets)):
                values_A = [counts[taxon][sample] for sample in all_sets[setname][subsets[iA]]]
                for iB in range(iA + 1, len(subsets)):
                    values_B = [counts[taxon][sample] for sample in all_sets[setname][subsets[iB]]]
                    data['classA'] = subsets[iA]
                    data['classB'] = subsets[iB]
                    data['MWU_U'], data['MWU_p'] = do_mannwhitneyu(values_A, values_B)
                    data['A_average'] = sum(values_A) / len(values_A)
                    data['B_average'] = sum(values_B) / len(values_B)
                    try:
                        data['fold_increase'] = log(data['B_average'], 2) - log(data['A_average'], 2)
                    except ValueError:  # zero values in average abundance
                        data['fold_increase'] = "NA"
                    output.write('\t'.join([f"{data[element]}" for element in header]) + '\n')
    output.close()
