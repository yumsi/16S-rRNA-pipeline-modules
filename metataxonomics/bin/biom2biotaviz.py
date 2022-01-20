#!/usr/bin/python3
# todoc
"""
biom2biotaviz
-----------

.. module:: biom2biotaviz
  :synopsis: Convert biom file to BiotaViz-style txt file

Script generates a BiotaViz-style tab-delimited txt file from a biom file (v1).
Output is written to <infile>.txt and <infile>.biotaviz.txt

Typical run::

    biom2biotaviz.py some_biom_file.biom1
    
Changes:
28/06/2020: JB: String formating update, added 'd' to label_replace dictionary
15/01/2021: TE: Script can now handle "Unassigned" OTU/ASV when these are not properly filtered out before the biom is created, includes warning and error file

Author: Jos Boekhorst
"""
import sys
import os


def usage():
    sys.stderr.write("Use: {sys.argv[0]} -i <infile>\n")


def read_txt(filename):
    infile = open(filename, 'r')
    text = infile.read().rstrip('\n')
    infile.close()
    return text


def remove_empty_terminals(taxon):
    """remove terminal taxonomy bits that are empty"""
    new_taxon = []
    taxon_split = taxon.split('; ')
    for element in taxon_split:
        if element.split('__')[-1] != "":
            new_taxon.append(element)
        else:
            break
    return '; '.join(new_taxon)


def read_OTU_table(filename):
    # first line is "from biom" comment, second line is header
    unassigned_taxa = 0
    unassigned_read = 0
    taxonomy = {}
    counts = {}
    collapsed = {}
    lines = read_txt(filename)
    lines = lines.split('\n')
    samples = lines[1].split('\t')[1:-1]
    # print(samples)
    for line in lines[2:]:
        lineg = line.split('\t')
        if lineg[-1] == "Unassigned":
            # print("Warning ! Unassigned Taxa !")
            unassigned_taxa += 1
            # print( type(lineg[2]) )
            inters = [ int(float(element)) for element in lineg[1:-1] ]
            # print( lineg[1:-1] )
            # print( inters )
            # print( sum(inters) )
            unassigned_read += sum(inters)
            # print("# Total sum of unassigned reads:", unassigned_read)
            # for i in lineg[1:-1]:
                # print( i, float(i), int(float(i)) )
            continue
        OTU = lineg[0]
        taxon = 'r__Root; '+lineg[-1]
        taxon = taxon.replace("NA;", "k__;")  # NG_Tax weirdness
        taxon = remove_empty_terminals(taxon)
        counts[OTU] = {}
        taxonomy[OTU] = taxon
        if taxon not in collapsed:
            collapsed[taxon] = {}
        for i, sample in enumerate(samples):
            count = float(lineg[i + 1])
            counts[OTU][sample] = count
            if sample not in collapsed[taxon]:
                collapsed[taxon][sample] = 0
            collapsed[taxon][sample] += count
    return collapsed, taxonomy, samples, unassigned_taxa, unassigned_read


def get_new_number(taxon_to_trace, trace):
    """get a new number for a taxon trace"""
    parent = trace[:-1]
    # get current highest child
    current_children = [0]
    for element in taxon_to_trace.keys():
        element_parent = element[:-1]
        if element_parent == parent:
            current_children.append(int(taxon_to_trace[element].rstrip('.').split('.')[-1]))
    return max(current_children) + 1


def traces_from_taxonomy(collapsed):
    taxon_to_trace = {tuple(['r__Root']): '1.'}
    for taxon in collapsed:
        taxon_gs = taxon.split('; ')
        full_name = []
        for taxon_g in taxon_gs:
            label = taxon_g
            full_name.append(label)

        # Now get a code. Parents need to be inferred.
        for i in range(1, len(full_name)):
            tmp = full_name[:i + 1]
            parent = full_name[:i]
            if not tuple(tmp) in taxon_to_trace:
                new_nr = get_new_number(taxon_to_trace, tuple(tmp))
                new_trace36 = f"{taxon_to_trace[tuple(parent)]}{str(new_nr)}."  # Python 3.6.12
                # print("Python 3.6:", new_trace36)
                new_trace35 = taxon_to_trace[tuple(parent)] + str(new_nr) + "." # Python 3.5.2
                # print("Python 3.5:", new_trace35)
                taxon_to_trace[tuple(tmp)] = new_trace36
    return taxon_to_trace


def reverse_dictionary(input_dict):
    # generate new dictionary by swapping key and value
    new_dict = {}
    for element in input_dict.keys():
        element_value = input_dict[element]
        new_dict[element_value] = element
    return new_dict


def infer_internal_counts(collapsed):
    """data was collapsed per taxon, now infer parent counts"""
    new_counts = {}
    for taxon in collapsed.keys():
        taxon_g = taxon.split('; ')
        for i in range(len(taxon_g)):
            taxon_sub = tuple(taxon_g[:i+1])
            if taxon_sub not in new_counts.keys():
                new_counts[taxon_sub] = {}
            for sample in collapsed[taxon].keys():
                if sample not in new_counts[taxon_sub].keys():
                    new_counts[taxon_sub][sample] = 0
                new_counts[taxon_sub][sample] += collapsed[taxon][sample]
    return new_counts


# settings
convert_command = "biom convert -i #infile# -o #infile#.txt --to-tsv --header-key taxonomy"

label_replace = {'r': 'no',
                 'k': 'domain',
                 'd': 'domain',
                 'p': 'phylum',
                 'c': 'class',
                 'o': 'order',
                 'f': 'family',
                 'g': 'genus',
                 's': 'species',
                 't': 'variant'}

if __name__ == "__main__":
    try:
        infile = sys.argv[1]
    except IndexError:
        usage()
        sys.exit()
    outfile = infile + '.txt'
    if os.path.isfile(outfile):
        sys.stderr.write("Outfile {outfile} exists, aborting\n".format(outfile=outfile))
        sys.exit()

    # generate tab-delimited OTU matrix with taxonomy in final column
    # alterantive would be the Python biom functions
    command = convert_command.replace('#infile#', infile)
    sys.stderr.write("# Executing: {command}\n".format(command=command))
    os.system(command)


    # read the tab-delimited data & collapse
    sys.stderr.write("# Reading input data from the biom.txt ..\n")
    collapsed, taxonomy, samples, unassigned_taxa, unassigned_read = read_OTU_table(outfile)
    taxon_to_trace = traces_from_taxonomy(collapsed)
    trace_to_taxon = reverse_dictionary(taxon_to_trace)
    counts = infer_internal_counts(collapsed)

    # printing the results
    if unassigned_taxa > 0:
        warning = f"# Warning: The dataset contains {unassigned_taxa} unclassified taxa with a total number of {unassigned_read} reads !\n# These unassigned OTU/ASV have been removed, but please have a look at them in the original biom or biom.txt file \n"
        sys.stderr.write( warning )
        with open('error.log', 'w') as f:
            f.write(warning)
            f.close()
    sys.stderr.write("# Printing output to BiotaViz ..\n")
    traces = list(trace_to_taxon.keys())
    traces.sort()
    samples.sort()
    print("#class\tclass id\t"+"\t".join(samples))

    for trace in traces:
        line = [trace]
        nice_taxon = trace_to_taxon[trace][-1]
        nice_taxon = nice_taxon.split('__')
        nice_taxon = label_replace[nice_taxon[0]] + " - " + nice_taxon[-1]
        line.append(nice_taxon)
        taxon = trace_to_taxon[trace]
        for sample in samples:
            line.append("%f" % (counts[taxon][sample]))
        print("\t".join(line))
