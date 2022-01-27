# Author: Harm Laurense
# Last edited: 18-10-2021
# Function: This script is for creating a biotaviz test file to test max size on 1% filter (e.g. 100 times 1% species)

import csv

taxonomic_ranks_dict_perc = {
    "domain": 1,
    "phylum": 1,
    "class": 1,
    "order": 1,
    "family": 1,
    "genus": 1,
    "species": 0.01}


def main():
    filename = "relative-table_kopie_test.biotaviz_test.txt"
    lines = load_txt(filename).strip().split('\n')
    write_testfile(lines)


def load_txt(filename):
    fileinput = open(filename)
    text = fileinput.read()
    fileinput.close()
    return text


def write_testfile(lines):
    tax_ranks_skip = []
    with open("biotaviz_maxsize_test.txt", "w", newline="") as f:
        firstline = 0
        wr = csv.writer(f, delimiter='\t')
        for line in lines:
            if firstline < 2:
                line = line.split("\t")
                wr.writerow(line)
                firstline += 1
            else:
                line = line.split("\t")
                tax_rank_1sthalf = line[1].split('-')[0].rstrip()
                tax_rank_2ndhalf = line[1].split('-')[1].rstrip()
                for index, element in enumerate(line[2:]):
                    if tax_rank_1sthalf not in tax_ranks_skip:
                        line[index + 2] = taxonomic_ranks_dict_perc.get(tax_rank_1sthalf)
                        if index == 3:
                            wr.writerow(line)
                            tax_ranks_skip.append(tax_rank_1sthalf)
                    else:
                        tax_rank_1sthalf = "species"
                        line[1] = "species - " + str(tax_rank_2ndhalf)
                        line[index + 2] = taxonomic_ranks_dict_perc.get(tax_rank_1sthalf)
                        if index == 3:
                            wr.writerow(line)
                            tax_ranks_skip.append(tax_rank_1sthalf)


main()
