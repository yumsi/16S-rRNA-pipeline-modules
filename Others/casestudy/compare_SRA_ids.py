# Author: Harm Laurense
# Last edited: 23-11-2021
# Function: Retrieve SRA ids by comparing to the metadata ids

def main():
    metadata_file = "NIHMS569508-supplement-02.txt"
    sra_file = "SraRunTable.txt"
    metadata_ids_list = get_metadata_ids(metadata_file)
    sra_dict = get_sra_ids(sra_file)
    overlapping_ids= compare_id_lists(metadata_ids_list, sra_dict)
    write_new_metadata_file(metadata_file, overlapping_ids, sra_dict)


def read_txt(filename):
    infile = open(filename, 'r')
    headers = infile.readline()
    text = infile.readlines()
    infile.close()
    return text, headers


def get_metadata_ids(metadata_file):
    id_list = []
    content, headers = read_txt(metadata_file)
    for line in content:
        line = line.split("\t")
        id_list.append(line[0].rstrip())
    return id_list


def get_sra_ids(sra_file):
    content2, headers = read_txt(sra_file)
    count = 0
    id_dict = {}
    chosen_sample_regions = ["stool"]
    for line in content2:
        count += 1
        line = line.split("\t")
        # 42 is column AQ (sample name) and 52 is for column AZ (region) ; region=area sample is taken from
        if line[42] and line[51]:
            sample_region = line[51]
            if sample_region in chosen_sample_regions:
                sample_name = line[42]
                if "." in sample_name:
                    sample_name = sample_name.split(".")[1]
                if "SKBTI-" in sample_name:
                    sample_name = sample_name.replace("-", "")
                if "_S" in sample_name:
                    sample_name = sample_name.replace("_", ".")
                if id_dict.get(sample_name):
                    id_dict[sample_name.rstrip()].append(line[0].rstrip())
                else:
                    id_dict[sample_name.rstrip()] = [line[0].rstrip()]
    print(id_dict)
    print(len(id_dict))
    return id_dict


def compare_id_lists(metadata_ids_list, dict):
    overlapping_metadata_ids = []
    overlapping_sra_ids = []
    no_matches = []
    merged_ids_list = []
    # print(metadata_ids_list)
    with open('sra_overlapping_ids.txt', 'w') as f:
        with open('sample_overlapping_ids.txt', 'w') as f2:
            for id in metadata_ids_list:
                if dict.get(id):
                    for sra in dict.get(id):
                        overlapping_sra_ids.append(sra)
                        f.write(sra + "\n")
                        if id not in overlapping_metadata_ids:
                            overlapping_metadata_ids.append(id)
                            f2.write(id + "\n")
                else:
                    no_matches.append(id)

    f.close()
    f2.close()
    # overlapping_ids = sorted(overlapping_ids)
    # no_matches = sorted(no_matches)
    print(overlapping_sra_ids)
    print(len(overlapping_sra_ids))
    print(overlapping_metadata_ids)
    print(len(overlapping_metadata_ids))
    # print(no_matches)
    print(len(no_matches))
    return overlapping_metadata_ids


def write_new_metadata_file(original_metadata_file, overlapping_metadata_ids, dict):
    content3, headers = read_txt(original_metadata_file)
    with open('metadata_overlapping.txt', 'w') as new_metadata_file:
        new_metadata_file.write("sra_samples\t"+headers)
        for line in content3:
            line_list = line.split("\t")
            metadata_id = line_list[0]
            if line_list[0] in overlapping_metadata_ids and len(dict.get(metadata_id)) == 1:
                new_metadata_file.write(dict.get(metadata_id)[0]+"\t"+str(line))
    new_metadata_file.close()


main()
