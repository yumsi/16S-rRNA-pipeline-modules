#!/usr/bin/python3
# Jos Boekhorst
import sys


def load_txt(file):
    fileinput = open(file, 'rU')
    text = fileinput.read()
    fileinput.close()
    return text


# main program
if __name__ == '__main__':
    sys.stderr.write(
        "WARNING: relative abundance is relative to domain Bacteria, not all reads! The value no - Root is forced to 1. All non-bacterial taxa were REMOVED!\n")

    if not len(sys.argv) == 2:
        sys.stderr.write('Use: %s <Biotaviz_style_file>\n' % sys.argv[0])
        sys.exit(1)

    infile = sys.argv[1]
    outfile = infile.replace('.txt', '_relative_abundance.txt')
    if infile == outfile:
        sys.stderr.write("Outfile same as infile, aborting (infile should be .txt)\n")
        sys.exit(1)

    output = open(outfile, 'w')
    lines = load_txt(infile).strip().split('\n')
    output.write(lines[0] + '\n')

    root_i = None
    root_trace = None

    lines = [element.split('\t') for element in lines]
    for i, line in enumerate(lines):  # find the line with "domain_bacteria"
        if line[1] == "domain - Bacteria" or line[1] == "level0 - Root":  # last bit: BiotaViz PICRUSt - style file
            root_i = i
            root_trace = line[0]
            break

    if root_i is None:
        sys.stderr.write("Did not find Bacteria root taxon, aborting\n")
        sys.exit()

    root_values = [float(element) for element in lines[root_i][2:]]

    for line in lines[1:]:
        relative_abundances = []
        for i, count in enumerate(line[2:]):
            relative_abundances.append(min("1", ('%.20f' % (float(count) / root_values[i])).rstrip('0').rstrip('.')))

        if line[0][:len(root_trace)] == root_trace:
            output.write('\t'.join([line[0]] + [line[1]] + relative_abundances) + '\n')

    output.close()
    print('Done, data written to %s' % outfile)
