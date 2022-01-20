import pandas as pd
from optparse import OptionParser

def main():
    options, verbose = parse_args()
    df = open_file(options, verbose)
    filter_unclassified(df, options, verbose)

def parse_args():
    """
    Parse CLI options.

    Returns
    -------
    dict
        User CLI options

    """

    usage = "usage: %prog [options]"
    required = ["input_file", "output_prefix"]
    parser = OptionParser(usage)
    parser.add_option('-i', '--input', 
                      dest="input_file", help="Raw absolute abundance feature table with taxonomy for a specific taxonomic level.")
    parser.add_option('-o', '--output',
                      dest="output_prefix", default="", help="Prefix to standard output filename.")
    parser.add_option('-v', '--verbose',
                      dest="verbose",
                      default=False,
                      action="store_true"
                      )                  
    options, args = parser.parse_args()

    for r in required:
        if options.__dict__[r] is None:
            parser.print_help()
            parser.error("parameter '%s' required"%r)
            sys.exit(1)
        
    if options.verbose:
        print ('INPUT FILE     :', options.input_file)
        print ('OUTPUT PREFIX    :', options.output_prefix)
        print ('VERBOSE        :', options.verbose)
        verbose = options.verbose
    else:
        verbose = False
    return options, verbose

def open_file(options, verbose):
    """
    Opens an abundance table of a specific taxa level. 
    
    Skip first row.

    Parameters
    ----------
    options : dict
        User CLI options
    verbose : bool
        Print logs

    Returns
    -------
    Dataframe
        Formatted for downstream usage.
    """
    if verbose:
        print("Reading %s..." % options.input_file)
    df = pd.read_csv(options.input_file, sep="\t", skiprows=1)
    return df

def filter_unclassified(df, options, verbose):
    """
    Generates both absolute and relative abundance files using Bacteria as root. Records are removed if they are misplaced in a specific taxa level raw abundance file.

    Parameters
    ----------
    df : Dataframe
        Dataframe containing the taxa and raw counts at a specific taxa level
    output_file : str
        Name of the output file
    verbose : bool
        Print logs

    """

    if verbose:
        print("Creating absolute and relative abundance feature tables...")

    df = df[~df['#OTU ID'].str.contains(";__")]
    df.set_index('#OTU ID', inplace=True)
    abs_df = df
    rel_df = df.div(df.sum())

    if len(options.output_prefix) > 0:
        if options.output_prefix[-1] != "_":
            prefix = options.output_prefix + "_"
    else:
        prefix = options.output_prefix
    abs_df.to_csv(prefix + 'absolute_feature_table.tsv', sep='\t', index=True, header=True)
    rel_df.to_csv(prefix + 'relative_feature_table.tsv', sep='\t', index=True, header=True)

main()
