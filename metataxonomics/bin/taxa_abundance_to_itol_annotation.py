import pandas as pd
from optparse import OptionParser
import sys
import random
import re

def main():
	options, verbose = parse_args()
	df = open_file(options, verbose)
	meta_df = extract_treatment_metadata(options, verbose)
	create_itol_label_file(meta_df, options.output_filename, verbose, options.condition_header)
	create_itol_annotation_file(df, options.output_filename, verbose)
	
def parse_args():
	"""
    Parse CLI options.

    Returns
    -------
    dict
        User CLI options

    """

	usage = "usage: %prog [options]"
	required=["input_file", "output_filename"]
	parser = OptionParser(usage)
	parser.add_option('-i', '--input', 
		              dest="input_file", help="Filtered relative abundance feature table with taxonomy of a specific taxonomic level.")
	parser.add_option('--metadata', 
		              dest="metadata", help="Meatadata with column 'treatment' for every sample.", default="")
	parser.add_option('-o', '--output',
		              dest="output_filename", help="iTOL annotation file for visualization.")
	parser.add_option('--condition_header',
		              dest="condition_header",
		              default="treatment",
		              help="Header name in metadata file containing sample conditions to be used within the iTOL plot."
		              )     
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
		print ('OUTPUT FILE    :', options.output_filename)
		print ('CONDITION_HEADER    :', options.condition_header)
		print ('VERBOSE        :', options.verbose)
		verbose = options.verbose
	else:
		verbose = False
	return options, verbose


def open_file(options, verbose=False):
	"""
    Opens a user submitted TSV file.

    Parameters
    ----------
    options : dict
        Provided CLI arguments
    verbose : bool
        Print logs

    Returns
    -------
    Dataframe
        Dataframe with input file contents

    """

	if verbose:
		print("Reading %s..." % options.input_file)

	df = pd.read_csv(options.input_file, sep="\t") # TE : Warning : If a column has only integer values, this leads to an error in string joining at line 156 or 164 !
	df.set_index('#OTU ID', inplace=True)
	print ('\n#here!\n')
	print (df.T)
	return df.T

def extract_treatment_metadata(options, verbose=False):
	"""
    Transposes the metadata dataframe and sets the sample ID as index.
	
	Parameters
    ----------
    options : dict
        Provided CLI arguments
    verbose : bool
        Print logs
    
	Returns
    -------
    DataFrame
        Transposed Metadata Dataframe

    """

	if verbose:
		print("Opening & extracting treatment column from %s..." % options.metadata)

	df = pd.read_csv(options.metadata, sep="\t")
	df.set_index('ID', inplace=True)
	return df.T

def generate_random_hex():
	"""
    Generates a random hex color code for all the taxa available.

    Returns
    -------
    str
        Hex color code

    """

	return "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

def create_itol_label_file(df, output_file, verbose=False, metadata_condition_header="treatment"):
	"""
    Creates a label file to be used with the iTOL rooted tree and iTOL annotation file.

    Parameters
    ----------
    df : Dataframe
        Metadata dataframe containing the sample ID and conditions
    output_file : str
        Name of the output file
    verbose : bool
        Print logs
	metadata_condition_header: str
		Header name containing the sample conditions
    """
	meta = 'DATASET_STYLE'
	separator = 'SEPARATOR COMMA' 
	dataset_label = "DATASET_LABEL,conditions"
	condition_metadata_colors = []
	condition_color_dict = {}
	sample_condition_dict = {}
	sample_list = []
	for index in df.head().index:
		if index == metadata_condition_header:
			conditions = df.loc[metadata_condition_header]
			for sample in conditions.index:
				sample_list.append(sample)
				if sample not in sample_condition_dict:
					sample_condition_dict[sample] = conditions[sample]
				if conditions[sample] not in condition_color_dict: 
					condition_color_dict[conditions[sample]] = generate_random_hex()
				temp_data = '\n'.join([','.join([sample, ','.join(["label", "node", "#000000","1", "normal",condition_color_dict[sample_condition_dict[sample]]])])])
				condition_metadata_colors.append(temp_data)
	sample_colors = [condition_color_dict[sample_condition_dict[sample]] for sample in sample_list]
	legend_title = ','.join(['LEGEND_TITLE', 'Sample Conditions'])
	legend_position_x = ','.join(['LEGEND_POSITION_X', '350'])
	legend_position_y = ','.join(['LEGEND_POSITION_Y', '100'])
	legend_shapes = ','.join(['LEGEND_SHAPES', ','.join([str(1)]*len(condition_color_dict.keys()))])
	legend_colors = ','.join(['LEGEND_COLORS', ','.join(condition_color_dict.values())])
	legend_labels = ','.join(['LEGEND_LABELS', ','.join(condition_color_dict.keys())])
	# print(condition_metadata_colors)
	data = '\n'.join(['DATA', '\n'.join([str(line) for line in condition_metadata_colors])])
	final_conditions = '\n\n'.join([meta, separator, dataset_label, legend_title, legend_position_x, legend_position_y, legend_shapes,legend_colors, legend_labels, data])
	
	if verbose:
		print("Exporting label file as %s..." % ("conditions_" + output_file))

	try:
		with open("conditions_" + output_file, 'w') as out_file:
			out_file.write(final_conditions)
			out_file.close()
	except Exception:
		print("ERROR: One or more files could not be created!")
	


def create_itol_annotation_file(df, output_file, verbose=False):
	"""
    Creates annotation file to be used with the iTOL rooted tree.

    The annotation includes a stacked barplot and a legend.

    Parameters
    ----------
    df : Dataframe
        Dataframe containing the taxa and counts
    output_file : str
        Name of the output file
    verbose : bool
        Print logs

    """
    
	if verbose:
	    print("Generating iTOL annotation file...")
	meta = 'DATASET_MULTIBAR'
	label = ','.join(['DATASET_LABEL','sample genus stacked barplot'])
	separator = 'SEPARATOR COMMA' 
	df_column_names = df.columns.values
	field_labels = ','.join(['FIELD_LABELS', ','.join(df_column_names)])
	if verbose:
		print("Generating hex color codes...")
	colors = [generate_random_hex() for i in range(len(df_column_names))]
	dataset_scale = ','.join(['DATASET_SCALE', ','.join(str(i) for i in range(0, 110, 10))])
	field_colors = ','.join(['FIELD_COLORS', ','.join(colors)])
	legend_title = ','.join(['LEGEND_TITLE', 'Legend'])
	legend_shapes = ','.join(['LEGEND_SHAPES', ','.join([str(1)]*len(df_column_names))])
	legend_colors = ','.join(['LEGEND_COLORS', ','.join(colors)])
	legend_labels = ','.join(['LEGEND_LABELS', ','.join([re.split(r';[a-z]__', item)[-1] for item in df_column_names])])
	
	data_list = []
	# for index in df.head().index: # TE : I removed .head() because it resulted in printing only the first 5 samples, was probably used for testing ..
	for index in df.index:
		values = df.loc[index]*100
		temp_data = ','.join([index, ','.join(values.astype('str').values)])
		data_list.append(temp_data)
	data = '\n'.join(['DATA', '\n'.join([str(line) for line in data_list])])
	
	
	final = '\n\n'.join([meta, separator, label, dataset_scale, field_colors, field_labels, legend_title, legend_shapes, legend_labels, legend_colors, data])
	
	if verbose:
		print("Exporting annotation file as %s..." % output_file)
	
	with open(output_file, 'w') as out_file:
		out_file.write(final)
		out_file.close()
    
	if verbose:
		print("Done.")
	
main()
