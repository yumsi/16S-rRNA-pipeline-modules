# QIIME 2 16S Microbiomics Pipeline

  

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A518.10.1-brightgreen.svg)](https://www.nextflow.io/)

  

## Requirements

  

*  [Nextflow](https://www.nextflow.io/index.html#GetStarted) (>= 18.10.1)

*  [Docker](https://docs.docker.com/engine/install/)

  

## Usage

### Clone repository
Clone the repository:<br>
`git clone https://gitlab.cmbi.umcn.nl/rtc-bioinformatics/metataxonomics.git`

Navigate to the directory:<br>
`cd metataxonomics`

### Download dependencies
Download the dependencies:<br>
`bash bin/download_dependencies.sh`

### Test pipeline
A test run of the pipeline can be performed by executing the following command:

  

`nextflow run main.nf -profile docker,test --reference_reads_for_taxonomic_classifier lib/SILVA-138-SSURef-Full-Seqs.qza --reference_taxonomy_for_taxonomic_classifier lib/Silva-v138-full-length-seq-taxonomy.qza -resume`


## Parameters

### `-profile <str>`

Indicate which profiles needs to be used. If selecting multiple profiles, seperate the trings by a comma.

Profiles to select from:
<ul>
  <li>docker      --> use Docker to run pipeline</li>
  <li>conda       --> use Conda to run pipeline</li>
  <li>singularity --> use Singularity to run pipeline</li>
  <li>test        --> use the test configurations</li>
  <li>debug       --> print hostname for debugging purposes</li>
</ul>

### `--extension <str>`

The expected .fastq file extensions. <br>

For example: <br>
`--extension "/*_R{1,2}*.{fastq.gz,fastq}"`

  
### `--reads <str>` 
Path to directory containing your reads.<br>
For example: <br>
Below an example directory containing test reads:
```
test_data
  │   test_R1.fastq
  │   test_R2.fastq
```
`--reads /test_data`

### `--FW_adapter <str>`
Forward adapters/primers to be trimmed. 

### `--RV_adapter <str>`
Reverse adapters/primers to be trimmed. 

### `--outdir <str>`<br>
Directory where results are saved.
For example: <br>
`--outdir "./results"`
### `--keep_untrimmed <bool>`
Indicate whether reads that were not able to be trimmed are kept.

### `--adapters <str>`
Path to file containing universal adapters for a global trimming of reads. `--adapters "/res/adapters.fasta"`
### `--classifier <str>`
Path to .qza file containing a trained classifier.


### `--train_classifier <bool>` 
Indicate whether a new classifier should be trained with the files provided with the `--reference_reads_for_taxonomic_classifier` and `--reference_taxonomy_for_taxonomic_classifier` parameters.

### `--use_taxonomy_classifier <bool>`
Indicate whether to use the classifier for the taxonomic classification. VSEARCH is the default classification method.

### `--reference_reads_for_taxonomic_classifier <str>`
Path to a .qza reference <i>reads</i> file for VSEARCH and the training of a new classifier (if applicable).

### `--reference_taxonomy_for_taxonomic_classifier <str>`
Path to a .qza reference <i>taxonomy</i> file for VSEARCH and the training of a new classifier (if applicable).

### `--taxa_level_for_itol <int>`
Indicate which taxonomic level (1 through 7) should be used for the iTOL processes.
### `--condition_header <str>`
Column name in the metadata file used to specify the conditions for the various samples.
### `--skip_qc <bool>`
Skip the quality control steps (FastQC, MultiQC).
### `--skip_tree <bool>`
Skip the tree building process.
### `--skip_barplot <bool>`
Skip the QIIME2 barplot creation process.
### `--skip_diversity <bool>`
Skip all processes that calculate diversity metrics.

## Flowchart
[Click here for flowchart](https://app.creately.com/diagram/S6a55prj4I1/view "Flowchart of pipeline")
