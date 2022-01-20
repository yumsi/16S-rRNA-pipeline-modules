#!/usr/bin/env nextflow

/*
 * Main script for the microbiomics pipeline v1.0.
 *
 * Center for Molecular and Biomolecular Informatics (CMBI), Radboudumc, Nijmegen.
 *
 */

log.info "Microbiomics - N F"
log.info "====================================="
log.info "job id                 : $workflow.runName"
log.info "reads                  : ${params.reads}"
log.info "output                 : ${params.outdir}"
log.info "metadata file          : ${params.metadata}"
log.info "adapters               : ${params.adapters}"
log.info "classifer              : ${params.classifier}"
log.info "font arial             : ${params.font}"
log.info "====================================="
log.info "\n"

/*
 * Import input files (reads, metadata)
 */
if(params.readPaths && params.reads == "<data>${params.extension}"){
	// This channel is only made when the test profile is selected with "-profile test" on the CLI.
    Channel
        .from(params.readPaths)
        .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
        .ifEmpty { exit 1, "No read files supplied!" }
        .into { read_pairs; read_pairs_fastqc }
} else {
    Channel
        .fromFilePairs( params.reads + params.extension, size: 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}${params.extension}\nNB: Path needs to be enclosed in quotes!" }
        .into { read_pairs; read_pairs_fastqc }
}
if (params.metadata) {
	Channel.fromPath("${params.metadata}", checkIfExists: true)
	.into { ch_metadata_for_barplot; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_beta_diversity; metadata_for_merge; metadata_for_itol; metadata_for_rankstat1; metadata_for_rankstat2}
} else {
	Channel.from()
		.into { ch_metadata_for_barplot; ch_metadata_for_diversity_core; ch_metadata_for_alpha_diversity; ch_metadata_for_beta_diversity; metadata_for_merge; metadata_for_itol}
}

if (params.font) {
	Channel.fromPath("${params.font}", checkIfExists: true)
	.set { font_for_sankeypng}
	}

if((params.FW_adapter) && (params.RV_adapter)){
    params.adapters_provided = true
}
else{
    params.adapters_provided = false
}

/*
 * Import taxonomy classifier .qza files
 */
seqs_train = Channel.empty()
taxon_train = Channel.empty()
if (params.train_classifier && params.use_taxonomy_classifier) {
    try { 
        Channel.fromPath("${params.classifier}", checkIfExists: true)
                .set { qiime_classifier }
    } catch(Exception ex) {
        qiime_classifier = Channel.of("empty") 
    }
}
if (params.train_classifier) {
    if (params.classifier && params.use_taxonomy_classifier) {
        Channel.empty()
            .ifEmpty { exit 1, "Attempting to train a taxonomic classifier while a classifier has already been provided. Do not attempt to use --train_classifier with --classifier." }
    }
    if (params.reference_reads_for_taxonomic_classifier && params.reference_taxonomy_for_taxonomic_classifier) {
        Channel.fromPath("${params.reference_reads_for_taxonomic_classifier}", checkIfExists: true).set{ seqs_train }
        Channel.fromPath("${params.reference_taxonomy_for_taxonomic_classifier}", checkIfExists: true).set{ taxon_train }
    }else {
        seqs_train = Channel.empty()
        taxon_train = Channel.empty()
        Channel.empty()
            .ifEmpty { exit 1, "No reference .qza files or .qza classifier file provided for the training of the classifier. e.g. --reference_reads_for_taxonomic_classifier <path to file> --reference_taxonomy_for_taxonomic_classifier <path to file>" }
    }
}

if (params.classifier && params.use_taxonomy_classifier) {
	Channel.fromPath("${params.classifier}", checkIfExists: true)
		   .set { qiime_classifier }
} else if (params.use_taxonomy_classifier == false) {
    if (params.reference_reads_for_taxonomic_classifier && params.reference_taxonomy_for_taxonomic_classifier) {
        Channel.fromPath("${params.reference_reads_for_taxonomic_classifier}", checkIfExists: true).set{ seqs }
        Channel.fromPath("${params.reference_taxonomy_for_taxonomic_classifier}", checkIfExists: true).set{ taxon }
         qiime_classifier = Channel.of("empty")
    } else {
    Channel.empty()
        .ifEmpty { exit 1, "No reference .qza files or .qza classifier file provided for the taxonomic classification process. e.g. --reference_reads_for_taxonomic_classifier <path to file> --reference_taxonomy_for_taxonomic_classifier <path to file>" }
    }   
}

Channel
    params.makesankeyhtml = "${params.outdir}/sankey/preparation-files/*.csv"
    sankey_prep_datasets = Channel.fromPath(params.makesankeyhtml)
    params.makesankeypng = "${params.outdir}/sankey/html-files/*.html"
    sankey_html_datasets = Channel.fromPath(params.makesankeypng)

/*
 * Cutadapt: Trim non-biological sequences (primers, adapters) from 5' side if needed.
 *
 * Input: Paired-end reads.
 * Output: Trimmed paired-end reads.
 *         Trimming statistics.
 */
process trimming {
    tag "${pair_id}"
    
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: 'true',
        saveAs: {filename -> 
        if (filename.indexOf(".gz") == -1) "logs/$filename"
        else if(true==true) filename 
        else null}

    input:
    tuple val(pair_id), file(reads) from read_pairs

    output:
    tuple val(pair_id), file("trimmed/*.*") into trimmed_read_pairs_fastqc
    file "trimmed/*.*" into fastq_trimmed_for_qiime_manifest
    file "cutadapt_log_*.txt" into fastq_cutadapt_log

    script:
    def adapters = params.adapters_provided ? "-g ${params.FW_adapter} -G ${params.RV_adapter}" : "-a file:${params.adapters} -A file:${params.adapters}"
    def discard_untrimmed = params.keep_untrimmed ? '' : '--discard-untrimmed'
    def fwd_out = reads[0].toString().split("\\.", 2)[0] + "-trimmed." + reads[0].toString().split("\\.", 2)[1]
    def rev_out = reads[1].toString().split("\\.", 2)[0] + "-trimmed." + reads[1].toString().split("\\.", 2)[1]
    """
    mkdir -p trimmed
    cutadapt ${adapters} ${discard_untrimmed} \
        -o trimmed/${fwd_out} -p trimmed/${rev_out} \
        ${reads[0]} ${reads[1]} > cutadapt_log_${pair_id}.txt
    """
}

// Join raw (untrimmed) reads and trimmed reads into one channel for FastQC. Format: [<sample_id>, [<untrimmed_files*n>, <trimmed_files*n>]]
before_and_after_trim_reads_fastqc = read_pairs_fastqc.join(trimmed_read_pairs_fastqc)
                                        .map { id, untrimmed_read_list, trimmed_read_list -> [id, [untrimmed_read_list, trimmed_read_list].flatten()] }

/*
 * FastQC; quality control analysis of sequencing data
 *
 * Input: List with the sample ID and the associated paths to both untrimmed and trimmed reads. 
 *        Example: [<sample_id>, [<untrimmed_files>*n, <trimmed_files>*n]] where n is the amount of files.
 * Output: .html (report) and .zip (figures) for every individual read file provided.
 */
process fastqc {
    tag "${pair_id}"
    publishDir "${params.outdir}/fastQC", mode: 'copy', overwrite: 'true',
    saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    tuple val(pair_id), file(reads) from before_and_after_trim_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    when:
    !params.skip_qc
    
    script: 
    """
    fastqc -d . -o . -quiet ${reads}
    """
}

/*
 * MultiQC; aggregate results and/or logs of fastQC and Cutadapt of all samples into one HTML report.
 *
 * Input: Two lists with all the results and/or logs from the fastQC and Cutadapt processes. MultiQC will only be executed when the before-mentioned processes are completed. 
 *        Example: [<output>*n] where n is the amount of files.
 * Output: .html (report) and log files for every individual read file provided.
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy', overwrite: 'true'

    input:
    file("fastqc/*") from fastqc_results.collect()
    file("cutadapt/logs/*") from fastq_cutadapt_log.collect()

    output:
    file("*multiqc_report.html") into multiqc_report
    file("*_data")
    file("multiqc_data/multiqc_general_stats.txt")

    when:
    !params.skip_qc

    script:
    """
    multiqc --force --interactive .
    """
}

/*
 * Create a tab-seperated QIIME manifest file containing per sample the paths for the forward and reverse reads.
 * Example: 
 * """
 * sample-id	forward-absolute-filepath	reverse-absolute-filepath
 * SRR8698828	SRR8698828_R1-trimmed.fastq.gz	SRR8698828_R2-trimmed.fastq.gz 
 * """
 */
fastq_trimmed_for_qiime_manifest
            .map { forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] } 
            .map { name, forward, reverse -> [ name.toString().take(name.toString().indexOf("_")), forward, reverse ] } 
            .map { name, forward, reverse -> [ name + "\t" + forward + "\t" + reverse ] } 
            .flatten()
            .collectFile(name: 'manifest.tsv', newLine: true, storeDir: "${params.outdir}/demux", seed: "sample-id\tforward-absolute-filepath\treverse-absolute-filepath")
            .set { qiime_manifest }

/*
 * Import trimmed reads into QIIME2 artifact (=object)
 *
 * Input: .tsv manifest file (see above for format)
 * Output: .qza (QIIME2 artifact), .qzv (QIIME2 visualization, terminal output), and metrics about the sequences (e.g. counts per sample, quality). 
 *         The QIIME2 artifact is to be used in further QIIME2 analyses.
 */ 
process qiime_import {
    publishDir "${params.outdir}/demux", mode: 'copy', overwrite: 'true'
    
    input:
    file(manifest) from qiime_manifest
    
    output:
    file("*.qza") into trimmed_artifacts
    file("*.{qzv, html, csv}")
    file("demux_resources/*")
    
    script:
    """
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${manifest} \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path demux.qza
    
    qiime demux summarize \
        --i-data demux.qza \
        --o-visualization demux.qzv
    
    qiime tools export \
        --input-path demux.qzv \
        --output-path demux_resources
    """
}

/*
 * DADA2 denoising of paired-end reads.
 * --p-trim-left-f and --p-trim-left-r trim 5' side on quality. Not the big idea with 16S I understand.
 *
 * Input: .qza with (trimmed) reads.
 * Output:  QIIME2 artifact containing the representative (ASV) sequences.
 *          QIIME2 artifact containing the feature table.
 *          QIIME2 artifact containing the denoising statistics.
 */
process dada2_denoise {
    publishDir "${params.outdir}/denoise", mode: 'copy', overwrite: 'true', saveAs: {filename -> 
        if (filename.indexOf("dada_stats/stats.tsv") == 0)              "abundance_table/unfiltered/dada_stats.tsv"
        else if (filename.indexOf("dada_report.txt") == 0)              "abundance_table/unfiltered/dada_report.txt"
        else if (filename.indexOf("table.qza") == 0)                    "abundance_table/unfiltered/$filename"
        else if (filename.indexOf("rel-table/feature-table.biom") == 0) "abundance_table/unfiltered/rel-feature-table.biom"
        else if (filename.indexOf("abs-table/feature-table.biom") == 0)     "abundance_table/unfiltered/abs-feature-table.biom"
        else if (filename.indexOf("rel-feature-table.tsv") > 0)         "abundance_table/unfiltered/rel-feature-table.tsv"
        else if (filename.indexOf("feature-table.tsv") > 0)             "abundance_table/unfiltered/abs-feature-table.tsv"
        else if (filename.indexOf("rep-seqs.qza") == 0)                 "representative_sequences/unfiltered/rep-seqs.qza"
        else if (filename.indexOf("unfiltered/*"))                      "representative_sequences/$filename"
        else null}
    
    input:
    file(artifact) from trimmed_artifacts

    output:
    file("*")
    file("table.qza") into (table_for_boxplot, table_for_dada_export, diversity_metrics_table, diversity_alpha_table, diversity_beta_table)
    file("rep-seqs.qza") into (repseq_raw_for_classifier, repseq_raw_for_classifier_vsearch, repseq_raw_for_dada_export, repseq_for_tree)
    file("denoising-stats.qza") into denoising_stats_for_metadata_merge
    tuple val("relative"), file("rel-table/feature-table.biom") into relative_biom_table_for_biotaviz
    tuple val("absolute"), file("abs-table/feature-table.biom") into absolute_biom_table_for_biotaviz
//     file("relative-table.qza") into rel_table_for_dada_export

    script:
    """
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${artifact} \
        --p-trim-left-f 0 `# TE : Trims first n-bases from 5' !` \
        --p-trim-left-r 0 `# TE` \
        --p-chimera-method "consensus" \
        --p-trunc-len-f 0 `# TE: Truncates everything after n-bases from 5', NB reads shorter than this number are discarded !` \
        --p-trunc-len-r 0 `# TE` \
        --p-trunc-q 2 `# TE: This value set to 2 was used by older Illumina software as threshold for inserting an N-nucleotide (in this case also truncates the read thereafter), default = 2 !` \
		--p-max-ee-f 3 `# TE: QIIME2 default = 2, DADA2 default = Inf !`\
		--p-max-ee-r 3 `# TE`\
		--p-min-fold-parent-over-abundance 2 `# TE: default = 1`\
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats denoising-stats.qza \
        --verbose \
		> dada_report.txt
    
    #convert absolute (raw) count table to biom format
    qiime tools export \
        --input-path table.qza \
        --output-path abs-table

    #convert absolute (raw) count table to relative abundances
    qiime feature-table relative-frequency \
        --i-table table.qza \
        --o-relative-frequency-table relative-table.qza

    #convert relative count table to biom format
    qiime tools export \
        --input-path relative-table.qza \
        --output-path rel-table
    
    #create representative sequences fasta file in qiime artifact
    qiime feature-table tabulate-seqs \
		--i-data rep-seqs.qza \
		--o-visualization rep-seqs.qzv
    
    #convert artifact to .tsv file
	qiime tools export \
        --input-path rep-seqs.qzv \
		--output-path unfiltered
    """
}

/*
 * Create a phylogenetic tree using MAFFT for the alignment afterwhich masking is performed, Fasttree for the unrooted tree and finally a rooted tree.
 *
 * Input: QIIME2 artifact containing the representative (ASV) sequences.
 * Output: .qza (QIIME2 artifact) for the rooted, .nwk (Newick file) for the created tree artifact.
 */ 
process phylogeny {
    publishDir "${params.outdir}/phylogeny", mode: 'copy', overwrite: 'true'

    input:
    file(repseqs) from repseq_for_tree
    
    output:
    file("rooted-tree.qza") into (diversity_metrics_tree, diversity_alpha_tree, diversity_beta_tree)
    file("phylogenetic_tree/*.nwk") into tree_for_itol
    file("phylogenetic_tree/*")

    when:
    !params.skip_tree

    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${repseqs} \
        --o-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza

    qiime tools export --input-path rooted-tree.qza  \
		--output-path phylogenetic_tree
    """
}

/*
    * Trains a sklearn Naive Bayes classifier using user-provided reference reads and -taxonomy in the QIIME artifact format.
    *
    * Input: QIIME2 artifact containing the reference reads.
    *        QIIME2 artifact containing the reference taxonomy.
    * Output: Trained classifier.qza (QIIME2 artifact) for the rooted, .nwk (Newick file) for the created tree artifact.
    */
if (params.train_classifier){
    process train_taxonomic_classifier {
        publishDir "${params.outdir}/trained_classifier/", mode: 'copy', overwrite: 'true'

        input:
        file(ref_reads) from seqs_train
        file(ref_taxonomy) from taxon_train

        output:
        file ("new_classifier.qza") into new_qiime_classifier

        when:
        params.train_classifier

        script:
        """
        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-reads ${ref_reads} \
            --i-reference-taxonomy ${ref_taxonomy} \
            --o-classifier new_classifier.qza \
            --p-verbose \
            --verbose
        """
    }
    all_classifiers = qiime_classifier.merge(new_qiime_classifier)
} else{
    new_qiime_classifier = Channel.of("empty")
    all_classifiers = qiime_classifier.merge(new_qiime_classifier)
}

if (params.use_taxonomy_classifier){
    /*
    * Classify representative (ASV) sequences using a pre-trained classifier based on the SILVA or GreenGenes 99% OTUs full-length sequences databases.
    *
    * Input: QIIME2 artifact containing the representative (ASV) sequences.
    *        QIIME2 artifact containing the pre-trained classifier.
    * Output: .qza (QIIME2 artifact) of the classified sequences, .tsv (tab-seperated values file) with the classification with it's associated confidence score.
    */ 
    process taxonomic_classifier {
        publishDir "${params.outdir}/classifier/sklearn", mode: 'copy', overwrite: 'true'

        input:
        file(repseq) from repseq_raw_for_classifier
        tuple file(classifier), file(trained_classifier) from all_classifiers
        
        output:
        file("taxonomy_sklearn.qza") into (taxonomy_for_boxplot, taxonomy_for_dada_export, taxonomy_for_metadata_merge)
        file("taxonomy_sklearn.qzv") into taxonomy_for_alpha_diversity
        file("taxonomy/biotaviz_taxonomy.tsv") into (taxonomy_for_biotaviz_rel, taxonomy_for_biotaviz_abs)
	

        when:
        !params.skip_diversity

        script:
        def current_classifier = classifier
        if (trained_classifier != "empty" && trained_classifier ==~ "new_classifier.qza") {
            current_classifier = trained_classifier
        }
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier ${current_classifier} \
            --i-reads ${repseq} \
            --o-classification taxonomy_sklearn.qza

        qiime metadata tabulate \
            --m-input-file taxonomy_sklearn.qza \
            --o-visualization taxonomy_sklearn.qzv
        
        qiime tools export --input-path taxonomy_sklearn.qza  \
            --output-path taxonomy
        qiime tools export --input-path taxonomy_sklearn.qzv  \
            --output-path taxonomy

        cat taxonomy/taxonomy.tsv | sed 's/Feature ID/#OTUID/' | sed 's/Taxon/taxonomy/' | sed 's/Confidence/confidence/' > taxonomy/biotaviz_taxonomy.tsv
        """
    }
} else {
    /*
     * Classify representative (ASV) sequences using a VSEARCH alignment with the SILVA-138 99% OTUs full-length sequences as reference.
     *
     * Input: QIIME2 artifact containing the representative (ASV) sequences.
     *        QIIME2 artifact containing the reference sequences.
     *        QIIME2 artifact containing the reference taxonomy.
     * Output: .qza (QIIME2 artifact) of the classified sequences, .tsv (tab-seperated values file) with the classification with it's associated confidence score.
     */ 
    process taxonomy_vsearch {
        publishDir "${params.outdir}/classifier/vsearch", mode: 'copy', overwrite: 'true'
        storeDir "${params.outdir}/classifier/vsearch"

        input:
        file(repseq) from repseq_raw_for_classifier_vsearch
        file(ref) from seqs
        file(tax) from taxon

        output:
        file("taxonomy_vsearch.qza") into (taxonomy_for_boxplot, taxonomy_for_dada_export, taxonomy_for_metadata_merge)
        file("taxonomy_vsearch.qzv") into taxonomy_for_alpha_diversity
        file("taxonomy/biotaviz_taxonomy.tsv") into (taxonomy_for_biotaviz_rel, taxonomy_for_biotaviz_abs)

        when:
        !params.skip_diversity

        script:
        """
        qiime feature-classifier classify-consensus-vsearch \
            --i-query ${repseq} \
            --i-reference-reads ${ref} \
            --i-reference-taxonomy ${tax} \
            --o-classification taxonomy_vsearch.qza

        qiime metadata tabulate \
            --m-input-file taxonomy_vsearch.qza \
            --o-visualization taxonomy_vsearch.qzv
        
        qiime tools export --input-path taxonomy_vsearch.qza  \
            --output-path taxonomy
        qiime tools export --input-path taxonomy_vsearch.qzv  \
            --output-path taxonomy
        
        cat taxonomy/taxonomy.tsv | sed 's/Feature ID/#OTUID/' | sed 's/Taxon/taxonomy/' | sed 's/Consensus/consensus/' > taxonomy/biotaviz_taxonomy.tsv
        """
    }
}

/*
 * Creates a barplot of the taxonomy distribution.
 *
 * Input: QIIME2 artifact containing classified taxonomy.
 *        QIIME2 artifact containing the feature table.
 *        QIIME2 artifact containing the metadata.
 * Output: .qzv (QIIME2 visualization).
 */ 
process barplot {
    publishDir "${params.outdir}/barplot", mode: 'copy', overwrite: 'true'

    input:
    file (taxonomy) from taxonomy_for_boxplot
    file(table) from table_for_boxplot
    file(metadata) from ch_metadata_for_barplot
    
    output:
    file("taxa-barplots.qzv")
    file("barplot/*")

    when:
    !params.skip_barplot

    script:
    """
     qiime taxa barplot \
        --i-table ${table} \
        --i-taxonomy ${taxonomy} \
        --m-metadata-file ${metadata} \
        --o-visualization taxa-barplots.qzv \
        --verbose

    qiime tools export \
        --input-path taxa-barplots.qzv \
		--output-path barplot
    """
}

/*
 * Creates a BiotaViz file from a taxonomy-annotated (absolute) biom file.
 *
 * Input: Absolute biom table without taxonomy
 *        Taxonomy file with the headers transformed into the right format.
 * Output: .biotaviz to use in downstream visualizations.
 */ 
process biotaviz_absolute {
    publishDir "${params.outdir}/BiotaViz", mode: 'copy', overwrite: 'true'

    input:
    tuple val(label), file(feature_table) from absolute_biom_table_for_biotaviz
    file(clean_taxonomy) from taxonomy_for_biotaviz_abs
     
    output:
    path("*.{biom,biotaviz}*")
    path("absolute-table.biotaviz.txt") into biotaviz_absolute_alldomains

    script:
    """
    biom add-metadata \
        -i ${feature_table} \
        -o ${label}-table-with-taxonomy.biom \
        --observation-metadata-fp ${clean_taxonomy} \
        --sc-separated taxonomy

    python3 $baseDir/bin/biom2biotaviz.py ${label}-table-with-taxonomy.biom > ${label}-table.biotaviz.txt
    """
}

/*
 * Creates a BiotaViz file from a taxonomy-annotated (relative) biom file.
 *
 * Input: Relative biom table without taxonomy
 *        Taxonomy file with the headers transformed into the right format.
 * Output: .biotaviz with absolute and relative counts that may contain multiple domains
 */ 
process biotaviz_relative {
    publishDir "${params.outdir}/BiotaViz", mode: 'copy', overwrite: 'true'

    input:
    tuple val(label), file(feature_table) from relative_biom_table_for_biotaviz
    file(clean_taxonomy) from taxonomy_for_biotaviz_rel

    output:
    path("*.{biom,biotaviz}*")
//     path("relative-table.biotaviz.txt") into biotaviz_for_rankstat1, biotaviz_for_rankstat2

    script:
    """
    biom add-metadata \
        -i ${feature_table} \
        -o ${label}-table-with-taxonomy.biom \
        --observation-metadata-fp ${clean_taxonomy} \
        --sc-separated taxonomy

    python3 $baseDir/bin/biom2biotaviz.py ${label}-table-with-taxonomy.biom > ${label}-table.biotaviz.txt
    """
}

/*
 * Creates a BiotaViz file from a taxonomy-annotated (absolute) biom file with only the Bacteria domain.
 *
 * Input: Absolute biotaviz counts file with taxonomy
 * Output: .biotaviz to use in downstream visualizations.
 */
process biotaviz_relative_bacteria {
    publishDir "${params.outdir}/BiotaViz", mode: 'copy', overwrite: 'true'

    input:
    file(biotaviz_absolute_counts_file) from biotaviz_absolute_alldomains

    output:
    path("absolute-table.biotaviz_relative_abundance.txt") into biotaviz_for_rankstat1, biotaviz_for_rankstat2

    script:
    """
    python3 $baseDir/bin/JOS_Biotaviz_counts_to_abundance.py ${biotaviz_absolute_counts_file}
    """
}

/*
 * Do rankstats on groups defined in mappingfile.
 *
 * Input: BiotaViz file, mapping file with columns of interest starting with "rankstat" 
*/

process rankstat {
    publishDir "${params.outdir}/rankstat", mode: 'copy', overwrite: 'true'

    input:
    file(clean_taxonomy) from biotaviz_for_rankstat1
    file(metadata) from metadata_for_rankstat1
    output:
    file("rankstat.txt")

    script:
    """
    python3 $baseDir/bin/JOS_mapping_rankstat.py -b ${clean_taxonomy} -m ${metadata} -o rankstat.txt
    """
}



/*
 * Create a phylogenetic tree using MAFFT for the alignment afterwhich masking is performed, Fasttree for the unrooted tree and finally a rooted tree.
 *
 * Input: QIIME2 artifact containing the representative (ASV) sequences.
 * Output: .qza (QIIME2 artifact) for the rooted, .nwk (Newick file) for the created tree artifact.
 */ 
process export_dada_output {
    publishDir "${params.outdir}/dada_export", mode: 'copy', overwrite: 'true',
    saveAs: {filename -> 
				 if (filename.indexOf("table/feature-table.biom") == 0)  "abundance_table/filtered/feature-table.biom"
			else if (filename.indexOf("table/feature-table.tsv") == 0)   "abundance_table/filtered/feature-table.tsv"
			else if (filename.indexOf("abs-abund-table-") == 0)          "abundance_table/filtered/$filename"
			else if (filename.indexOf("filtered/*"))                     "representative_sequences/$filename"
			else null} 
    
    input:
    file(table) from table_for_dada_export
    file(repseq) from repseq_raw_for_dada_export
    file(taxonomy) from taxonomy_for_dada_export

    output:
    file("*")
    file("abs-abund-table-${params.taxa_level_for_itol}.tsv") into filter_abundance_feature_table
    file("table/feature-table.tsv") into (data_output_for_minmax, data_output_for_minmax_alpha, data_output_for_minmax_beta)
    
    script:
    """
    #Raw count table in biom format "table/feature-table.biom"
	qiime tools export --input-path ${table}  \
		--output-path table
	#Raw count table "table/feature-table.tsv"
	biom convert -i table/feature-table.biom \
		-o table/feature-table.tsv  \
		--to-tsv
	#Representative sequence fasta file "${params.outdir}/representative_sequences/sequences.fasta"
	qiime feature-table tabulate-seqs  \
		--i-data ${repseq}  \
		--o-visualization rep-seqs.qzv
	qiime tools export --input-path rep-seqs.qzv  \
		--output-path filtered

    ##on several taxa level
	array=( 2 3 4 5 6 7 )
	for i in \${array[@]}
	do
		#collapse taxa
		qiime taxa collapse \
			--i-table ${table} \
			--i-taxonomy ${taxonomy} \
			--p-level \$i \
			--o-collapsed-table table-\$i.qza
		#export to biom
		qiime tools export --input-path table-\$i.qza \
			--output-path table-\$i
		#convert to tab separated text file
		biom convert \
			-i table-\$i/feature-table.biom \
			-o abs-abund-table-\$i.tsv --to-tsv
	done
    """
}

process filter_feature_table{
    publishDir "${params.outdir}/filtered_feature_table_taxonomic_level", mode: 'copy', overwrite: 'true'

    input:
    file(abs_feature_table_level) from filter_abundance_feature_table

    output:
    tuple file("*absolute_feature_table.tsv"), file("*relative_feature_table.tsv") into filtered_abs_and_rel_genus_abundance_for_itol

    script:
    """
    python3 $baseDir/bin/filter_feature_table_taxa.py -i ${abs_feature_table_level} -v
    """
}

/*
 * Calculate non-phylogenetic & phylogenetic alpha- and beta diversity metrics. A python script is used to determine the sample depth.
 *
 * Input: QIIME2 artifact containing the rooted tree.
 *        QIIME2 artifact containing the feature table.
 *        QIIME2 artifact containing the metadata.
 *        Sample depth
 * Output: Multiple .qza (QIIME2 visualization) files with diversity metrics.
 */ 
process all_diversity_metrics {
    publishDir "${params.outdir}/diversity_metrics", mode: 'copy', overwrite: 'true'
    
    input:
    file(metadata) from ch_metadata_for_diversity_core
    file(tree) from diversity_metrics_tree
    file(table) from diversity_metrics_table
    file(stats) from data_output_for_minmax

    output:
    file("diversity_core/*.{qza,qzv}")
    tuple file("diversity_core/shannon_vector.qza"), file("diversity_core/faith_pd_vector.qza") into alpha_metrics_merge_metadata
    file("diversity_core/*")

    when:
    !params.skip_tree || !params.skip_diversity

    script:
    """
    mindepth=\$(python3 $baseDir/bin/table_minmax.py $stats minimum)

    if [ \"\$mindepth\" -gt \"10000\" ]; then echo \"\nUse the sampling depth of \$mindepth for rarefaction\" ; fi
	if [ \"\$mindepth\" -lt \"10000\" -a \"\$mindepth\" -gt \"5000\" ]; then echo \"\n######## WARNING! The sampling depth of \$mindepth is quite small for rarefaction!\" ; fi
	if [ \"\$mindepth\" -lt \"5000\" -a \"\$mindepth\" -gt \"1000\" ]; then echo \"\n######## WARNING! The sampling depth of \$mindepth is very small for rarefaction!\" ; fi
	if [ \"\$mindepth\" -lt \"1000\" ]; then echo \"\n######## ERROR! The sampling depth of \$mindepth seems too small for rarefaction!\" ; fi
    
    qiime diversity core-metrics-phylogenetic \
		--m-metadata-file ${metadata} \
		--i-phylogeny ${tree} \
		--i-table ${table} \
        --p-sampling-depth \$mindepth \
		--output-dir diversity_core \
		--p-n-jobs 1 \
		--quiet
    """
}

/*
 * Merge alpha diversity metrics into one file.
 *
 * Input: Alpha diversity QIIME2 artifact.
 * Output: Combined .qza (QIIME2 artifact), .tsv with all the input metrics.
 */ 
process alpha_rarefaction {
    publishDir "${params.outdir}/diversity_metrics", mode: 'copy', overwrite: 'true'

    input:
    file(table) from diversity_alpha_table
    file(stats) from data_output_for_minmax_alpha
    file(metadata) from ch_metadata_for_alpha_diversity

    output:
    file("alpha_rarefaction.qzv")
    file("alpha_rarefaction/*.csv")
    
    when:
    !params.skip_tree || !params.skip_diversity
    
    script:
    """
    maxdepth=\$(python3 $baseDir/bin/table_minmax.py $stats maximum)

    qiime diversity alpha-rarefaction \
        --i-table ${table} \
        --m-metadata-file ${metadata} \
        --p-max-depth \$maxdepth \
        --o-visualization alpha_rarefaction.qzv
    
    qiime tools export \
        --input-path alpha_rarefaction.qzv \
		--output-path alpha_rarefaction
    """
}

/*
 * Merge alpha diversity metrics into one file.
 *
 * Input: Alpha diversity QIIME2 artifact.
 * Output: Combined .qza (QIIME2 artifact), .tsv with all the input metrics.
 */ 
process beta_rarefaction {
    publishDir "${params.outdir}/diversity_metrics", mode: 'copy', overwrite: 'true'

    input:
    file(table) from diversity_beta_table
    file(metadata) from ch_metadata_for_beta_diversity
    file(stats) from data_output_for_minmax_beta
    file(tree) from diversity_beta_tree
    
    output:
    file("beta_rarefaction.qzv")
    file("beta_rarefaction/*upgma.tre") into clustered_tree_for_itol
    file("beta_rarefaction/*")

    when:
    !params.skip_tree || !params.skip_diversity

    script:
    """
    mindepth=\$(python3 $baseDir/bin/table_minmax.py $stats minimum)

    qiime diversity beta-rarefaction \
        --i-table ${table} \
        --i-phylogeny ${tree} \
        --p-metric 'weighted_unifrac' \
        --p-clustering-method 'upgma' \
        --m-metadata-file ${metadata} \
        --p-sampling-depth \$mindepth \
        --o-visualization beta_rarefaction.qzv
    
    qiime tools export \
        --input-path beta_rarefaction.qzv \
		--output-path beta_rarefaction
    """
}
//
// /*
//  * Creates an iTOL visualization (PDF) from a newick file and additional annotation files.
//  *
//  * Input: Rooted tree (.nwk)
//  *        Additional annotation files
//  * Output: .pdf with an annotated iTOL tree.
//  */
// process itol {
//     publishDir "${params.outdir}/iTOL", mode: 'copy', overwrite: 'true',
//     saveAs: {filename ->
// 			if (filename.indexOf("itol_annotation") == 0)          "annotation/$filename"
// 			else if (filename.indexOf("iTOL-results_") == 0)     "results/$filename"
//             else "misc/$filename"
//             }
//
//
//     input:
//     file(tree) from clustered_tree_for_itol
//     file(metadata) from metadata_for_itol
//     tuple file(absolute_filtered_genus_abundance), file(relative_filtered_genus_abundance) from filtered_abs_and_rel_genus_abundance_for_itol
//
//     output:
//     file("*.{txt,pdf,log,zip,tree}")
//
//     when:
//     !params.skip_tree
//
//     script:
//     """
//     #Tree file must have a .tree extension
//     mv ${tree} ${tree}.tree
//
//     #Create iTOL annotation file
//     python3 $baseDir/bin/taxa_abundance_to_itol_annotation.py -i ${relative_filtered_genus_abundance} --metadata ${metadata} --condition_header ${params.condition_header} -o itol_annotation.txt -v
//
//     #All the data must be zipped. Add annotation files after the tree.
//     zip data.zip ${tree}.tree itol_annotation.txt conditions_itol_annotation.txt
//
//     #Perform curl form request
// 	curl -s --form "zipFile=@data.zip" --form "APIkey=XKkfw0H5bTlS8ftQqWzFtg" --form "projectName=QIIME2" https://itol.embl.de/batch_uploader.cgi > iTOL.log
//     IDtree=`perl $baseDir/bin/get-TreeID.pl iTOL.log`
//     curl -s "https://itol.embl.de/batch_downloader.cgi?tree=\$IDtree&format=pdf&display_mode=1&datasets_visible=0,1,2,3,4&tree_x=220&tree_y=-100" > iTOL-results_\$IDtree.pdf
//     """
// }

/*
 * Merge alpha diversity metrics into one file.
 *
 * Input: Alpha diversity QIIME2 artifact.
 * Output: Combined .qza (QIIME2 artifact), .tsv with all the input metrics.
 */ 
process merge_alpha_diversity_metrics {
    publishDir "${params.outdir}/merged_metrics/alpha_diversity", mode: 'copy', overwrite: 'true'

    input:
    tuple file(shannon), file(faith) from alpha_metrics_merge_metadata
    
    output:
    file("*.qzv")
    file("data/*.tsv")

    script:
    """
    qiime metadata tabulate \
        --m-input-file ${shannon} \
        --m-input-file ${faith} \
        --o-visualization combined-alpha-metadata.qzv

    qiime tools export \
        --input-path combined-alpha-metadata.qzv \
		--output-path data
    """
}

/*
 * Merge QC statistics into one file.
 *
 * Input: QC stats.
 * Output: Combined .qza (QIIME2 artifact), .tsv with all the input metrics.
 */ 
process merge_qc_metrics {
    publishDir "${params.outdir}/merged_metrics/QC", mode: 'copy', overwrite: 'true'

    input:
    file(denoising_stats) from denoising_stats_for_metadata_merge
    
    output:
    file("*/*.{qzv,tsv}")

    script:
    """
    qiime metadata tabulate \
        --m-input-file ${denoising_stats} \
        --o-visualization combined-qc-metadata.qzv
    
    qiime tools export \
        --input-path combined-qc-metadata.qzv \
		--output-path data
    """
}

/*
 * Create .csv files necessary for generating sankey diagrams
 *
 * Input: BiotaViz file, mapping file with columns of interest starting with "rankstat"
*/

process sankey_prep {
    publishDir "${params.outdir}/sankey/preparation-files", mode: 'copy', overwrite: 'true'
    storeDir "${params.outdir}/sankey/preparation-files"

    input:
    file(clean_taxonomy) from biotaviz_for_rankstat2
    file(metadata) from metadata_for_rankstat2

    output:
    file("*.csv") into sankey_prep_files

    script:
    """
    python3 $baseDir/bin/sankey-file-prep.py 0.01 false ${metadata} false ${clean_taxonomy}
    """
}

/*
 * Create sankey diagrams (saved in html files)
 *
 * Input: sankey_prep_datasets (.csv files containing the necessary information); created by process sankey_prep
*/


process sankey {
    publishDir "${params.outdir}/sankey/html-files", mode: 'copy', overwrite: 'true'
    storeDir "${params.outdir}/sankey/html-files"

    input:
    file(prep_file) from sankey_prep_datasets

    output:
    file("${prep_file.baseName}.html")

    script:
    """
    Rscript $baseDir/bin/sankey-diagram-html-generator.R ${prep_file}
    """
}

/*
 * Create sankey diagrams (saved in png files)
 *
 * Input: sankey_html_datasets (.html files containing the sankey diagram); created by process sankey
*/

process sankey_png {
    publishDir "${params.outdir}/sankey/png-files", mode: 'copy', overwrite: 'true'

    input:
    file(html_file) from sankey_html_datasets

    output:
    file("${html_file.baseName}.png")

    script:
    """
    Rscript $baseDir/bin/sankey-diagram-png-generator.R ${html_file}
    """
}

