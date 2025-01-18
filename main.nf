#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define input parameters
params.inputFile = '/data/volume/bacterial_dna.fasta' // Path to the input FASTA file in the mounted volume
params.cutoff = 0.5 // Default GC content cutoff
params.outputDir = '/data/volume' // Path to the output directory in the mounted volume

// Process to calculate GC content
process calculate_gc_content {
    // Specify output directory
    publishDir params.outputDir, mode: 'copy' // Copy output files to the output directory

    // Define inputs
    input:
    path fasta_file 
    val gc_cutoff 

    // Define outputs
    output:
    path 'output.txt'

    // Script to process the FASTA file
    script:
    """
    awk '/^>/ {if (seq) {print header; print seq} header=\$0; seq=""; next} {seq=seq\$0} END {if (seq) {print header; print seq}}' $fasta_file | \
    awk -v cutoff=$gc_cutoff '{
        if (NR % 2 == 1) {header=\$0; next}
        seq=\$0;
        g_count=gsub(/G/, "G", seq);
        c_count=gsub(/C/, "C", seq);
        total_count=length(seq);
        gc_content=(g_count + c_count) / total_count;
        if (gc_content > cutoff) {
            print header;
            print seq;
        }
    }' > output.txt
    """
}

// Main workflow definition
workflow {
    inputChannel = Channel.fromPath(params.inputFile)
    calculate_gc_content(inputChannel, params.cutoff)
}