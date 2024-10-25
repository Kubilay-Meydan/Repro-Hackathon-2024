#!/usr/bin/env nextflow

// Define a list of SRR numbers to download
def SRR_LIST = Channel.of("SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726")

// Define a process to download a single FASTQ file using fasterq-dump
process downloadFastq {

    // Input each SRR ID directly from the channel
    input:
    val srr_id

    // Output each downloaded file to the fastq_files directory
    output:
    path "fastq_files/${srr_id}.fastq"

    // Script to download FASTQ file for each SRR ID
    script:
    """
    # Create output directory
    mkdir -p fastq_files

    # Download FASTQ file for the given SRR ID
    echo "Starting download for ${srr_id}..."
    fasterq-dump ${srr_id} --split-spot --skip-technical --outdir fastq_files
    echo "Finished downloading ${srr_id}."
    """
}

// Workflow execution
workflow {
    downloadFastq(SRR_LIST)
}
