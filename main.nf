#!/usr/bin/env nextflow

// Create a channel from the list of SRR numbers
SRR_LIST = Channel.from(["SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726"])

// Define a process to download FASTQ files using fasterq-dump
process downloadFastq {
    
    // Input is the SRR number
    input:
    val srr_id

    // Output the FASTQ files
    output:
    file("./fastq_files/${srr_id}_1.fastq")
    file("./fastq_files/${srr_id}_2.fastq")

    // Define the script to be executed
    script:
    """
    # Create the fastq_files directory if it doesn't exist
    mkdir -p ./fastq_files/

    # Check if files are already downloaded
    if [[ ! -f ./fastq_files/${srr_id}_1.fastq || ! -f ./fastq_files/${srr_id}_2.fastq ]]; then
        echo "Starting download for ${srr_id}..."
        fasterq-dump ${srr_id} --split-files --outdir ./fastq_files/ &> ./fastq_files/${srr_id}_log.txt
        echo "Finished downloading ${srr_id}."
    else
        echo "Files for ${srr_id} already exist. Skipping download."
    fi
    """
}

// Workflow execution
workflow {
    downloadFastq(SRR_LIST)
}
