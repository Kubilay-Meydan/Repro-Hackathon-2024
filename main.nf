#!/usr/bin/env nextflow

// Define a process to download FASTQ files using fasterq-dump
process downloadFastq {

    // Output directory where all files will be saved
    output:
    path "fastq_files"

    // Script to download FASTQ files for each SRR ID
    script:
    """
    # Define the list of SRR numbers directly within the process
    SRR_LIST=("SRR10379721" "SRR10379722" "SRR10379723" "SRR10379724" "SRR10379725" "SRR10379726")

    # Create output directory
    mkdir -p fastq_files

    # Loop through each SRR and download concurrently, redirecting logs
    for srr_id in \${SRR_LIST[@]}; do
        echo "Starting download for \${srr_id}..."
        fasterq-dump \${srr_id} --split-spot --skip-technical --outdir fastq_files > fastq_files/\${srr_id}.log 2>&1 &
    done

    # Wait for all background tasks to complete
    wait
    echo "All downloads completed."

    # Remove all .log files after downloads are complete
    rm -f fastq_files/*.log
    """
}

// Workflow execution
workflow {
    downloadFastq()
}
