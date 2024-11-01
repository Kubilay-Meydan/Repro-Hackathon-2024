#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define a list of SRR numbers to download
workflow {
    // List of SRR IDs as a channel
    def SRR_LIST = Channel.of("SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726")

    // Run downloadFastq and store output in fastq_files channel
    fastq_files = downloadFastq(SRR_LIST)

    // Download reference genome
    reference_fasta = downloadReference()

    // Pass fastq_files to trimReads
    trimmed_files = trimReads(fastq_files)

    // Display or save trimmed files as needed
    trimmed_files.view()
}

// Process to download and gzip FASTQ files using SRA Toolkit
process downloadFastq {
    container 'sra-toolkit-docker'
    
    input:
    val srr_id

    output:
    path "fastq_files/${srr_id}.fastq.gz"

    script:
    """
    mkdir -p fastq_files
    if [ ! -f fastq_files/${srr_id}.fastq.gz ]; then
        echo "Starting download for ${srr_id}..."
        fasterq-dump ${srr_id} --split-spot --skip-technical --outdir fastq_files
        gzip fastq_files/${srr_id}.fastq
        echo "Finished downloading and gzipping ${srr_id}."
    else
        echo "File for ${srr_id} already exists. Skipping download."
    fi
    """
}

// Process to download reference genome
process downloadReference{
    // si wget dispo sur sra-toolkit
    container "sra-toolkit_docker"

    output:
    path "reference/reference.fasta"

    script:
    """
    mkdir -p reference
    wget -q -O reference/reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
    """
}

// Process to trim reads using cutadapt
process trimReads {
    container 'cutadapt-docker'

    input:
    path fastq_file

    output:
    path "trimmed_reads/${fastq_file.baseName}_trimmed.fastq.gz"

    script:
    """
    cutadapt --version
    mkdir -p trimmed_reads
    cutadapt -q 20 -m 25 -o trimmed_reads/${fastq_file.baseName}_trimmed.fastq.gz ${fastq_file}
    """
}
