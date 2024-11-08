#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define a list of SRR numbers to download
workflow {
    // List of SRR IDs as a channel
    def SRR_LIST = Channel.of("SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726")

    // Run downloadFastq and store output in fastq_files channel
    fastq_files = downloadFastq(SRR_LIST)

    // Pass fastq_files to trimReads
    trimmed_files = trimReads(fastq_files)

    // Display or save trimmed files as needed
    trimmed_files.view()

    // Run downloadGenome process
    genome_fasta = downloadGenome()

    // Run indexGenome process with genome_fasta as input
    indexed_genome_dir = indexGenome(genome_fasta)  // Récupère le dossier de sortie d'index

    // Mapping des lectures avec les fichiers indexés et trimmed
    mapped_files = mapReads(indexed_genome_dir, trimmed_files)    

    // Afficher les fichiers d'alignement générés
    mapped_files.view()

    // Download genome annote
    genome_annot = dwnldAnnotationGTF()

    // Mapped file to featureCount
    feature_count_results = featureCount(mapped_files,genome_annot)
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

// New process to download genome using wget
process downloadGenome {
    container 'bowtie-docker'

    output:
    path "CP000253.1.fasta"

    script:
    """
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta&retmode=text" -O CP000253.1.fasta
    """
}

// New process to index genome using Bowtie
process indexGenome {
    container 'bowtie-docker'

    input:
    path genome_fasta

    output:
    path "bowtie_index/"

    script:
    """
    mkdir -p bowtie_index
    bowtie-build ${genome_fasta} bowtie_index/CP000253.1
    """
}

process mapReads {
    container 'bowtie-docker'

    input:
    path indexed_genome_dir // Chemin du dossier contenant les fichiers d'index (bowtie_index)
    path trimmed_files

    output:
    path "mapped_reads/${trimmed_files.baseName}.sam"

    script:
    """
    mkdir -p mapped_reads
    zcat ${trimmed_files} | bowtie -q ${indexed_genome_dir}/CP000253.1 --sam - > mapped_reads/${trimmed_files.baseName}.sam
    """
}

process dwnldAnnotationGTF {
    container 'subreads-docker'

    output:
    path "GCF_000013425.1_ASM1342v1_genomic.gtf"

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gtf.gz
    gunzip GCF_000013425.1_ASM1342v1_genomic.gtf.gz
    """
}

process featureCount {
    container 'subreads-docker'

    input:
    path mapped_file
    path genome_annot

    output:
    path "featurecount_files/${mapped_file.baseName}.txt"

    script:
    """
    mkdir -p featurecount_files
    featureCounts -t gene -g gene_id -s 1 -a ${genome_annot} -o featurecount_files/${mapped_file.baseName}.txt ${mapped_file}
    """
}


