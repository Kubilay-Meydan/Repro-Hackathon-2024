#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define a list of SRR numbers to download
workflow {
    def SRR_LIST = Channel.of("SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726")
    
    fastq_files = downloadFastq(SRR_LIST)
    trimmed_files = trimReads(fastq_files)
    genome_fasta = downloadGenome()
    indexed_genome_dir = indexGenome(genome_fasta)
    mapped_files = mapReads(indexed_genome_dir, trimmed_files)
    genome_annot = dwnldAnnotationGTF()
    feature_count_results = featureCount(mapped_files.collect(), genome_annot)


    //ajouter le fichier local comme entrée
    gene_info_file = file("GeneSpecificInformation_NCTC8325.tsv")

     //geneDiffplot added : 
    de_results = diffGenePlot(feature_count_results, gene_info_file)

    // Generate the MA plot dynamically
    ma_plot = generateMAPlot(feature_count_results)
    // ma_plot.view()

}

// Process to download and gzip FASTQ files using SRA Toolkit
process downloadFastq {
    container 'kubilaymeydan/sra-toolkit-docker:latest'
    
    input:
    val srr_id

    output:
    path "fastq_files/${srr_id}.fastq.gz"

    script:
    """
    mkdir -p fastq_files
    if [ ! -f fastq_files/${srr_id}.fastq.gz ]; then
        echo "Starting download for ${srr_id}..."
        fasterq-dump ${srr_id} --split-spot --skip-technical --threads ${task.cpus} --outdir fastq_files
        gzip fastq_files/${srr_id}.fastq
        echo "Finished downloading and gzipping ${srr_id}."
    else
        echo "File for ${srr_id} already exists. Skipping download."
    fi
    """
}

// Process to trim reads using cutadapt
process trimReads {
    container 'kubilaymeydan/cutadapt-docker:latest'

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

// Process to download genome using wget
process downloadGenome {
    container 'kubilaymeydan/bowtie-docker:latest'

    output:
    path "GCF_000013425.1_ASM1342v1_genomic.fna"

    script:
    """
    wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz" -O GCF_000013425.1_ASM1342v1_genomic.fna.gz
    gunzip GCF_000013425.1_ASM1342v1_genomic.fna.gz
    """
}


// New process to index genome using Bowtie
process indexGenome {
   container 'kubilaymeydan/bowtie-docker:latest'

   input:
   path genome_fasta

   output:
   tuple val(genome_fasta.baseName), path("bowtie_index/")

   script:
   """
   mkdir -p bowtie_index
   bowtie-build ${genome_fasta} bowtie_index/${genome_fasta.baseName}
   """
}


process mapReads {
   container 'kubilaymeydan/bowtie-docker:2.0'

   input:
   tuple val(name), path(indexed_genome_dir) // Chemin du dossier contenant les fichiers d'index (bowtie_index)
   path trimmed_files

   output:
   path "mapped_reads/${trimmed_files.baseName}.bam"

   script:
   """
   mkdir -p mapped_reads
   zcat ${trimmed_files} | bowtie -q ${indexed_genome_dir}/${name} -p ${task.cpus} --sam - | samtools view -bS - > mapped_reads/${trimmed_files.baseName}.bam   
   """
}

process dwnldAnnotationGTF {
    container 'kubilaymeydan/subreads-docker:latest'

    output:
    path "GCF_000013425.1_ASM1342v1_genomic.gtf"

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gtf.gz
    gunzip GCF_000013425.1_ASM1342v1_genomic.gtf.gz
    """
} 

process featureCount {
    container 'kubilaymeydan/subreads-docker:latest'

    input:
    path mapped_files
    path genome_annot

    output:
    path "featurecount_files/merged_feature_counts.txt"

    publishDir "${baseDir}/results", mode: 'copy'

    script:
    """
    mkdir -p featurecount_files
    featureCounts -t gene -g gene_id -s 1 -a ${genome_annot} -o featurecount_files/merged_feature_counts.txt ${mapped_files.join(' ')}
    """
}

process generateMAPlot {
    container 'kubilaymeydan/r-docker:latest'

    input:
    path feature_counts

    output:
    path "*.pdf"
    path "*.png"

    publishDir "${baseDir}/results", mode: 'copy' 

    script:
    """
    Rscript --vanilla /workspace/ma_plot.r ${feature_counts}
    """
}


process diffGenePlot {
    container 'kubilaymeydan/r-docker:latest'

    input:
    path feature_counts
    path gene_info_file

    output:
    path "plot_translation.png"

    publishDir "${baseDir}/results", mode: 'copy'
    

    script:
    """
    Rscript --vanilla /workspace/diff_gene_plot.r ${feature_counts} ${gene_info_file} 
    """
}





