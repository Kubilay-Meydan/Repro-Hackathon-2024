<a name="top"></a>
```
  _____    ______   _____    _____     ____                                        
 |  __ \  |  ____| |  __ \  |  __ \   / __ \                                       
 | |__) | | |__    | |__) | | |__) | | |  | |    ______                            
 |  _  /  |  __|   |  ___/  |  _  /  | |  | |   |______|                           
 | | \ \  | |____  | |      | | \ \  | |__| |                                      
 |_|  \_\ |______| |_| _____|_|_ \_\  \____/     _______   _    _    ____    _   _ 
 | |  | |     /\      / ____| | |/ /     /\     |__   __| | |  | |  / __ \  | \ | |
 | |__| |    /  \    | |      | ' /     /  \       | |    | |__| | | |  | | |  \| |
 |  __  |   / /\ \   | |      |  <     / /\ \      | |    |  __  | | |  | | | . ` |
 | |  | |  / ____ \  | |____  | . \   / ____ \     | |    | |  | | | |__| | | |\  |
 |_|  |_| /_/    \_\  \_____| |_|\_\ /_/    \_\    |_|    |_|  |_|  \____/  |_| \_|
                                                                                   
                                                                                   
```                                   
                                                                                   
                                                                                   
[![Nextflow Version](https://img.shields.io/badge/Nextflow-24.10-brightgreen)](https://www.nextflow.io)
[![Docker Version](https://img.shields.io/badge/Docker-20.10.21-blue)](https://www.docker.com)
[![R Version](https://img.shields.io/badge/R-3.4.1-darkblue)](https://www.r-project.org/)
[![Bowtie Version](https://img.shields.io/badge/Bowtie-0.12.7-green)](http://bowtie-bio.sourceforge.net/index.shtml)
[![DESeq2 Version](https://img.shields.io/badge/DESeq2-1.46.0-red)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
[![SRA Toolkit Version](https://img.shields.io/badge/SRA%20Toolkit-3.1.1-lightgreen)](https://github.com/ncbi/sra-tools)
[![FastQC Version](https://img.shields.io/badge/FastQC-0.12.1-darkblue)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
[![featureCounts Version](https://img.shields.io/badge/featureCounts-1.4.6-red)](http://bioinf.wehi.edu.au/featureCounts/)
[![Samtools Version](https://img.shields.io/badge/Samtools-1.21-darkred)](http://www.htslib.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-turquoise.svg)](https://fr.wikipedia.org/wiki/Licence_MIT)

## 📑 Table of Contents
- [Introduction](#Introduction)
- [Setup and Installation](#setup-and-installation)
- [Contributors](#contributors)
- [Contacts](#contacts)

## 🧬 Introduction <a name="Introduction"></a>

**ReproHackathon RNASeq Pipeline** reproduces findings from [*Peyrusson et al.* (Nature Communications, 2020)](https://www.nature.com/articles/s41467-020-15966-7) on the transcriptomic profile of *Staphylococcus aureus* persisters under antibiotic stress. The study highlights how *S. aureus* persisters—cells that tolerate antibiotics without genetic resistance—may contribute to recurring infections.

Built with [Nextflow](https://www.nextflow.io) and [Docker](https://www.docker.com), this pipeline employs RNA-Seq to study transcriptomic changes, identifying factors influencing antibiotic persistence and tolerance. It covers genome mapping, read counting, and differential expression analysis, capturing the complex adaptation of *S. aureus* under stress.

### Key Features
- **High Reproducibility**: Leveraging Docker containers ensures easy deployment and consistent results across systems.
- **Data Analysis**: Includes genome mapping, read counting, and statistical analysis for identifying differentially expressed genes (DEGs).
- **Modular Workflow**: Nextflow's flexible design allows for easy customization and expansion of the pipeline. Prebuilt Docker images ensure rapid setup.

## 🔧 Setup and Installation <a name="setup-and-installation"></a>

### Requirements
- **Nextflow** `v24.10.0+`
- **Docker** `v20.10.21+`
- **Hardware**: 8 cores & 8 GB RAM minimum

### Installation Instructions
```bash
# Clone this repository
git clone https://github.com/Kubilay-Meydan/ReproHackathon.git
cd ReproHackathon

# Install Nextflow if not already installed
curl -s https://get.nextflow.io | bash

# Ensure Docker is installed and running: https://docs.docker.com/get-docker/

# Launch the pipeline
Nextflow run main.nf
```

## 👥 Contributors <a name="contributors"></a>
This project was developed by:
- [Kubilay Meydan](https://github.com/Kubilay-Meydan)
- [Nathan Carré](https://github.com/Nathan-Carre)
- [Youna Maillé](https://github.com/YounaMKr)
- [Emma Le Roi Pardonche](https://github.com/emmaleroyp)

## 🗨️ Contacts <a name="contacts"></a>

For questions, feel free to reach out to us through GitHub or connect on LinkedIn.

[Back to top](#top)
