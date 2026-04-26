# TopHat2/Cufflinks scRNA-seq Pipeline  
## Replication of Camp et al. (2015) Figure 3D  
### University of Leicester — BS7120 Steered Research Project Coursework  

---

## Overview

This repository contains the **original RNA-seq pipeline** used to erreplicate Figure 3D from:

> Camp, J.G. et al. (2015). *Human cerebral organoids recapitulate gene expression programs of fetal neocortex development*.  
> PNAS, 112(51), 15672–15677.  
> https://doi.org/10.1073/pnas.1520760112  

The original paper used **TopHat2 for alignment** and **Cufflinks for transcript quantification**. This pipeline closely follows that approach, with downstream clustering and visualisation performed using **Seurat**.

---

## Data

Raw sequencing data: GEO accession [GSE75140](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75140)

- 734 single cells (508 organoid, 226 fetal neocortex)  
- Paired-end reads (96–100 bp)  
- Illumina HiSeq 2500  
- Fluidigm C1 platform  

Reference genome: **GENCODE GRCh38 release 22**

- Genome FASTA: `GRCh38.primary_assembly.genome.fa`  
- GTF annotation: `gencode.v22.primary_assembly.annotation.gtf`  
- https://www.gencodegenes.org/human/release_22.html  

---

## Scripts 
Run in order of top to bottom:
| Script | Description |
|--------|------------|
| `fastq_convert.sh` | Convert the downloaded SRA files to FASTQ format |
| `retry_missing_fastq.sh` | Retry failed or incomplete FASTQ conversions |
| `nochr.py` | Remove chromosome prefixes from GTF for compatibility |
| `tophat_all2.slurm` | Align reads to GRCh38 using TopHat2 (HPC job script) |
| `cufflinks_top.sh` | Quantify gene expression (FPKM values) |
| `expression_matrix.py` | Merge Cufflinks outputs into a unified expression matrix |
| `v5_steered.R` | Seurat analysis (PCA, clustering, t-SNE, marker genes) | 

## Software Versions

| Tools | Version |
|------|--------|
| SRA Toolkit | latest |
| TopHat2 | 2.1.x |
| Cufflinks | 2.2.x |
| Python | 3. and 2.7x |
| R | 4.5.x |
| Seurat | v5 |

