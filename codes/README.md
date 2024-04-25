# BusIER

## Table of Content
- <a href="#analyses">Analyses</a>
    - <a href="#download">Downloading Raw Data</a>
    - <a href="#check">Checking Reference Bias on BUSCO</a>

## <a id="analyses">Analyses</a>

### <a id="download">Downloading Raw Data</a>
In this step, we download the reference assembly and short reads from NCBI. The parameters for this step is set in `1_data_download/1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `BusIER/codes`                                                                                                   |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelization                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`            | Metadata file for reference assembly (e.g., `BusIER/data/eucs_refseq.txt`)                                                            |
| `file_shortreads`        | Metadata file for short reads (e.g., `BusIER/data/eucs_shortreads.txt`)                                                               |
| `exe_datasets`           | Executable for NCBI Datasets                                                                                                          |
| `bin_sratoolkit`         | `bin` directory for SRA-Toolkit                                                                                                       |

#### Output
Running the code will create the following folders:
- `refseq/`: folder with all reference assemblies
- `short_reads/`: folder with all raw short reads

### <a id="check">Checking Reference Bias on BUSCO</a>
In this step, we run correlation analysis to check for the presence of reference bias in BUSCO genes. The parameters for this step is set in `2_busco_check/1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `BusIER/codes`                                                                                                   |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelization                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`            | Metadata file for reference assembly (e.g., `BusIER/data/eucs_refseq.txt`)                                                            |
| `file_shortreads`        | Metadata file for short reads (e.g., `BusIER/data/eucs_shortreads.txt`)                                                               |
| `file_adapters`          | Metadata file for sequencing adapters (e.g., `BusIER/data/eucs_adapter.txt`)                                                          |
| `file_refseq_chr`        | Metadata file for reference assembly chromosomes (e.g., `BusIER/data/eucs_refseq_chr.txt`)                                            |
| `exe_adapterremoval`     | Executable for AdapterRemoval                                                                                                         |
| `exe_bwamem2`            | Executable for BWA-MEM2                                                                                                               |
| `exe_samtools`           | Executable for Samtools                                                                                                               |
| `exe_bcftools`           | Executable for Bcftools                                                                                                               |
| `exe_qualimap`           | Executable for QualiMap                                                                                                               |
| `exe_treeshrink`         | Executable for TreeShrink                                                                                                             |
| `exe_gffread`            | Executable for GffRead                                                                                                                |
| `exe_gff2bed`            | Executable for Gff2Bed from BEDOPS                                                                                                    |
| `exe_busco`              | Executable for BUSCO                                                                                                                  |
| `exe_iqtree2`            | Executable for IQ-TREE 2                                                                                                              |
| `exe_mafft`              | Executable for MAFFT                                                                                                                  |
| `busco_lineage`          | Lineage for BUSCO pipeline                                                                                                            |
| `busco_mode`             | Mode for BUSCO pipeline. Options: genome, transcriptome, or protein                                                                   |
| `type`                   | Method to extract BUSCOs from mapped reads. Options: coordinate or pipeline                                                           |

---
*Last update: 25 April 2024 by Jeremias Ivan*