# ReSelect BUSCO

## Table of Content
- <a href="#analyses">Analyses</a>
    - <a href="#prepare">Data Preparation</a>
    - <a href="#check">Checking Reference Bias on BUSCO</a>

## <a id="analyses">Analyses</a>

### <a id="prepare">Data Preparation</a>
In this step, we download the reference assembly and short reads from NCBI. Then, we performed quality-control (QC) on the short reads and map them to all of the available references. The parameters for this step is set in `1_data_download/1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `BusIER/codes`                                                                                                   |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelisation                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`            | Metadata file for reference assembly (e.g., `BusIER/data/eucs_refseq.txt`)                                                            |
| `file_shortreads`        | Metadata file for short reads (e.g., `BusIER/data/eucs_shortreads.txt`)                                                               |
| `file_adapters`          | Metadata file for sequencing adapters (e.g., `BusIER/data/eucs_adapter.txt`)                                                          |
| `exe_datasets`           | Executable for NCBI Datasets                                                                                                          |
| `bin_sratoolkit`         | `bin` directory for SRA-Toolkit                                                                                                       |
| `exe_adapterremoval`     | Executable for AdapterRemoval                                                                                                         |
| `exe_bwamem2`            | Executable for BWA-MEM2                                                                                                               |
| `exe_samtools`           | Executable for Samtools                                                                                                               |
| `exe_bcftools`           | Executable for Bcftools                                                                                                               |
| `exe_qualimap`           | Executable for QualiMap                                                                                                               |
| `min_read_quality`       | Minimum quality for short reads                                                                                                       |
| `thread_adapterremoval`  | Number of threads for AdapterRemoval                                                                                                  |
| `thread_bwamem2`         | Number of threads for BWA-MEM2                                                                                                        |
| `thread_samtools`        | Number of threads for Samtools and Bcftools                                                                                           |
| `thread_qualimap`        | Number of threads for QualiMap                                                                                                        |

#### Output
Running the code will create the following folders in `outdir/prefix`:
- `refseq/`: folder with all reference assemblies
- `short_reads/`: folder with all raw short reads
    - `filtered/`: folder with all filtered short reads from AdapterRemoval
- `readmap/`: folder with the output of mapping (i.e., BAM, VCF, and consensus FASTA sequence)
    - `qualimap/`: folder with the output of QualiMap for each pair of short reads and reference

### <a id="check">Checking Reference Bias on BUSCO</a>
In this step, we run correlation analysis to check for the extent of reference bias in BUSCO and assess if it changes the BUSCO tree topology. The parameters for this step is set in `2_busco_check/1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `BusIER/codes`                                                                                                   |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelisation                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`            | Metadata file for reference assembly (e.g., `BusIER/data/eucs_refseq.txt`)                                                            |
| `file_shortreads`        | Metadata file for short reads (e.g., `BusIER/data/eucs_shortreads.txt`)                                                               |
| `file_genome_treefile`   | Tree file for the reference sequences (e.g., `BusIER/data/eucs.treefile`)                                                             |
| `file_buscotree`         | Configuration file for the BUSCO tree analyses (e.g., `BusIER/data/eucs_buscotree.txt`)                                               |
| `exe_busco`              | Executable for BUSCO                                                                                                                  |
| `exe_gffread`            | Executable for GffRead                                                                                                                |
| `exe_gff2bed`            | Executable for Gff2Bed from BEDOPS                                                                                                    |
| `exe_samtools`           | Executable for Samtools                                                                                                               |
| `exe_mafft`              | Executable for MAFFT                                                                                                                  |
| `exe_iqtree2`            | Executable for IQ-TREE 2                                                                                                              |
| `exe_treeshrink`         | Executable for TreeShrink                                                                                                             |
| `exe_astral`             | Executable for ASTRAL                                                                                                                 |
| `busco_lineage`          | Lineage for BUSCO pipeline                                                                                                            |
| `busco_mode`             | Mode for BUSCO pipeline. Options: genome, transcriptome, or protein                                                                   |
| `type`                   | Method to extract BUSCOs from mapped reads. Options: coordinate or pipeline                                                           |
| `thread_busco`           | Number of threads for BUSCO                                                                                                           |
| `min_busco_depth`        | Minimum BUSCO depth for mapped reads                                                                                                  |
| `max_busco_depth`        | Maximum BUSCO depth for mapped reads                                                                                                  |
| `busco_tree_mode`        | Mode for BUSCO tree analyses. Options: control, random, or oneref                                                                     |
| `outgroup`               | Outgroup(s) from the list of references (optional)                                                                                    |

#### Output
Running the code will create the following folders in `outdir/prefix`:
- `busco_lineage/`: folder with the lineage dataset for running BUSCO pipeline
- `busco_check/`
    - `busco_refseq/`: folder with all BUSCO runs on individual reference assembly
        - `fasta/`: folder with all BUSCO sequences inferred from BUSCO GFF files. Applicable only for `type==coordinate`.
    - `busco_coordinate/` and/or `busco_pipeline/`
        - `short_reads/`: folder with all BUSCO sequences for all mapped reads
        - `trees/`: folder with all BUSCO alignments and trees for reference assemblies and mapped reads
        - `summary/`: folder with summary files from correlation analyses
            - `correlation_figs/`: folder with scatter plots from the correlation analysis
        - `busco_tree/`: folder with output from the BUSCO tree analyses
            - `control/`, `random/`, and/or `oneref/`: folder with individual analysis across number of replicates

---
*Last update: 09 December 2024 by Jeremias Ivan*