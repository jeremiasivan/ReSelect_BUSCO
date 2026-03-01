# ReSelect BUSCO

## Table of Content
- <a href="#prepare">Data Preparation</a>
    - <a href="#owndata">Using Your Own Datasets</a>
- <a href="#check">Checking Reference Bias on BUSCO</a>

## <a id="prepare">Data Preparation</a>
In this step, we download the reference genomes and short reads from NCBI. Then, we performed quality-control (QC) on the short reads and map them to all of the available references. The parameters for this step is set in `1_data_download/1_main.Rmd`.

| Parameters                 | Definition                                                                                                                            |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                  | Directory for folder `ReSelect_BUSCO/codes/`                                                                                          |
| `prefix`                   | Prefix for output files and folder                                                                                                    | 
| `outdir`                   | Output directory                                                                                                                      |
| `thread`                   | Number of threads for parallelisation                                                                                                 |
| `redo`                     | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`              | Metadata file for reference assembly (e.g., `ReSelect_BUSCO/data/eucs_refseq.txt`)                                                    |
| `file_shortreads`          | Metadata file for short reads (e.g., `ReSelect_BUSCO/data/eucs_shortreads.txt`)                                                       |
| `file_adapters`            | Metadata file for sequencing adapters (e.g., `ReSelect_BUSCO/data/eucs_adapter.txt`)                                                  |
| `skip_refseq_download`     | Skip download step for reference genomes                                                                                              |
| `skip_shortreads_download` | Skip download step for short reads data                                                                                               |
| `skip_shortreads_qc`       | Skip QC step for short reads data                                                                                                     |
| `exe_datasets`             | Executable for NCBI Datasets                                                                                                          |
| `bin_sratoolkit`           | `bin/` directory for SRA-Toolkit                                                                                                      |
| `exe_adapterremoval`       | Executable for AdapterRemoval                                                                                                         |
| `exe_bwamem2`              | Executable for BWA-MEM2                                                                                                               |
| `exe_samtools`             | Executable for Samtools                                                                                                               |
| `exe_bcftools`             | Executable for Bcftools                                                                                                               |
| `exe_qualimap`             | Executable for QualiMap                                                                                                               |
| `min_read_quality`         | Minimum quality for short reads                                                                                                       |
| `thread_adapterremoval`    | Number of threads for AdapterRemoval                                                                                                  |
| `thread_bwamem2`           | Number of threads for BWA-MEM2                                                                                                        |
| `thread_samtools`          | Number of threads for Samtools and Bcftools                                                                                           |
| `thread_qualimap`          | Number of threads for QualiMap                                                                                                        |

### Output
Running the code will create the following folders in `outdir/prefix`:
- `refseq/`: folder with all reference assemblies
- `short_reads/`: folder with all raw short reads
    - `filtered/`: folder with all filtered short reads from AdapterRemoval
- `readmap/`: folder with the output of mapping (i.e., BAM, VCF, and consensus FASTA sequence)
    - `qualimap/`: folder with the output of QualiMap for each mapped reads
    - `metadata.tsv`: file with the FASTA directories of all references
    - `summary.tsv`: file with the summary coverage for all mapped reads
- `logs/`: folder that contains the log files for each run

### <a id="owndata">Using Your Own Datasets</a>
If you have your own dataset, please follow these steps (*I am sorry in advance for the complicated folder structures as I follow default outputs from NCBI Datasets, SRA Toolkit, and AdapterRemoval*):
1. For example, I have a paired-end short-read data (`newsp_R1.fastq` and `newsp_R2.fastq`) and three reference genomes (`ref1.fna`, `ref2.fna`, `ref3.fna`).
2. For `file_shortreads` parameter, add each short-read dataset (one per line):

    | id          | species          | is_check          |
    | ----------- | ---------------- | ----------------- |
    | `newsp`     | `newspecies`     | `TRUE`            |


3. For `file_refseq` parameter, add each reference genome (one per line):

    | id          | species          | grouping          |
    | ----------- | ---------------- | ----------------- |
    | `ref1`      | `reference1`     | `group1`          |
    | `ref2`      | `reference2`     | `group1`          |
    | `ref3`      | `reference3`     | `group2`          |

4. For short-read dataset:
    - If you have *raw* sequencing reads, store them in: `$outdir/short_reads/$species/`. <br><br>
    Following the example, if the `$outdir` parameter is `/home/jeremias/output`, then the directories of the short-read files are:
        - `/home/jeremias/output/short_reads/newspecies/newsp_R1.fastq`
        - `/home/jeremias/output/short_reads/newspecies/newsp_R2.fastq`
        - **Ensure all files end with `.fastq`**

    - If you have *cleaned* sequencing reads, store them in: `$output/short_reads/filtered/`. <br><br>
    Following the example, if the `$outdir` parameter is `/home/jeremias/output`, then the directories of the short-read files are:
        - `/home/jeremias/output/short_reads/filtered/newspecies.pair1.truncated`
        - `/home/jeremias/output/short_reads/filtered/newspecies.pair2.truncated`
        - **Ensure all files end with `.pairX.truncated`, where `X` represents either 1 or 2**

5. For reference genomes, store each in: `$outdir/refseq/$species/ncbi_dataset/data/$id/`. <br><br>
   Following the example, if the `$outdir` parameter is `/home/jeremias/output`, then the directories of the reference genomes are:
    - `/home/jeremias/output/refseq/reference1/ncbi_dataset/data/ref1/ref1.fna`
    - `/home/jeremias/output/refseq/reference2/ncbi_dataset/data/ref2/ref2.fna`
    - `/home/jeremias/output/refseq/reference3/ncbi_dataset/data/ref3/ref3.fna`
    - **Ensure all files end with `.fna`**

6. Set `skip_refseq_download`, `skip_shortreads_download`, and/or `skip_shortreads_qc` parameters to be `TRUE`

## <a id="check">Checking Reference Bias on BUSCO</a>
In this step, we run correlation analysis to check for the extent of reference bias in BUSCO and assess if it changes the BUSCO tree topology. The parameters for this step is set in `2_busco_check/1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `ReSelect_BUSCO/codes/`                                                                                          |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelisation                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`            | Metadata file for reference genomes (e.g., `ReSelect_BUSCO/data/eucs_refseq.txt`)                                                     |
| `file_shortreads`        | Metadata file for short reads (e.g., `ReSelect_BUSCO/data/eucs_shortreads.txt`)                                                       |
| `file_genome_treefile`   | Tree file for the reference genomes (e.g., `ReSelect_BUSCO/data/eucs.treefile`)                                                       |
| `exe_busco`              | Executable for BUSCO                                                                                                                  |
| `exe_gff2bed`            | Executable for Gff2Bed from BEDOPS                                                                                                    |
| `exe_samtools`           | Executable for Samtools                                                                                                               |
| `exe_qualimap`           | Executable for QualiMap                                                                                                               |
| `exe_mafft`              | Executable for MAFFT                                                                                                                  |
| `exe_trimal`             | Executable for TrimAl                                                                                                                 |
| `exe_iqtree2`            | Executable for IQ-TREE 2                                                                                                              |
| `exe_treeshrink`         | Executable for TreeShrink                                                                                                             |
| `exe_astral`             | Executable for ASTRAL                                                                                                                 |
| `busco_lineage`          | Lineage for BUSCO pipeline                                                                                                            |
| `busco_mode`             | Mode for BUSCO pipeline. Options: genome, transcriptome, or protein                                                                   |
| `thread_busco`           | Number of threads for BUSCO                                                                                                           |
| `min_busco_depth`        | Minimum BUSCO depth for mapped reads                                                                                                  |
| `is_astral_constrained`  | Constrain the ASTRAL tree topology based on `file_genome_treefile`                                                                    |
| `outgroup`               | Outgroup(s) from the list of references (optional)                                                                                    |
| `focal_species`          | Focal species for the correlation analyses (optional)                                                                                 |

#### Output
Running the code will create the following folders in `outdir/prefix`:
- `busco_lineage/`: folder with the lineage dataset for running BUSCO pipeline
- `busco_check/`
    - `busco_refseq/`: folder with all BUSCO runs on individual reference genome
        - `fasta/`: folder with all BUSCO sequences inferred from BUSCO GFF files. Applicable only for `type==coordinate`.
        - `metadata.tsv`: file with the error status for each BUSCO for each reference
    - `short_reads/`: folder with all BUSCO sequences for all mapped reads
        - `metadata.tsv`: file with the error status for each BUSCO for each mapped reads
    - `trees/`: folder with all BUSCO alignments and trees for reference genomes and mapped reads
    - `summary/`: folder with summary files from correlation analyses
        - `correlation_figs/`: folder with scatter plots from the correlation analysis
        - `prefix.cor.sumtable`: file with the results of the correlation analyses
        - `prefix.dist.tiff`: file with the phylogenetic distances of mapped reads that come from the same species
    - `busco_tree/`: folder with output from the BUSCO tree analyses
        - `all/`, `bias/`, and `nonbias/`: folder with individual analysis for different sets of loci
            - `alignment/`: folder with locus sequences
            - `per_reference/`: folder with one reconstruction per run
            - `trees`: folder with locus trees and ASTRAL tree

---
*Last update: 01 March 2026 by Jeremias Ivan*