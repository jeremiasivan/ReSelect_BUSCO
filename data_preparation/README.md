# PhyloRBT Data Preparation

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#readmap">Mapping Short Reads</a>
- <a href="#busco">Extracting BUSCO Loci</a>

## <a id="prereqs">Prerequisites</a>
In order to prepare BUSCO locus alignments and trees used in PhyloRBT, additional software are required. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also provide the executable paths on `run_pipeline.R`

### Software
|    Name    |                                    Website                                     |                             Anaconda                             |
| ---------- |:------------------------------------------------------------------------------:|:----------------------------------------------------------------:|
| BCFtools   | <a href="https://github.com/samtools/bcftools">Link</a>                        | <a href="https://anaconda.org/bioconda/bcftools">Link</a>        |
| BUSCO      | <a href="https://busco.ezlab.org">Link</a>                                     | <a href="https://anaconda.org/bioconda/busco">Link</a>           |
| BWA-MEM2   | <a href="https://github.com/bwa-mem2/bwa-mem2">Link</a>                        | <a href="https://anaconda.org/bioconda/bwa-mem2">Link</a>        |
| Gff2BED    | <a href="https://bedops.readthedocs.io/en/latest/index.html">Link</a>          | <a href="https://anaconda.org/bioconda/gff2bed">Link</a>         |
| MAFFT      | <a href="https://github.com/GSLBiotech/mafft">Link</a>                         | <a href="https://anaconda.org/bioconda/mafft">Link</a>           |
| QualiMap   | <a href="http://qualimap.conesalab.org/">Link</a>                              | <a href="https://anaconda.org/bioconda/qualimap">Link</a>        |
| SAMtools   | <a href="https://github.com/samtools/samtools">Link</a>                        | <a href="https://anaconda.org/bioconda/samtools">Link</a>        |
| TreeShrink | <a href="https://github.com/uym2/TreeShrink">Link</a>                          | <a href="https://anaconda.org/bioconda/treeshrink">Link</a>      |
| TrimAl     | <a href="https://github.com/inab/trimal">Link</a>                              | <a href="https://anaconda.org/bioconda/trimal">Link</a>          |

If using `conda`, you can use the following command to install the software:
```
conda install -c bioconda bcftools busco bwa-mem2 gff2bed mafft qualimap samtools treeshrink trimal
```
    
If `bcftools` returns a `libgsl.so.25` error, you can either download the software <a href="https://www.htslib.org/download/">here</a>, or try to set the `conda` channel priorities before installing any package:
```
conda config --prepend channels r
conda config --prepend channels bioconda
conda config --prepend channels conda-forge
conda config --set channel_priority strict
```

## <a id="readmap">Mapping Short Reads</a>
In this step, we download the reference genomes and short reads from NCBI. Then, we performed quality-control (QC) on the short reads and map them to all of the available references. The parameters for this step is set in `1_data_download/1_main.Rmd`.

| Parameters                 | Definition                                                                                                                            |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                  | Directory for folder `PhyloRBT/`                                                                                                      |
| `prefix`                   | Prefix for output files and folder                                                                                                    | 
| `outdir`                   | Output directory                                                                                                                      |
| `thread`                   | Number of threads for parallelisation                                                                                                 |
| `redo`                     | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`              | Metadata file for reference assembly (e.g., `ReSelect_BUSCO/data/eucs_refseq.txt`)                                                    |
| `file_shortreads`          | Metadata file for short reads (e.g., `ReSelect_BUSCO/data/eucs_shortreads.txt`)                                                       |
| `exe_bwamem2`              | Executable for BWA-MEM2                                                                                                               |
| `exe_samtools`             | Executable for Samtools                                                                                                               |
| `exe_bcftools`             | Executable for Bcftools                                                                                                               |
| `exe_qualimap`             | Executable for QualiMap                                                                                                               |
| `thread_bwamem2`           | Number of threads for BWA-MEM2                                                                                                        |
| `thread_samtools`          | Number of threads for Samtools and Bcftools                                                                                           |
| `thread_qualimap`          | Number of threads for QualiMap                                                                                                        |

### Output
Running the code will create the following folders in `outdir/prefix`:
- `readmap/`: folder with the output of mapping (i.e., BAM, VCF, and consensus FASTA sequence)
    - `qualimap/`: folder with the output of QualiMap for each mapped reads
    - `metadata.tsv`: file with the FASTA directories of all references
    - `summary.tsv`: file with the summary coverage for all mapped reads

## <a id="busco">Extracting BUSCO Loci</a>
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
| `exe_busco`              | Executable for BUSCO                                                                                                                  |
| `exe_gff2bed`            | Executable for Gff2Bed from BEDOPS                                                                                                    |
| `exe_samtools`           | Executable for Samtools                                                                                                               |
| `exe_qualimap`           | Executable for QualiMap                                                                                                               |
| `exe_mafft`              | Executable for MAFFT                                                                                                                  |
| `exe_trimal`             | Executable for TrimAl                                                                                                                 |
| `exe_iqtree2`            | Executable for IQ-TREE 2                                                                                                              |
| `exe_treeshrink`         | Executable for TreeShrink                                                                                                             |
| `busco_lineage`          | Lineage for BUSCO pipeline                                                                                                            |
| `busco_mode`             | Mode for BUSCO pipeline. Options: genome, transcriptome, or protein                                                                   |
| `thread_busco`           | Number of threads for BUSCO                                                                                                           |
| `min_busco_depth`        | Minimum BUSCO depth for mapped reads                                                                                                  |

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

---
*Last update: 23 July 2026 by Jeremias Ivan*