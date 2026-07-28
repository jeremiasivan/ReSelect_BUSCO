# PhyloRBT Data Preparation

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>
    - <a href="#readmap">Mapping Short Reads</a>
    - <a href="#busco">Extracting BUSCO Loci</a>
    - <a href="#runpipeline">Running Both Analyses</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
For PhyloRBT data preparation, additional software are required. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also provide the executable paths on `run_pipeline.R`

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
conda install -c conda-forge -c bioconda bcftools bedops busco bwa-mem2 mafft qualimap samtools treeshrink trimal
```
    
If `bcftools` returns a `libgsl.so.25` error, you can either download the software <a href="https://www.htslib.org/download/">here</a>, or try to set the `conda` channel priorities before installing any package:
```
conda config --prepend channels r
conda config --prepend channels bioconda
conda config --prepend channels conda-forge
conda config --set channel_priority strict
```

## <a id="#genpipe">General Pipeline</a>

### <a id="readmap">Mapping Short Reads</a>
In this step, we map each set of short-reads to all available reference genomes. The parameters for this step is set in `1_readmap/1_main.Rmd`.

| Parameters                 | Definition                                                                                                                            |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                  | Directory for folder `PhyloRBT/`                                                                                                      |
| `prefix`                   | Prefix for output files and folder                                                                                                    | 
| `outdir`                   | Output directory                                                                                                                      |
| `thread`                   | Number of threads for parallelisation                                                                                                 |
| `redo`                     | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`              | Metadata file for reference assembly (e.g., `PhyloRBT/files/eucs_reference.txt`)                                                      |
| `file_shortreads`          | Metadata file for short reads (e.g., `PhyloRBT/files/eucs_shortreads.txt`)                                                            |
| `exe_bwamem2`              | Executable for BWA-MEM2                                                                                                               |
| `exe_samtools`             | Executable for Samtools                                                                                                               |
| `exe_bcftools`             | Executable for Bcftools                                                                                                               |
| `exe_qualimap`             | Executable for QualiMap                                                                                                               |
| `thread_bwamem2`           | Number of threads for BWA-MEM2                                                                                                        |
| `thread_samtools`          | Number of threads for Samtools and Bcftools                                                                                           |
| `thread_qualimap`          | Number of threads for QualiMap                                                                                                        |

#### Output
Running the code will create the following folders in `outdir/prefix`:
- `readmap/`: folder with individual folders for each set of mapped reads (including the BAM, VCF, and consensus FASTA sequence)
    - `mapped_reads/qualimap/`: folder with the QualiMap output for each set of mapped reads
    - `summary.tsv`: file with the summary coverage for all mapped reads

### <a id="busco">Extracting BUSCO Loci</a>
In this step, we run correlation analysis to check for the extent of reference bias in BUSCO and assess if it changes the BUSCO tree topology. The parameters for this step is set in `2_busco_check/1_main.Rmd`.

| Parameters               | Definition                                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `codedir`                | Directory for folder `PhyloRBT/`                                                                                                      |
| `prefix`                 | Prefix for output files and folder                                                                                                    | 
| `outdir`                 | Output directory                                                                                                                      |
| `thread`                 | Number of threads for parallelisation                                                                                                 |
| `redo`                   | If `FALSE`, skip analysis if output files exist; if `TRUE`, overwrite previous results                                                |
| `file_refseq`            | Metadata file for reference assembly (e.g., `PhyloRBT/files/eucs_reference.txt`)                                                      |
| `file_shortreads`        | Metadata file for short reads (e.g., `PhyloRBT/files/eucs_shortreads.txt`)                                                            |
| `exe_busco`              | Executable for BUSCO                                                                                                                  |
| `exe_gff2bed`            | Executable for Gff2Bed from BEDOPS                                                                                                    |
| `exe_samtools`           | Executable for Samtools                                                                                                               |
| `exe_qualimap`           | Executable for QualiMap                                                                                                               |
| `exe_mafft`              | Executable for MAFFT                                                                                                                  |
| `exe_trimal`             | Executable for TrimAl                                                                                                                 |
| `exe_iqtree2`            | Executable for IQ-TREE2                                                                                                               |
| `exe_treeshrink`         | Executable for TreeShrink                                                                                                             |
| `busco_lineage`          | Lineage for BUSCO pipeline                                                                                                            |
| `busco_mode`             | Mode for BUSCO pipeline. Options: genome, transcriptome, or protein                                                                   |
| `thread_busco`           | Number of threads for BUSCO                                                                                                           |
| `min_busco_depth`        | Minimum BUSCO depth for mapped reads                                                                                                  |
| `run_trimal`             | Run `TrimAl` for individual locus alignments                                                                                          |
| `run_treeshrink`         | Run `Treeshrink` for individual gene trees                                                                                            |

#### Output
Running the code will create the following folders in `outdir/prefix`:
- `busco_lineage/`: folder with the lineage dataset for running BUSCO pipeline
- `busco_extraction/`
    - `refseq/`: folder with all BUSCO runs on individual reference genomes
        - `list_busco.txt`: list of single-copy, complete BUSCO loci that are shared between reference genomes
    - `short_reads/`: folder with all BUSCO runs on individual mapped reads
        - `metadata.tsv`: file with the error status for individual BUSCO loci for each set of mapped reads
    - `trees/`: folder with individual BUSCO alignments and trees comprising all reference genomes and mapped reads
        - `unfiltered/`: unfiltered locus alignments and trees from BUSCO pipeline
        - `trimal/`: TrimAl-filtered locus alignments and their associated trees (only if `run_trimal==TRUE`)
        - `unfiltered_treeshrink/`: Treeshrink-filtered gene trees (only if `run_treeshrink==TRUE`)

### <a id="runpipeline">Running Both Analyses</a>
In order to run both analyses on the same set of input files, you should update <a href="./config.yaml">`config.yaml`</a> and run the following command:
```
Rscript run_data_preparation.R --config config.yaml
Rscript run_data_preparation.R --config config.yaml --redo
```

---
## <a id="refs">References</a>
1. Danecek, P., et al. (<a href="https://doi.org/10.1093/gigascience/giab008">2021</a>). **Twelve years of SAMtools and BCFtools**. *GigaScience*, *10*(2), giab008.

2. Manni, M., et al. (<a href="https://doi.org/10.1002/cpz1.323">2021</a>). **BUSCO: Assessing Genomic Data Quality and Beyond**. *Current Protocols*, *1*(12), e323.

3. Vasimuddin, M., et al. (<a href="https://doi.org/10.1109/IPDPS.2019.00041">2019</a>). **Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems**. *IEEE Parallel and Distributed Processing Symposium*.

4. Neph, S., et al. (<a href="https://doi.org/10.1093/bioinformatics/bts277">2012</a>). **BEDOPS: high-performance genomic feature operations**. *Bioinformatics*, *28*(14), 1919-1920.

5. Katoh, K. & Standley, D.M. (<a href="https://doi.org/10.1093/molbev/mst010">2013</a>). **MAFFT multiple sequence alignment software version 7: Improvements in performance and usability**. *Molecular Biology and Evolution*, *30*(4), 772–780.

6. García-Alcalde, F., et al. (<a href="https://doi.org/10.1093/bioinformatics/bts503">2012</a>). **Qualimap: evaluating next-generation sequencing alignment data**. *Bioinformatics*, *28*(20), 2678-2679.

7. Mai, U. & Mirarab, S. (<a href="https://doi.org/10.1186/s12864-018-4620-2">2018</a>). **TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees**. *BMC Genomics*, *19*(272).

8. Capella-Gutiérrez, S., et al. (<a href="https://doi.org/10.1093/bioinformatics/btp348">2009</a>). **trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses**. *Bioinformatics*, *25*(15), 1972-1973.

9. Anthropic. (<a href="https://claude.ai/">2026</a>). Claude 4.6 Sonnet was used to generate `config.yaml` and `run_data_preparation.R`. 

---
*Last update: 28 July 2026 by Jeremias Ivan*