# PhyloRBT

**PhyloRBT (Phylogenetic Reference Bias Test)** is an R pipeline to assess reference bias in individual loci based on pairwise phylogenetic distances. It consists of two main steps: an optional data preparation and reference bias checking. It is mainly developed and tested using Linux, so there might be incompatibilities using MacOS and Windows.

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
PhyloRBT requires several software and R packages to run. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also use `install.packages()` built-in function in R or RStudio.

### Software
|    Name    |                                    Website                                     |                             Anaconda                             |
| ---------- |:------------------------------------------------------------------------------:|:----------------------------------------------------------------:|
| ASTRAL     | <a href="https://github.com/smirarab/ASTRAL">Link</a>                          | <a href="https://anaconda.org/bioconda/astral-tree">Link</a>     |
| BCFtools   | <a href="https://github.com/samtools/bcftools">Link</a>                        | <a href="https://anaconda.org/bioconda/bcftools">Link</a>        |
| BUSCO      | <a href="https://busco.ezlab.org">Link</a>                                     | <a href="https://anaconda.org/bioconda/busco">Link</a>           |
| BWA-MEM2   | <a href="https://github.com/bwa-mem2/bwa-mem2">Link</a>                        | <a href="https://anaconda.org/bioconda/bwa-mem2">Link</a>        |
| Gff2BED    | <a href="https://bedops.readthedocs.io/en/latest/index.html">Link</a>          | <a href="https://anaconda.org/bioconda/gff2bed">Link</a>         |
| IQ-TREE    | <a href="http://www.iqtree.org">Link</a>                                       | <a href="https://anaconda.org/bioconda/iqtree">Link</a>          |
| MAFFT      | <a href="https://github.com/GSLBiotech/mafft">Link</a>                         | <a href="https://anaconda.org/bioconda/mafft">Link</a>           |
| QualiMap   | <a href="http://qualimap.conesalab.org/">Link</a>                              | <a href="https://anaconda.org/bioconda/qualimap">Link</a>        |
| SAMtools   | <a href="https://github.com/samtools/samtools">Link</a>                        | <a href="https://anaconda.org/bioconda/samtools">Link</a>        |
| TreeShrink | <a href="https://github.com/uym2/TreeShrink">Link</a>                          | <a href="https://anaconda.org/bioconda/treeshrink">Link</a>      |
| TrimAl     | <a href="https://github.com/inab/trimal">Link</a>                              | <a href="https://anaconda.org/bioconda/trimal">Link</a>          |

### R packages
|    Name      |                                 CRAN / Bioconductor                                     |                                   Anaconda                               |
| ------------ |:---------------------------------------------------------------------------------------:|:------------------------------------------------------------------------:|
| Biostrings   | <a href="https://bioconductor.org/packages/Biostrings">Link</a>                         | <a href="https://anaconda.org/bioconda/bioconductor-biostrings">Link</a> |
| data.table   | <a href="https://cran.r-project.org/package=data.table">Link</a>                        | <a href="https://anaconda.org/conda-forge/r-data.table">Link</a>         |
| doSNOW       | <a href="https://cran.r-project.org/package=doSNOW">Link</a>                            | <a href="https://anaconda.org/conda-forge/r-dosnow">Link</a>             |
| log4r        | <a href="https://cran.r-project.org/package=log4r">Link</a>                             | <a href="https://anaconda.org/conda-forge/r-log4r">Link</a>              |
| optparse     | <a href="https://cran.r-project.org/package=optparse">Link</a>                          | <a href="https://anaconda.org/conda-forge/r-optparse">Link</a>           |
| phangorn     | <a href="https://cran.r-project.org/package=phangorn">Link</a>                          | <a href="https://anaconda.org/conda-forge/r-phangorn">Link</a>           |
| rmarkdown    | <a href="https://cran.r-project.org/package=rmarkdown">Link</a>                         | <a href="https://anaconda.org/conda-forge/r-rmarkdown">Link</a>          |
| tidyverse    | <a href="https://cran.r-project.org/package=tidyverse">Link</a>                         | <a href="https://anaconda.org/conda-forge/r-tidyverse">Link</a>          |
| yaml         | <a href="https://cran.r-project.org/package=yaml">Link</a>                              | <a href="https://anaconda.org/conda-forge/r-yaml">Link</a>               |

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/PhyloRBT.git
    ```

2. **Install the prerequisites** <br>
    - Create a new conda environment
        ```
        conda create -n phylorbt
        conda activate phylorbt
        ```
    - Installing prerequisites
        ```
        conda install -c conda-forge -c bioconda r-data.table r-doSNOW r-log4r r-optparse r-phangorn r-rmarkdown r-tidyverse r-yaml bioconductor-biostrings astral-tree bcftools busco bwa-mem2 gff2bed iqtree mafft qualimap samtools treeshrink trimal
        ```
    
    If `bcftools` returns a `libgsl.so.25` error, you can either download the software <a href="https://www.htslib.org/download/">here</a>, or try to set the `conda` channel priorities before installing any package:
    ```
    conda config --prepend channels r
    conda config --prepend channels bioconda
    conda config --prepend channels conda-forge
    conda config --set channel_priority strict
    ```

3. **Update the parameters in `config.yaml`** <br>

4. **Run PhyloRBT** <br>
    ```
    Rscript run_pipeline.R --config config.yaml
    Rscript run_pipeline.R --config config.yaml --redo
    ```

    In UNIX-based operating systems (e.g., Linux and MacOS), it is advisable to use `nohup` or `tmux` to run the whole pipeline. For Windows, you can use `psmux`. 

---
*Last update: 23 July 2026 by Jeremias Ivan*