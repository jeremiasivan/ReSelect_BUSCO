# ReSelect BUSCO

**ReSelect BUSCO (Reference Selection on BUSCO Inference)** is an R pipeline to assess the extent of bias introduced by reference selection in BUSCO inference. It consists of two main steps: data preparation and reference bias checking. It is mainly developed and tested using MacOS and Linux, so there might be incompatibilities using Windows.

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>

## <a id="prereqs">Prerequisites</a>
This pipeline requires several software and R packages to run. All software have to be executable (except for **SRA Toolkit** which points to `bin` directory), while the R packages should be installed either in your local directory or virtual environment. We recommend you to use environment management system (e.g. `conda`) to install the software and packages, but you can also use `install.packages()` built-in function in R or RStudio.

### Software
| Name                    | Website / Github                                                      |
| ----------------------- |:---------------------------------------------------------------------:|
| AdapterRemoval          | <a href="https://adapterremoval.readthedocs.io/en/stable/#">Link</a>  |
| ASTRAL                  | <a href="https://github.com/smirarab/ASTRAL">Link</a>                 |
| Bcftools                | <a href="https://github.com/samtools/bcftools">Link</a>               |
| BUSCO                   | <a href="https://busco.ezlab.org">Link</a>                            |
| BWA-MEM2                | <a href="https://github.com/bwa-mem2/bwa-mem2">Link</a>               |
| Gff2Bed (BEDOPS)        | <a href="https://bedops.readthedocs.io/en/latest/index.html">Link</a> |
| GffRead                 | <a href="https://github.com/gpertea/gffread">Link</a>                 |
| IQ-TREE 2               | <a href="http://www.iqtree.org">Link</a>                              |
| MAFFT                   | <a href="https://mafft.cbrc.jp/alignment/software/">Link</a>          |
| NCBI Datasets           | <a href="https://github.com/ncbi/datasets">Link</a>                   |
| QualiMap                | <a href="http://qualimap.conesalab.org/">Link</a>                     |
| Samtools                | <a href="https://github.com/samtools/samtools">Link</a>               |
| SRA Toolkit             | <a href="https://github.com/ncbi/sra-tools">Link</a>                  |
| TreeShrink              | <a href="https://github.com/uym2/TreeShrink">Link</a>                 |

### R packages
|    Name      |                                 CRAN / Bioconductor                                     |                                   Anaconda                               |
| ------------ |:---------------------------------------------------------------------------------------:|:------------------------------------------------------------------------:|
| Biostrings   | <a href="https://www.bioconductor.org/packages//2.7/bioc/html/Biostrings.html">Link</a> | <a href="https://anaconda.org/bioconda/bioconductor-biostrings">Link</a> |
| data.table   | <a href="https://cran.r-project.org/package=data.table">Link</a>                        | <a href="https://anaconda.org/conda-forge/r-data.table">Link</a>         |
| doSNOW       | <a href="https://cran.r-project.org/package=doSNOW">Link</a>                            | <a href="https://anaconda.org/conda-forge/r-dosnow">Link</a>             |
| log4r        | <a href="https://cran.r-project.org/package=log4r">Link</a>                             | <a href="https://anaconda.org/conda-forge/r-log4r">Link</a>              |
| phangorn     | <a href="https://cran.r-project.org/package=phangorn">Link</a>                          | <a href="https://anaconda.org/conda-forge/r-phangorn">Link</a>           |
| rmarkdown    | <a href="https://cran.r-project.org/package=rmarkdown">Link</a>                         | <a href="https://anaconda.org/conda-forge/r-rmarkdown">Link</a>          |
| tidyverse    | <a href="https://cran.r-project.org/package=tidyverse">Link</a>                         | <a href="https://anaconda.org/conda-forge/r-tidyverse">Link</a>          |

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/BusIER.git
    ```

2. **Install the prerequisites** <br>
    Please download the required software from the links above. For the `R` packages, I prefer to download them from `Anaconda` as below.

    - Setting up conda environment with R
        ```
        conda create -n busier
        conda activate busier
        ```
    -  Installing R packages
        ```
        conda install package-name
        ```
        Notes: Please install all of the R packages and their dependencies. A good starting point is to install <a href="https://anaconda.org/conda-forge/r-essentials">`r-essentials`</a> which includes commonly-used packages in R. 

3. **Update the parameters in the file** <br>

4. **Run the code file** <br>
    For running individual steps:
    ```
    Rscript -e "rmarkdown::render('~/SimNOW/codes/1_data_download/1_main.Rmd')"
    Rscript -e "rmarkdown::render('~/SimNOW/codes/2_busco_check/1_main.Rmd')"
    ```

    For running the whole pipeline:
    ```
    Rscript ~/SimNOW/codes/run_all.R
    ```

    In UNIX-based operating systems (e.g., Linux and MacOS), it is advisable to use `nohup` or `tmux` to run the whole pipeline. For Windows, you can use `start`, but I have never tried it before. 

---
*Last update: 04 September 2024 by Jeremias Ivan*