# PhyloRBT

**PhyloRBT (Phylogenetic Reference Bias Test)** is an R pipeline to assess reference bias in individual loci based on pairwise phylogenetic distances. It consists of two main steps: an optional data preparation and the reference bias test. It is mainly developed and tested using Linux, so there might be incompatibilities using MacOS and Windows.

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#inout">Input and Output Files</a>
- <a href="#genpipe">General Pipeline</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
PhyloRBT requires several software and R packages to run. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also use `install.packages()` built-in function in R or RStudio.

### Software
|    Name    |                                    Website                                     |                             Anaconda                             |
| ---------- |:------------------------------------------------------------------------------:|:----------------------------------------------------------------:|
| ASTRAL     | <a href="https://github.com/smirarab/ASTRAL">Link</a>                          | <a href="https://anaconda.org/bioconda/astral-tree">Link</a>     |
| IQ-TREE    | <a href="http://www.iqtree.org">Link</a>                                       | <a href="https://anaconda.org/bioconda/iqtree">Link</a>          |

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

## <a id="inout">Input and Output Files</a>

### Input Files
PhyloRBT requires a set of gene alignments and their corresponding trees as input. Each alignment should comprise different reconstructions of the same set of short-reads that are mapped to different reference genomes. For example, if we have short-reads data from `SpeciesA` and reference genomes from `SpeciesX`, `SpeciesY`, and `SpeciesZ`, each gene alignment should look like this:
```
> SpeciesA--SpeciesX
ATGACTAAACGTCACGGGTAGAAAAAA
> SpeciesA--SpeciesY
ATGACTAACCGTCACGGGTAGAAAAAA
> SpeciesA--SpeciesZ
ATGACTAATCGTCACGGGTAGAAAAAA
> SpeciesX
ATGACTAGGCGTCACGGGTAGAAAAAA
> SpeciesY
ATGACTACCCGTCACGGGTAGAAAAAA
> SpeciesZ
ATGACTATTCGTCACGGGTAGAAAAAA
```
where `SpeciesA--SpeciesX` refers to `SpeciesA` mapped to reference `SpeciesX`.

Each gene alignment and tree should be stored in individual folders, while the parent directory should be provided in <a href="./config.yaml">`config.yaml`</a>. For example:
```
parent_directory/
├── gene01/
│   ├── gene01.fa
│   └── gene01.fa.treefile
├── gene02/
│   ├── gene02.fa
│   └── gene02.fa.treefile
...
```

### Output Files
Running the code will create the following folders in `outdir/prefix`:
- `summary/`: folder with summary files from correlation analyses
    - `correlation_figs/`: folder with scatter plots from the correlation analysis
    - `prefix.cor.sumtable`: file with the results of the correlation analyses
    - `prefix.dist.tiff`: file with the phylogenetic distances of mapped reads that come from the same species
- `busco_tree/`: folder with output from the locus tree analyses
    - `all/`, `bias/`, and `nonbias/`: folder with individual analysis for different sets of loci
        - `alignment/`: folder with locus sequences
        - `per_reference/`: folder with one reconstruction per run
        - `trees`: folder with locus trees and ASTRAL tree

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
        conda install -c conda-forge -c bioconda r-data.table r-doSNOW r-log4r r-optparse r-phangorn r-rmarkdown r-tidyverse r-yaml bioconductor-biostrings astral-tree iqtree
        ```

3. **[Optional] Data Preparation**
    <br>If you have short-reads data and a set of reference genomes, you can refer to <a href="./data_preparation/">`data_preparation/`</a> to: (i) map each set of short-reads to individual reference genomes, and (ii) extract BUSCO loci from each reference genome and set of mapped reads. Please see the <a href="./data_preparation/README.md">README</a> file for more details.

    If you run `data_preparation`, the output folder that stores individual gene alignments and trees will be stored at `$outdir/busco_extraction/trees/unfiltered/`. If you set `$run_trimal==TRUE` or `$run_treeshrink==TRUE`, the post-trimming output files will be stored at `$outdir/busco_extraction/trees/trimal/` or `$outdir/busco_extraction/trees/unfiltered_treeshrink/`.

4. **Update the parameters in `config.yaml`** <br>

5. **Run PhyloRBT** <br>
    ```
    Rscript run_pipeline.R --config config.yaml
    Rscript run_pipeline.R --config config.yaml --redo
    ```

    In UNIX-based operating systems (e.g., Linux and MacOS), it is advisable to use `nohup` or `tmux` to run the whole pipeline. For Windows, you can use `psmux`. 

### <a id="refs">References</a>
1. Zhang, C., et al. (<a href="https://doi.org/10.1186/s12859-018-2129-y">2018</a>). **ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees**. *BMC Bioinformatics*, *19*(153).

2. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.

3. Pagès, H., et al. (<a href="https://doi.org/10.18129/B9.bioc.Biostrings">2026</a>). **Biostrings: Efficient manipulation of biological strings**. *R package*.

4. Barrett, T., et al. (<a href="https://doi.org/10.32614/CRAN.package.data.table">2026</a>). **data.table: Extension of 'data.frame'**. *R package*.

5. Daniel, F. (<a href="https://cran.r-project.org/package=doSNOW">2022</a>). **doSNOW: Foreach Parallel Adaptor for the 'snow' Package**. *R package*.

6. White, J.M., & Jacobs, A. (<a href="https://doi.org/10.32614/CRAN.package.log4r">2024</a>). **log4r: A Fast and Lightweight Logging System for R, Based on 'log4j'**. *R package*.

7. Davis, T.L. (<a href="https://doi.org/10.32614/CRAN.package.optparse">2026</a>). **optparse: Command Line Option Parser**. *R package*.

8. Schliep, K. (<a href="https://doi.org/10.1093/bioinformatics/btq706">2011</a>). **phangorn: Phylogenetic Analysis in R**. *Bioinformatics*, *27*(4), 592–593.

9. Allaire, J.J., et al. (<a href="https://doi.org/10.32614/CRAN.package.rmarkdown">2026</a>). **rmarkdown: Dynamic Documents for R**. *R package*.

10. Wickham, H., et al. (<a href="https://doi.org/10.21105/joss.01686">2019</a>). **Welcome to the tidyverse**. *Journal of Open Source Software*, *4*(43), 1686.

11. Stephens, J., et al. (<a href="https://doi.org/10.32614/CRAN.package.yaml">2025</a>). **yaml: Methods to Convert R Data to YAML and Back**. *R package*.

12. Anthropic. (<a href="https://claude.ai/">2026</a>). Claude 4.6 Sonnet was used to generate `config.yaml` and `run_pipeline.R`. 

---
*Last update: 24 July 2026 by Jeremias Ivan*