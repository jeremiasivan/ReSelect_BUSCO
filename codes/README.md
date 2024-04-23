# BusIER

## Table of Content
- <a href="#analyses">Analyses</a>
    - <a href="#download">Downloading Raw Data</a>
    - <a href="#check">Checking Reference Bias on BUSCOs</a>

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

### <a id="check">Checking Reference Bias on BUSCOs</a>
In this step, we run correlation analysis to check for the presence of reference bias in BUSCOs. The parameters for this step is set in `2_busco_check/1_main.Rmd`.

*to be updated*

---
*Last update: 23 April 2024 by Jeremias Ivan*