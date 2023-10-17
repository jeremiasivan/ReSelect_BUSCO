# functions for codes/1_data_download

# function: download references FASTA sequence using NCBI Datasets
f_refseq_download <- function(dir_datasets, accession, file_output) {
    cmd <- paste(dir_datasets, "download genome",
                "accession", accession,
                "--filename", file_output,
                "--include genome --no-progressbar")

    system(cmd)
    system(paste("unzip", file_output))
}

