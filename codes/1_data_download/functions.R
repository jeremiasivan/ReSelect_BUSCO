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

# function: download short reads FASTQ sequence using SRA Toolkit
f_shortreads_download <- function(dir_sratoolkit, accession, dir_output) {
    # prefetch from NCBI
    exe_prefetch <- paste0(dir_sratoolkit,"/bin/prefetch")
    cmd_download <- paste(exe_prefetch, "--output-directory", dir_output, accession)
    system(cmd_download)

    # output SRA file
    file_sra <- paste0(dir_output,"/",accession,"/",accession,".sra")
    
    # create fastq folder
    dir_fastq <- paste0(dir_output,"/",accession,"/fastq/")
    if (!dir.exists(dir_fastq)) {
        dir.create(dir_fastq, recursive=T)
    }

    # download fastq files
    exe_fastqdump <- paste0(dir_sratoolkit,"/bin/fastq-dump")
    cmd_fastq <- paste(exe_fastqdump, "--outdir", dir_fastq, "--split-files", file_sra)
    system(cmd_fastq)
}