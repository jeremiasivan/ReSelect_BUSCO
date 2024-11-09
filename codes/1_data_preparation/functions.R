# functions for codes/1_data_download

# function: add entry to logfile
f_write_log <- function(fn_log, msg) {
    write.table(msg, file=fn_log, quote=F, row.names=F, col.names=F, append=T)
}

# function: download references FASTA sequence using NCBI Datasets
f_refseq_download <- function(dir_datasets, accession, file_output) {
    cmd <- paste(dir_datasets, "download genome",
                "accession", accession,
                "--filename", file_output,
                "--include genome --no-progressbar")
    system(cmd)
}

# function: download short reads FASTQ sequence using SRA Toolkit
f_shortreads_download <- function(dir_sratoolkit, accession, dir_output) {
    # prefetch from NCBI
    exe_prefetch <- paste0(dir_sratoolkit,"/prefetch")
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
    exe_fastqdump <- paste0(dir_sratoolkit,"/fastq-dump")
    cmd_fastq <- paste(exe_fastqdump, "--outdir", dir_fastq, "--split-files", file_sra)
    system(cmd_fastq)
}

# function: quality control for short reads
f_qc_short_reads <- function(fastq_one, fastq_two, fn_adapters, prefix, min_quality, thread, exe_adapterremoval) {
    cmd_qc <- paste(exe_adapterremoval, "--file1", fastq_one)
    
    # check the reverse fastq file
    if (!is.null(fastq_two)) {
        cmd_qc <- paste(cmd_qc, "--file2", fastq_two)
    }

    # update the adapters
    if (fn_adapters != "" && file.exists(fn_adapters)) {
        cmd_qc <- paste(cmd_qc, "--adapter-list", fn_adapters)
    }

    # set the minimum quality score and merge overlapping reads
    cmd_qc <- paste(cmd_qc, "--basename", prefix,
                    "--trimqualities --trimns --minquality", min_quality)
    
    # run AdapterRemoval to get prefix.collapsed.truncated
    system(cmd_qc)
}

# function: readmap fastq to reference sequence using BWA-MEM2
f_read_mapping <- function(refseq, fastq_one, fastq_two, thread, exe_bwamem2, file_sam) {
    # index reference file
    cmd_index <- paste(exe_bwamem2, "index", refseq)
    system(cmd_index)

    # readmap
    cmd_readmap <- paste(exe_bwamem2, "mem -t", thread, refseq)

    if (is.null(fastq_two)) {
        cmd_readmap <- paste(cmd_readmap, fastq_one, ">", file_sam)
    } else {
        cmd_readmap <- paste(cmd_readmap, fastq_one, fastq_two, ">", file_sam)
    }
    system(cmd_readmap)
}

# function: convert SAM to BAM
f_sam_to_bam <- function(prefix, dir_output, thread, exe_samtools) {
    # initiate variables
    fn_sam <- paste0(dir_output, "/", prefix, ".sam")
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")

    # set the number of threads
    nthread <- paste("--threads", thread)

    # run samtools
    cmd_samtools <- paste(exe_samtools, "view", nthread, "-b -u", fn_sam, "|", # convert to uncompressed (-u) BAM (-b)
                          exe_samtools, "collate", nthread, "-O -u - |",       # group reads with the same name together, output as STDOUT (-O)
                          exe_samtools, "fixmate", nthread, "-m -u - - |",     # correct flags used in the file, adding mate score tags (-m)
                          exe_samtools, "sort", nthread, "-u - |",             # sort the reads based on their positions
                          exe_samtools, "markdup -r", nthread, "-", fn_bam)    # remove duplicates based on the mate score tags
    system(cmd_samtools)

    # index BAM file
    cmd_bam_index <- paste(exe_samtools, "index", nthread, fn_bam)
    system(cmd_bam_index)
}

# function: variant calling
f_variant_calling <- function(prefix, dir_output, thread, refseq, exe_bcftools) {
    # initiate variables
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")
    fn_vcf <- paste0(dir_output, "/", prefix, ".vcf.gz")
    fn_fas <- paste0(dir_output, "/", prefix, ".fa")

    nthread <- paste("--threads", thread)

    # run bcftools mpileup
    cmd_bcftools <- paste(exe_bcftools, "mpileup", nthread, "-Ou -f", refseq, fn_bam, "|",      # generate genotype likelihoods at each position with coverage
                          exe_bcftools, "call", nthread, "-Ou -mv |",                           # variant calling with default settings (-m) and output only variant sites (-v)
                          exe_bcftools, "view", nthread, "-V indels -i 'QUAL>15 & MQ>30' |",    # filter out variants with low quality score
                          exe_bcftools, "norm", nthread, "-f", refseq, "-Oz -o", fn_vcf)        # normalize variants
    system(cmd_bcftools)

    # index VCF file
    cmd_vcf_index <- paste(exe_bcftools, "index -t", nthread, fn_vcf)
    system(cmd_vcf_index)

    # generate consensus sequence
    cmd_consensus <- paste("cat", refseq, "|", exe_bcftools, "consensus", fn_vcf, ">", fn_fas)
    system(cmd_consensus)
}
