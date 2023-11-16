# functions for codes/3_read_mapping

# function: read map fastq to reference sequence using BWA-MEM2
f_read_map <- function(refseq, fastq, thread, exe_bwamem2, file_sam) {
    # index reference file
    cmd_index <- paste(exe_bwamem2, "index", refseq)
    system(cmd_index)

    # read-map
    cmd_readmap <- paste(exe_bwamem2, "mem -t", thread, refseq)

    len_fastq <- length(fastq)
    if (len_fastq == 1) {
        cmd_readmap <- paste(cmd_readmap, fastq[1], ">", file_sam)
    } else if (len_fastq == 2) {
        cmd_readmap <- paste(cmd_readmap, fastq[1], fastq[2], ">", file_sam)
    }
    system(cmd_readmap)
}

# function: combine multiple FASTA headers into concatenated sequence
f_combine_fasta_headers <- function(header, file_fasta, file_output) {
    # initialize variable
    combined_fasta <- paste0('>', header, "\n")

    # initiate connection
    con <- file(file_fasta, open = "r")

    # iterate through lines
    while (length(line <- readLines(con, n=1)) > 0) {
        if (!grepl("^>", line)) {
            combined_fasta <- paste0(combined_fasta, line, "\n")
        }
    }

    # close connection
    close(con)

    # write output file
    writeLines(combined_fasta, file_output)
}

# function: convert SAM to BAM and FASTA
f_convert_SAM <- function(prefix, dir_output, thread, exe_samtools) {
    fn_sam <- paste0(dir_output, "/", prefix, ".sam")
    fn_bam_fixmate <- paste0(dir_output, "/", prefix, ".fixmate.bam")
    fn_bam_sort <- paste0(dir_output, "/", prefix, ".sort.bam")
    fn_bam_markdup <- paste0(dir_output, "/", prefix, ".markdup.bam")
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")
    fn_fasta_one <- paste0(dir_output, "/", prefix, ".1.fa")
    fn_fasta_two <- paste0(dir_output, "/", prefix, ".2.fa")

    # run samtools fixmate
    cmd_fixmate <- paste(exe_samtools, "fixmate",
                         "--threads", thread,
                         "-O bam,level=1",
                         "-m", fn_sam, fn_bam_fixmate)
    system(cmd_fixmate)

    # run samtools sort
    cmd_sort <- paste(exe_samtools, "sort",
                      "--threads", thread,
                      "-l 1",
                      "-o", fn_bam_sort, fn_bam_fixmate)
    system(cmd_sort)                  
    
    # run samtools markdup
    cmd_markdup <- paste(exe_samtools, "markdup",
                         "--threads", thread,
                         "-O bam,level=1",
                         fn_bam_sort, fn_bam_markdup) 
    system(cmd_markdup)

    # run samtools view to BAM
    cmd_view <- paste(exe_samtools, "view",
                      "--threads", thread,
                      fn_bam_markdup, 
                      "-o", fn_bam)
    system(cmd_view)

    # run samtools fasta
    cmd_fasta <- paste(exe_samtools, "fasta",
                       "--threads", thread,
                       "-1", fn_fasta_one, "-2", fn_fasta_two,
                       fn_bam)                                          
    system(cmd_fasta)

    # run FASTA header concatenation
    fn_fasta_concat_one <- paste0(dir_output, "/", prefix, ".1.fa")
    fn_fasta_concat_two <- paste0(dir_output, "/", prefix, ".2.fa")
    f_combine_fasta_headers(prefix, fn_fasta_one, fn_fasta_concat_one)
    f_combine_fasta_headers(prefix, fn_fasta_two, fn_fasta_concat_two)
}

#f_index_reference <- function(refseq, exe_samtools) {
#    cmd_index <- paste(exe_samtools, "faidx", refseq)
#    system(cmd_index)
#}

#f_variant_call <- function(prefix, refseq, bam, chromosome, dir_vcf, exe_bcftools){
#    fn_bcf <- paste0(dir_vcf, "/", prefix, ".genotype_likelihood.bcf")
#    
#    cmd_bcftools <- paste(exe_bcftools, "mpileup -f", refseq)
#
#    if (chromosome != "" && !is.null(chromosome)) {
#        cmd_bcftools <- paste(cmd_bcftools, "-r", chromosome)
#    }
#    
#    cmd_bcftools <- paste(cmd_bcftools, fn_bam, "|",
#                          exe_bcftools, "call",
#                          "-m -Oz -f GQ",
#                          "-o", fn_output)
#
#    system(cmd_bcftools)
#}