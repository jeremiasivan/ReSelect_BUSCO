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

# function: convert SAM to BAM and FASTA
f_convert_SAM <- function(prefix, dir_sam, dir_bam, dir_fasta, thread, exe_samtools) {
    samtools_thread <- paste0("@", thread)

    fn_sam <- paste0(dir_sam, "/", prefix, ".sam")
    fn_bam_fixmate <- paste0(dir_bam, "/", prefix, ".fixmate.bam")
    fn_bam_sort <- paste0(dir_bam, "/", prefix, ".sort.bam")
    fn_bam_markdup <- paste0(dir_bam, "/", prefix, ".markdup.bam")
    fn_bam <- paste0(dir_bam, "/", prefix, ".bam")
    fn_fasta <- paste0(dir_fasta, "/", prefix, ".fa")

    # run samtools fixmate
    cmd_fixmate <- paste(exe_samtools, "fixmate",
                         samtools_thread,
                         "-O bam",
                         "-m", fn_sam, fn_bam_fixmate)
    system(cmd_fixmate)

    # run samtools sort
    cmd_sort <- paste(exe_samtools, "sort",
                      samtools_thread,
                      "-o", fn_bam_sort, fn_bam_fixmate)
    system(cmd_sort)                  
    
    # run samtools markdup
    cmd_markdup <- paste(exe_samtools, "markdup",
                         samtools_thread,
                         "-O bam",
                         fn_bam_sort, fn_bam_markdup) 
    system(cmd_markdup)

    # run samtools view to BAM
    cmd_view <- paste(exe_samtools, "view",
                      samtools_thread,
                      fn_bam_markdup, fn_bam)
    system(cmd_view)

    # run samtools fasta
    cmd_fasta <- paste(exe_samtools, "fasta",
                       samtools_thread,
                       fn_bam, ">", fn_fasta)                                          
    system(cmd_fasta)
}