# functions for codes/3_read_mapping

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