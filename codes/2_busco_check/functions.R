# functions for codes/2_busco_check

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
f_convert_sam2bam <- function(prefix, dir_output, thread, exe_samtools) {
    fn_sam <- paste0(dir_output, "/", prefix, ".sam")
    fn_bam_fixmate <- paste0(dir_output, "/", prefix, ".fixmate.bam")
    fn_bam_sort <- paste0(dir_output, "/", prefix, ".sort.bam")
    fn_bam_markdup <- paste0(dir_output, "/", prefix, ".markdup.bam")
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")

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
}

# function: run BUSCO pipeline
f_run_busco <- function(fn_fasta, lineage, prefix, dir_output, mode, thread, exe_busco) {
    cmd_busco <- paste(exe_busco,
                       "--offline",
                       "-i", fn_fasta,
                       "-l", lineage,
                       "-m", mode,
                       "-o", prefix,
                       "--out_path", dir_output,
                       "-c", thread,
                       "--quiet --force")
    system(cmd_busco)
}

# functions: extract BUSCO region
f_extract_busco_from_BAM <- function(fn_bam, fn_out_bam, coordinates, fn_out_1, fn_out_2, exe_samtools) {
    # run samtools view
    cmd_view <- paste(exe_samtools, "view",
                      "-o", fn_out_bam,
                      "-b", fn_bam,
                      coordinates)
    
    # check if there is an error message
    out_msg <- system(cmd_view, intern=TRUE)
    if (length(grep(".+specifies an invalid region or unknown reference*", out_msg)) != 0) {
        unlink(fn_out_bam)
        return(NULL)
    }

    # run samtools fasta
    cmd_fasta <- paste(exe_samtools, "fasta",
                       "-1", fn_out_1, "-2", fn_out_2,
                       fn_out_bam)                                          
    system(cmd_fasta)
}

# functions: combine individual FASTA as MSA
f_fasta2msa <- function(fn_input, header, fn_out) {
    # initiate variable
    first_sequence <- TRUE

    # open the FASTA file
    con <- file(fn_input, "r")

    # iterate over lines
    while (length(line <- readLines(con, n = 1)) > 0) {
        if (grepl("^>+", line)) {
            if (first_sequence) {
                write.table(paste0(">", header), file=fn_out, quote=F, row.names=F, col.names=F, append=T)
                first_sequence <- FALSE
            }
        } else {
            write.table(line, file=fn_out, quote=F, row.names=F, col.names=F, append=T)
        }
    }
    
    # close the file connection
    close(con)
}

# functions: run IQ-Tree 2
f_iqtree2 <- function(fn_input, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-T 1 --quiet -redo")
    system(cmd_iqtree2)
}