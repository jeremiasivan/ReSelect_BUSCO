# functions for codes/2_busco_check

# function: run BUSCO pipeline
f_run_busco <- function(fn_fasta, lineage, prefix, dir_output, mode, thread, exe_busco) {
    cmd_busco <- paste(exe_busco,
                       "-i", fn_fasta,
                       "-l", lineage,
                       "-m", mode,
                       "--metaeuk",
                       "-o", prefix,
                       "--download_path", dir_output,
                       "--out_path", dir_output,
                       "-c", thread,
                       "--quiet --force --offline")
    system(cmd_busco)
}

# function: check coverage
f_calculate_read_coverage <- function(fn_bam, fn_bed, exe_samtools) {
    # retrieve average coverage
    cmd_coverage <- paste(exe_samtools, "depth",
                          "-b", fn_bed,
                          fn_bam)
    ls_output <- system(cmd_coverage, intern=T)
    ls_coverage <- sapply(ls_output, function(x){ strsplit(x, split="\t")[[1]][3] })

    # return average coverage
    return(round(mean(as.numeric(ls_coverage)),3))
}

# function: combine individual FASTA as MSA
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

# function: run MAFFT
f_mafft <- function(fn_input, fn_output, params_mafft, exe_mafft) {
    cmd_mafft <- paste(exe_mafft, params_mafft,
                       fn_input, ">", fn_output)
    system(cmd_mafft)
}

# function: run TrimAl
f_trimal <- function(fn_input, fn_output, params_trimal, exe_trimal) {
    cmd_trimal <- paste(exe_trimal,
                        "-in", fn_input,
                        "-out", fn_output,
                        params_trimal)
    system(cmd_trimal)
}

# function: run IQ-Tree 2
f_iqtree2 <- function(fn_input, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-T 1 --quiet -redo")

    # check if bootstrap is required
    seq <- Biostrings::readBStringSet(fn_input)
    if (length(unique(as.character(seq))) > 3) {
        cmd_iqtree2 <- paste(cmd_iqtree2, "-bb 1000")
    }
    
    system(cmd_iqtree2)
}

# function: run TreeShrink
f_treeshrink <- function(fn_input, prefix, dir_output, exe_treeshrink) {
    cmd_treeshrink <- paste(exe_treeshrink,
                            "-t", fn_input,
                            "-O", prefix,
                            "-o", dir_output)
    system(cmd_treeshrink)
}
