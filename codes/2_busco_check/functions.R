# functions for codes/2_busco_check

# function: quality control for short reads
f_qc_short_reads <- function(fastq, fn_adapters, prefix, min_quality, thread, exe_adapterremoval) {
    cmd_qc <- paste(exe_adapterremoval, "--file1", fastq[1])
    
    # check the number of fastq files
    if (length(fastq) == 2) {
        cmd_qc <- paste(cmd_qc, "--file2", fastq[2])
    }

    # update the adapters
    if (fn_adapters != "" && file.exists(fn_adapters)) {
        cmd_qc <- paste(cmd_qc, "--adapter-list", fn_adapters)
    }

    # set the minimum quality score and merge overlapping reads
    cmd_qc <- paste(cmd_qc, "--basename", prefix,
                    "--trimqualities --minquality", min_quality,
                    "--collapse")
    
    # run AdapterRemoval to get prefix.collapsed.truncated
    system(cmd_qc)
}

# function: read map fastq to reference sequence using BWA-MEM2
f_read_mapping <- function(refseq, fastq, thread, exe_bwamem2, file_sam) {
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

# function: convert SAM to BAM and VCF
f_variant_calling <- function(prefix, dir_output, thread, refseq, exe_samtools, exe_bcftools) {
    # initiate variables
    fn_sam <- paste0(dir_output, "/", prefix, ".sam")
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")
    fn_vcf <- paste0(dir_output, "/", prefix, ".vcf.gz")
    fn_fas <- paste0(dir_output, "/", prefix, ".fa")

    nthread <- paste("--threads", thread)

    # run samtools
    cmd_samtools <- paste(exe_samtools, "view", nthread, "-b -u", fn_sam, "|",
                          exe_samtools, "collate", nthread, "-O -u - |",
                          exe_samtools, "fixmate", nthread, "-m -u - - |",
                          exe_samtools, "sort", nthread, "-u - |",
                          exe_samtools, "markdup", nthread, "-", fn_bam)
    system(cmd_samtools)

    # run bcftools mpileup
    cmd_bcftools <- paste(exe_bcftools, "mpileup", nthread, "-Ou -f", refseq, fn_bam, "|",
                          exe_bcftools, "call", nthread, "-Ou -mv |",
                          exe_bcftools, "view", nthread, "-i 'QUAL>20' |",
                          exe_bcftools, "norm", nthread, "-f", refseq, "-Oz -o", fn_vcf) 
    system(cmd_bcftools)

    # index VCF file
    cmd_vcf_index <- paste(exe_bcftools, "index -t", nthread, fn_vcf)
    system(cmd_vcf_index)

    # generate consensus sequence
    cmd_consensus <- paste("cat", refseq, "|", exe_bcftools, "consensus", fn_vcf, ">", fn_fas)
    system(cmd_consensus)
}

# function: run BUSCO pipeline
f_run_busco <- function(fn_fasta, lineage, prefix, dir_output, mode, thread, exe_busco) {
    cmd_busco <- paste(exe_busco,
                       "-i", fn_fasta,
                       "-l", lineage,
                       "-m", mode,
                       "-o", prefix,
                       "--download_path", dir_output,
                       "--out_path", dir_output,
                       "-c", thread,
                       "--quiet --force --offline")
    system(cmd_busco)
}

# functions: extract BUSCO region
f_extract_busco <- function(fn_vcf, refseq, busco_header, fn_out, exe_samtools, exe_bcftools) {
    # set up variable
    is_reverse <- FALSE

    # remove > sign
    no_header <- unlist(strsplit(busco_header, split=">"))[2]

    # extract start and stop coordinates
    ls_coordinates <- unlist(strsplit(no_header, split=":"))
    start_coordinate <- as.numeric(unlist(strsplit(ls_coordinates[2], split="-"))[1]) + 1
    stop_coordinate <- as.numeric(unlist(strsplit(ls_coordinates[2], split="-"))[2]) + 1

    # update the coordinates
    coordinates <- paste0(ls_coordinates[1], ":", start_coordinate, "-", stop_coordinate)
    if (start_coordinate > stop_coordinate) {
        coordinates <- paste0(ls_coordinates[1], ":", stop_coordinate, "-", start_coordinate)
        is_reverse <- TRUE
    }

    # extract consensus BUSCO sequence
    cmd_consensus <- paste(exe_samtools, "faidx", refseq, coordinates, "|",
                           exe_bcftools, "consensus", fn_vcf, "-o", fn_out)
    if (is_reverse) {
        cmd_consensus <- paste0(cmd_consensus, "_temp")
    }                                   
    system(cmd_consensus)

    # reverse complement the sequence
    if (is_reverse) {
        fasta_sequence <- Biostrings::readDNAStringSet(paste0(fn_out, "_temp"))
        rc_fasta_sequence <- Biostrings::reverseComplement(fasta_sequence)
        Biostrings::writeXStringSet(rc_fasta_sequence, filepath=fn_out)
        unlink(paste0(fn_out, "_temp"))
    }
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

# functions: run MAFFT
f_mafft <- function(fn_input, fn_output, params_mafft, exe_mafft) {
    cmd_mafft <- paste(exe_mafft, params_mafft,
                       fn_input, ">", fn_output)
    system(cmd_mafft)
}

# functions: run IQ-Tree 2
f_iqtree2 <- function(fn_input, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-bb 1000",
                         "-T 1 --quiet -redo")
    system(cmd_iqtree2)
}