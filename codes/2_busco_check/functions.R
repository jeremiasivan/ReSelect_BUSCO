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
    cmd_samtools <- paste(exe_samtools, "view", nthread, "-b -u", fn_sam, "|", # convert to uncompressed (-u) BAM (-b)
                          exe_samtools, "collate", nthread, "-O -u - |", # group reads with the same name together, output as STDOUT (-O)
                          exe_samtools, "fixmate", nthread, "-m -u - - |", # correct flags used in the file, adding mate score tags (-m)
                          exe_samtools, "sort", nthread, "-u - |", # sort the reads based on their positions
                          exe_samtools, "markdup", nthread, "-", fn_bam) # mark duplicates based on the mate score tags
    system(cmd_samtools)

    # run bcftools mpileup
    cmd_bcftools <- paste(exe_bcftools, "mpileup", nthread, "-Ou -f", refseq, fn_bam, "|", # generate genotype likelihoods at each position with coverage
                          exe_bcftools, "call", nthread, "-Ou -mv |", # variant calling with default settings (-m) and output only variant sites (-v)
                          exe_bcftools, "view", nthread, "-i 'QUAL>20' |", # filter out variants with low quality score
                          exe_bcftools, "norm", nthread, "-f", refseq, "-Oz -o", fn_vcf) # normalize variants
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

# function: extract coordinates from FASTA header
f_extract_coordinates <- function(fasta_header, busco, prefix) {
    # remove > sign
    no_header <- unlist(strsplit(fasta_header, split=">"))[2]
    ls_header <- unlist(strsplit(no_header, split=":"))
    ls_coordinates <- unlist(strsplit(ls_header[2], split="-"))

    seq_name <- ls_header[1]
    first_coordinate <- as.numeric(ls_coordinates[1]) + 1
    second_coordinate <- as.numeric(ls_coordinates[2]) + 1

    # set the start and stop coordinates
    strand <- "+"
    start_coordinate <- NULL
    stop_coordinate <- NULL

    if (first_coordinate < second_coordinate) {
        start_coordinate <- first_coordinate
        stop_coordinate <- second_coordinate
    } else if (first_coordinate > second_coordinate) {
        start_coordinate <- second_coordinate
        stop_coordinate <- first_coordinate
        strand <- "-"
    } else {
        return(list(errmsg=paste0("Error: ", busco, " coordinates for ", prefix, ". Skipped.")))
    }

    return(list(seqname=seq_name, start=start_coordinate, stop=stop_coordinate, strand=strand))
}

# function: generate GFF file
f_create_gff <- function(coordinates, busco, fn_out) {
    df_gff <- data.table::data.table(
        seqname=rep(coordinates$seqname, 4),
        source=rep("Metaeuk", 4),
        feature=c("gene", "mRNA", "exon", "CDS"),
        start=rep(coordinates$start, 4),
        end=rep(coordinates$stop, 4),
        score=rep(".", 4),
        strand=rep(coordinates$strand, 4),
        frame=rep(".", 4),
        attribute=c(paste0("ID=",busco),
                    paste0("ID=",busco,"_mRNA;Parent=",busco),
                    paste0("ID=",busco,"_exon_0;Parent=",busco,"_mRNA"),
                    paste0("ID=",busco,"_CDS_0;Parent=",busco,"_exon_0"))
    )

    # save the new GFF file
    data.table::fwrite(df_gff, file=fn_out, sep="\t", quote=F, row.names=F, col.names=F)
}

# function: manipulate and save GFF file
f_manipulate_gff <- function(fn_input, coordinates, busco, prefix, fn_out) {
    # read GFF file
    df_gff <- data.table::fread(fn_input, header=FALSE, select=1:9)
    data.table::setnames(df_gff, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

    # convert attribute to match with gffread namingconvention
    df_gff$attribute <- gsub("TCS_ID", "ID", df_gff$attribute)

    # extract relevant GFF entry
    df_gff <- df_gff[df_gff$seqname==coordinates$seqname & df_gff$strand==coordinates$strand]
    if (nrow(df_gff) == 0) {
        return(list(errmsg=paste0("Error: ", busco, " GFF extraction for ", prefix, ". Skipped.")))
    }

    # extract the respective gene index
    gene_idx <- which(df_gff$feature == "gene" & df_gff$start == coordinates$start & df_gff$end == coordinates$stop)
    if (length(gene_idx) != 1) {
        # extract all CDS based on the coordinates
        df_gff_cds <- df_gff[df_gff$feature=="CDS" & df_gff$start>=coordinates$start & df_gff$end<=coordinates$stop]
        if (nrow(df_gff_cds) == 0) {
            return(list(errmsg=paste0("Error: ", busco, " gene extraction for ", prefix, ". Skipped.")))
        }

        # save the new GFF file
        data.table::fwrite(df_gff_cds, file=fn_out, sep="\t", quote=F, row.names=F, col.names=F)
        return(list(warnmsg=paste0("Warn: ", busco, " GFF file for ", prefix, " is CDS-only.")))
    }

    # extract all gene indices
    ls_gene_idx <- which(df_gff$feature == "gene")
    entry_idx <- match(gene_idx, ls_gene_idx)
    
    # subset the GFF table
    df_gff_subset <- NULL
    if (entry_idx == length(ls_gene_idx)) {
        df_gff_subset <- df_gff[gene_idx:nrow(df_gff)]
    } else {
        df_gff_subset <- df_gff[gene_idx:ls_gene_idx[entry_idx+1]-1]
    }
    
    # save the new GFF file
    data.table::fwrite(df_gff_subset, file=fn_out, sep="\t", quote=F, row.names=F, col.names=F)
}

# function: extract all BUSCO alignments from GFF
f_extract_fasta_from_gff <- function(fn_input, fn_gff, fn_cds_out, fn_concat_out, exe_gffread){
    # run gffread
    cmd_gffread <- paste(exe_gffread, "-g", fn_input, "-x", fn_cds_out, fn_gff)
    system(cmd_gffread)

    # open the CDS sequences
    all_cds <- Biostrings::readDNAStringSet(fn_cds_out)
    all_headers <- stringr::str_sort(names(all_cds), numeric=T)

    # concat all CDS sequences
    header <- ""
    seq <- ""
    for (i in all_headers) {
        ls_header <- unlist(strsplit(i, split="\\|"))

        if (header == "") {
            header <- paste0(header, ls_header[length(ls_header)])
        } else {
            header <- paste0(header, "|", ls_header[length(ls_header)])
        }
        
        seq <- paste0(seq, all_cds[[i]])
    }

    # save the concatenated FASTA in a file
    concat_cds <- Biostrings::DNAStringSet(seq)
    names(concat_cds) <- header
    Biostrings::writeXStringSet(concat_cds, filepath=fn_concat_out)
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

# function: run IQ-Tree 2
f_iqtree2 <- function(fn_input, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-bb 1000",
                         "-T 1 --quiet -redo")
    system(cmd_iqtree2)
}