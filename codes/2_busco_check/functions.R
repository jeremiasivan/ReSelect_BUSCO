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

# run MNTD on BUSCOs
f_run_mntd <- function(ls_busco, ls_refseq, is_ref_included, dir_busco_tree, ls_species_name, prefix, thread) {
    # output files
    fn_mntd_z <- paste0(prefix, ".mntd.z.sumtable")
    fn_mntd_p <- paste0(prefix, ".mntd.p.sumtable")
    fn_mntd_summary <- paste0(prefix, ".mntd.sumtable")

    # create doSNOW cluster
    nwcl <- makeCluster(thread)
    doSNOW::registerDoSNOW(nwcl)

    # iterate over BUSCOs
    ls_output <- foreach (busco = ls_busco, .combine='c') %dopar% {
        # initiate variables
        df_presence <- data.frame(refs=character(), reads=character(), present=numeric())

        # output variables
        mntd_z <- list(busco=busco)
        mntd_p <- list(busco=busco)
        mntd_sum <- list(busco=busco)

        # check if treefile exists
        fn_tree <- paste0(dir_busco_tree, busco, "/", busco, "_aligned.fa.treefile")
        if (!file.exists(fn_tree)){
            return(NULL)
        }

        # read tree
        tre <- ape::read.tree(fn_tree)
        ls_tips <- tre$tip.label

        # iterate over reference sequences
        for (ref in ls_refseq) {
            # extract list of taxa with specific reference
            ls_taxa <- ls_tips[stringr::str_detect(ls_tips, ref)]
            
            # remove reference sequence if variable is FALSE
            if (!is_ref_included) {
                ls_taxa <- ls_taxa[ls_taxa != ref]
            }
            
            # add taxa to the presence/absence table
            for (taxon in ls_taxa) {
                df_presence <- rbind(df_presence, c(refs=ref, reads=taxon, present=1))
            }
        }

        # transform data.frame to community format used in picante
        df_presence <- labdsv::matrify(df_presence)

        # run MNTD
        mntd_result <- picante::ses.mntd(df_presence, ape::cophenetic.phylo(tre))

        # iterate over reference sequences
        for (i in 1:nrow(mntd_result)) {
            ref_name <- rownames(mntd_result)[i]
            mntd_z[[ls_species_name[ref_name]]] <- mntd_result$mntd.obs.z[i]
            mntd_p[[ls_species_name[ref_name]]] <- mntd_result$mntd.obs.p[i]

            # update variables
            if (mntd_result$mntd.obs.z[i] > 0 && mntd_result$mntd.obs.p[i] > 0.95) {
                mntd_sum[[ls_species_name[ref_name]]] <- "S"
            } else if (mntd_result$mntd.obs.z[i] < 0 && mntd_result$mntd.obs.p[i] < 0.05) {
                mntd_sum[[ls_species_name[ref_name]]] <- "C"
            } else {
                mntd_sum[[ls_species_name[ref_name]]] <- ""
            }
        }

        return(list(mntd_z=mntd_z, mntd_p=mntd_p, mntd_sum=mntd_sum))
    }

    stopCluster(nwcl)

    # save MNTD Z score
    ls_mntd_z_out <- ls_output[names(ls_output)=="mntd_z"]
    df_mntd_z_output <- data.table::as.data.table(do.call(rbind, ls_mntd_z_out), fill=TRUE)
    data.table::fwrite(df_mntd_z_output, file=fn_mntd_z, sep="\t", quote=F, row.names=F)

    # save MNTD p-value
    ls_mntd_p_out <- ls_output[names(ls_output)=="mntd_p"]
    df_mntd_p_output <- data.table::as.data.table(do.call(rbind, ls_mntd_p_out), fill=TRUE)
    data.table::fwrite(df_mntd_p_output, file=fn_mntd_p, sep="\t", quote=F, row.names=F)

    # save significant clusters and spreads
    ls_mntd_sum_out <- ls_output[names(ls_output)=="mntd_sum"]
    df_mntd_sum_output <- data.table::as.data.table(do.call(rbind, ls_mntd_sum_out), fill=TRUE)
    data.table::fwrite(df_mntd_sum_output, file=fn_mntd_summary, sep="\t", quote=F, row.names=F)

    # plot data.frame
    f_mntd_visualization(fn_mntd_summary, prefix)
}

# function: plot MNTD results
f_mntd_visualization <- function(fn_mntd_summary, prefix) {
    # output files
    fn_mntd_cluster_tiff <- paste0(prefix, ".mntd.cluster.tiff")
    fn_mntd_spread_tiff <- paste0(prefix, ".mntd.spread.tiff")

    # open file
    df_mntd <- data.table::fread(fn_mntd_summary)

    # generate plots
    df_mntd_melt <- reshape2::melt(df_mntd, id.vars="busco")
    df_mntd_melt_cluster <- df_mntd_melt[df_mntd_melt$value == "C" | df_mntd_melt$value == "",]
    df_mntd_melt_spread <- df_mntd_melt[df_mntd_melt$value == "S" | df_mntd_melt$value == "",]

    # plot significant clusters
    tiff(file=fn_mntd_cluster_tiff, units="px", width=2880, height=1800)
    ggplot(df_mntd_melt_cluster, aes(x=variable, y=busco)) +
        geom_tile(aes(fill=value), color="white") +
        ggtitle("BUSCOs with Reference-based Clusters") + ylab("BUSCO") +
        scale_fill_manual(values=c("white","red")) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme(plot.title = element_text(hjust = 0.5, size = 50),
            plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size=30),
            axis.text.y = element_blank(),
            axis.title.y = element_text(size=40),
            legend.position = "none")
    dev.off()

    # save significnat spreads
    tiff(file=fn_mntd_spread_tiff, units="px", width=2880, height=1800)
    ggplot(df_mntd_melt_spread, aes(x=variable, y=busco)) +
        geom_tile(aes(fill=value), color="white") +
        ggtitle("BUSCOs with Reference-based Spreads") + ylab("BUSCO") +
        scale_fill_manual(values=c("white","blue")) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme(plot.title = element_text(hjust = 0.5, size = 50),
            plot.margin = margin(1.25, 1.25, 1.25, 1.25, "cm"),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size=30),
            axis.text.y = element_blank(),
            axis.title.y = element_text(size=40),
            legend.position = "none")
    dev.off()
}

# function: run entrez-direct to extract species from assembly accession number
f_extract_species_from_assembly <- function(exe_efetch, accession) {
    cmd_entrez <- paste(exe_efetch, "-db assembly -id", accession, "-format docsum")
    metadata <- system(cmd_entrez, intern=T)

    # extract the species name
    species_name <- metadata[grep("<Organism>", metadata)]
    species_name <- gsub(".*<Organism>(.+)<\\/Organism>.*", "\\1", species_name)

    return(species_name)
}

# function: run entrez-direct to extract species from SRA accession number
f_extract_species_from_reads <- function(exe_efetch, accession) {
    cmd_entrez <- paste(exe_efetch, "-db sra -id", accession, "-format runinfo")
    metadata <- system(cmd_entrez, intern=T)

    # extract index of the species name
    ls_metadata <- strsplit(metadata, split=",") 
    idx_species <- which(ls_metadata[[1]] == "ScientificName")
    species_name <- ls_metadata[[2]][idx_species]

    return(species_name)
}