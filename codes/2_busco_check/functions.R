# functions for codes/2_busco_check

# function: check if vector is all NULL or NAs
f_all_null_or_na <- function(vector) {
    all_null_or_na <- all(is.null(vector) | is.na(vector))
    return (all_null_or_na)
}

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
        seqname=coordinates$seqname,
        source="Metaeuk",
        feature="CDS",
        start=coordinates$start+1,
        end=coordinates$stop+1,
        score=".",
        strand=coordinates$strand,
        frame=".",
        attribute=paste0("ID=",busco,"_CDS_0")
    )

    # save the new GFF file
    data.table::fwrite(df_gff, file=fn_out, sep="\t", quote=F, row.names=F, col.names=F)
}

# function: manipulate and save GFF file
f_manipulate_gff <- function(fn_input, coordinates, busco, prefix, fn_out) {
    # create a child function to check if CDS comes from multiple genes
    f_check_cds_target_id <- function(df_gff) {
        # extract attributes
        ls_attribute <- df_gff$attribute
        ls_target_id <- sapply(ls_attribute, function(x){
            stringr::str_match(x, "Target_ID=([^;]+)")[2]
        })

        if (length(unique(ls_target_id)) != 1) {
            # extract the most common target ID
            tb_target_id <- table(ls_target_id)
            target_id <- names(tb_target_id[tb_target_id==max(tb_target_id)])[1]

            # extract the GFF entries
            df_gff_filter <- df_gff[grepl(target_id, df_gff$attribute),]
            return(df_gff_filter)
        }

        return(df_gff)
    }

    # read GFF file
    df_gff <- data.table::fread(fn_input, header=FALSE, select=1:9)
    data.table::setnames(df_gff, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

    # convert attribute to match with gffread naming convention
    df_gff$attribute <- gsub("TCS_ID", "ID", df_gff$attribute)

    # extract relevant GFF entry
    df_gff <- df_gff[df_gff$seqname==coordinates$seqname & df_gff$strand==coordinates$strand]
    if (nrow(df_gff) == 0) {
        return(list(errmsg=paste0("Error: ", busco, " GFF entry is not found for ", prefix, ". Skipped.")))
    }

    # extract the respective gene index
    gene_idx <- which(df_gff$feature == "gene" & df_gff$start == coordinates$start & df_gff$end == coordinates$stop)
    if (length(gene_idx) != 1) {
        # extract all CDS based on the coordinates
        df_gff_cds <- df_gff[df_gff$feature=="CDS" & df_gff$start>=coordinates$start & df_gff$end<=coordinates$stop]
        if (nrow(df_gff_cds) == 0) {
            return(list(errmsg=paste("Error:", busco, "for", prefix, "has invalid CDS coordinates. Skipped.")))
        }

        # remove entry with length == 1
        df_gff_cds <- df_gff_cds[abs(df_gff_cds$start-df_gff_cds$end)!=1,]

        # remove multiple genes
        df_gff_cds <- f_check_cds_target_id(df_gff_cds)

        # save the new GFF file
        data.table::fwrite(df_gff_cds, file=fn_out, sep="\t", quote=F, row.names=F, col.names=F)
        return(list(warnmsg=paste("Warn:", busco, "GFF file for", prefix, "is based on CDS coordinates.")))
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
    
    # remove entry with length == 1
    df_gff_subset <- df_gff_subset[abs(df_gff_subset$start-df_gff_subset$end)!=1,]

    # remove multiple genes
    df_gff_subset <- f_check_cds_target_id(df_gff_subset)

    # save the new GFF file
    df_gff_subset_cds <- df_gff_subset[df_gff_subset$feature=="CDS",]
    data.table::fwrite(df_gff_subset_cds, file=fn_out, sep="\t", quote=F, row.names=F, col.names=F)
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

# function: extract all BUSCO alignments from GFF
f_extract_fasta_from_gff <- function(fn_input, fn_gff, fn_cds_out, fn_concat_out, exe_gffread){
    # run gffread
    cmd_gffread <- paste(exe_gffread, "-g", fn_input, "-y", fn_cds_out, fn_gff)
    system(cmd_gffread)

    # open the CDS sequences
    all_cds <- Biostrings::readAAStringSet(fn_cds_out)
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
    concat_cds <- Biostrings::AAStringSet(seq)
    names(concat_cds) <- header
    Biostrings::writeXStringSet(concat_cds, filepath=fn_concat_out)
}

# function: compare two FASTA files
f_compare_fasta <- function(fn_fasta_one, fn_fasta_two) {
    # open the two FASTA files
    fasta_one <- Biostrings::readAAStringSet(fn_fasta_one)
    fasta_two <- Biostrings::readAAStringSet(fn_fasta_two)

    # check if there is only one sequence
    header_one <- names(fasta_one)
    header_two <- names(fasta_two)

    if (length(header_one) != 1 || length(header_two) != 1) {
        return(list(errmsg="Non-identical header"))
    }

    # compare if the two are equal
    is_identical <- identical(fasta_one[[header_one]], fasta_two[[header_two]])
    if (is_identical) {
        return(list(is_identical=TRUE))
    }

    # if not equal, run pairwise distance calculation
    fasta_msa <- Biostrings::pairwiseAlignment(fasta_one, fasta_two, substitutionMatrix = "BLOSUM100")
    return(list(is_identical=FALSE, score=Biostrings::pid(fasta_msa)))
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

# function: run TreeShrink
f_treeshrink <- function(fn_input, prefix, dir_output, exe_treeshrink) {
    cmd_treeshrink <- paste("python", exe_treeshrink,
                            "-t", fn_input,
                            "-O", prefix,
                            "-o", dir_output)
    system(cmd_treeshrink)
}

# function: change tips from accesion ID to species name
f_tips_to_species <- function(fn_treefile, fn_treefile_species, ls_species_name) {
    # read input file
    tre <- readLines(fn_treefile)

    # iterate over tips
    for (id in names(ls_species_name)) {
        tre <- gsub(id, ls_species_name[id], tre)
    }

    # save the file
    writeLines(tre, con=fn_treefile_species)
}

# function: change tips from species name to accesion ID for reference sequences
f_ref_tips_to_id <- function(fn_treefile, df_refs) {
    # read input file
    tre <- readLines(fn_treefile)

    # iterate over tips
    for (i in 1:nrow(df_refs)) {
        tre <- gsub(df_refs$species[i], df_refs$id[i], tre)
    }

    return(tre)
}

# function: extract R2 and p-value
f_extract_summary_lm <- function(lm_result) {
    # extract R2
    rsquared <- tryCatch({
        round(summary(lm_result)$r.squared, 3)
    }, error = function(err) {
        return(NULL)
    })

    # extract p-value
    pvalue <- tryCatch({
        round(summary(lm_result)$coefficients[2,4], 3)
    }, error = function(err) {
        return(NULL)
    })

    # extract slope
    slope <- round(coef(lm_result)[2],3)

    return(list(rsquared=rsquared, pvalue=pvalue, slope=slope))
}

# function: run Spearman correlatio test
f_spearman_test <- function(x, y) {
    corr <- cor.test(x, y, method='spearman')

    return(list(rho=round(corr$estimate, 3), pvalue=round(corr$p.value, 3)))
}

# function: check the closest reference given reads
f_check_closest_ref <- function(dist_matrix, species_read, ls_species_ref) {
    # initiate variables
    closest_ref <- ls_species_ref[1]
    closest_ref_dist <- dist_matrix[rownames(dist_matrix) == species_read, colnames(dist_matrix) == ls_species_ref[1]]

    # return the closest reference if there is only one reference
    len_ls_species_ref <- length(ls_species_ref)
    if (len_ls_species_ref == 1) {
        return (list(ref=closest_ref, dist=round(closest_ref_dist, 3)))
    }

    # iterate over reference sequence
    for (ref in ls_species_ref[2:len_ls_species_ref]) {
        # check the phylogenetic distance
        ref_dist <- dist_matrix[rownames(dist_matrix) == species_read, colnames(dist_matrix) == ref]

        # update the closest reference is distance is smaller
        if (ref_dist < closest_ref_dist) {
            closest_ref <- ref
            closest_ref_dist <- ref_dist
        }
    }

    # return the closest reference
    return (list(ref=closest_ref, dist=round(closest_ref_dist, 3)))
}

# function: collapse branch with low bootstrap value [NOT USED]
f_collapse_branch <- function(fn_tree, bootstrap, fn_out, exe_nwed) {
    newick_cmd <- paste0(exe_nwed, " ", fn_tree, " 'i & b<", bootstrap, "' o > ", fn_out)
    system(newick_cmd)
}

# function: calculate normalised RF distance between two trees
f_calculate_nRF <- function(fn_tree_one, fn_tree_two) {
    # open the two treefiles
    tree_one <- ape::read.tree(fn_tree_one)
    tree_two <- ape::read.tree(fn_tree_two)

    # calculate nRF
    nrf_dist <- phangorn::RF.dist(tree_one, tree_two, normalize=TRUE)

    return(round(nrf_dist,3))
}

# function: extract distance value based on the number of highly-supported branches
f_calculate_treedist <- function(fn_gene_tree, fn_refs_tree, min_bootstrap) {
    # child function to get the tips under a given node (source: chatGPT)
    f_get_tips <- function(tree, node) {
        descendants <- phangorn::Descendants(tree, node, type = "tips")
        return(tree$tip.label[unlist(descendants)])
    }

    # child function to check if tips are monophyletic
    f_is_monophyletic <- function(tips, min_bootstrap, tree) {
        # initial check
        is_monophyletic <- ape::is.monophyletic(tree, tips)
        if (is_monophyletic) {
            return(TRUE)
        }

        # extract bootstrap values from the tree
        df_node_bs <- data.frame(node = (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), 
                                 bootstrap = as.numeric(tree$node.label))
        df_node_bs <- na.omit(df_node_bs)

        # extract the MRCA node
        node_mrca <- ape::getMRCA(tree, tips)

        # iterate over tips
        ls_nodepath <- c()
        for (tip in tips) {
            tip_node <- which(tree$tip.label==tip)
            ls_nodepath <- c(ls_nodepath, ape::nodepath(tree, tip_node, node_mrca))
        }

        # extract unique nodes
        ls_nodepath <- unique(ls_nodepath)

        # extract bootstrap values
        ls_bootstrap_values <- df_node_bs$bootstrap[df_node_bs$node%in%c(node_mrca,ls_nodepath)]
        if (all(ls_bootstrap_values < min_bootstrap)) {
            return(TRUE)
        }

        return(FALSE)
    }

    # open the trees
    gene_tree <- ape::read.tree(fn_gene_tree)
    refs_tree <- ape::read.tree(fn_refs_tree)

    # extract bootstrap values and their respective nodes from gene tree (source: ChatGPT)
    df_node_bootstrap <- data.frame(node = (length(gene_tree$tip.label) + 1):(length(gene_tree$tip.label) + gene_tree$Nnode), 
                                    bootstrap = as.numeric(gene_tree$node.label))

    # remove rows with NA bootstrap values (if any)
    df_node_bootstrap <- na.omit(df_node_bootstrap)

    # nRF value
    nRF <- 0
    n_infbranch <- 0
    nRF_increment <- 1 / (length(gene_tree$tip.label)-3)

    # iterate over branching events with high bootstrap value
    for (i in 1:nrow(df_node_bootstrap)) {
        # check if bootstrap value is low
        if (df_node_bootstrap$bootstrap[i] < min_bootstrap) {
            next
        }

        # get the child nodes (source: ChatGPT)
        child_nodes <- gene_tree$edge[gene_tree$edge[,1] == df_node_bootstrap$node[i], 2]

        # get the tips under each child node (source: ChatGPT)
        group1_tips <- f_get_tips(gene_tree, child_nodes[1])
        group2_tips <- f_get_tips(gene_tree, child_nodes[2])

        # check if both groups are monophyletic
        group1_monophyletic <- f_is_monophyletic(group1_tips, min_bootstrap, refs_tree)
        group2_monophyletic <- f_is_monophyletic(group2_tips, min_bootstrap, refs_tree)
        if (!group1_monophyletic || !group2_monophyletic) {
            nRF <- nRF + nRF_increment
        }

        n_infbranch <- n_infbranch + 1
    }

    # return value for nRF
    return(list(dist=round(nRF,3), n_infbranch=n_infbranch))
}

# function: extract average bootstrap value
f_calculate_mean_bs <- function(fn_tree) {
    # open treefile
    tree <- ape::read.tree(fn_tree)

    # calculate mean bootstrap values
    mean_bs <- mean(as.numeric(tree$node.label[tree$node.label!=""]))

    return(mean_bs)
}

# function: run ASTRAL-III 
f_astral <- function(fn_input, fn_output, fn_log, exe_astral) {
    cmd_astral <- paste("java -jar", exe_astral,
                    "-i", fn_input,
                    "-o", fn_output,
                    "-t 2 2>", fn_log)
    system(cmd_astral)
}