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
                         "-bb 1000",
                         "-T 1 --quiet -redo")
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

# function: run Spearman correlation test
f_spearman_test <- function(x, y) {
    corr <- cor.test(x, y, method='spearman')

    return(list(rho=round(corr$estimate, 3), pvalue=round(corr$p.value, 3)))
}

# function: extract mapped reads from an alignment
f_extract_fasta <- function(fn_input, ls_header, fn_output) {
    # open alignment
    seq <- Biostrings::readAAStringSet(fn_input)

    # extract reference sequences and mapped reads
    seq_subset <- seq[names(seq)%in%ls_header]

    # save the file
    Biostrings::writeXStringSet(seq_subset, filepath=fn_output)
}

# function: generate window trees
f_multiple_tree <- function(dir_aln, prefix, thread, dir_iqtree2) {
    iqtree_cmd <- paste(dir_iqtree2,
                        "-S", dir_aln,
                        "-pre", prefix,
                        "-T", thread,
                        "--quiet -redo")
    system(iqtree_cmd)
}

# function: run ASTRAL-III 
f_astral <- function(fn_input, fn_output, fn_log, exe_astral) {
    cmd_astral <- paste("java -jar", exe_astral,
                    "-i", fn_input,
                    "-o", fn_output,
                    "-t 2 2>", fn_log)
    system(cmd_astral)
}

# function: run ASTRAL-III (constrained)
f_astral_constrained <- function(fn_input, fn_output, fn_log, fn_sptree, exe_astral) {
    cmd_astral <- paste("java -jar", exe_astral,
                    "-i", fn_input,
                    "-j", fn_sptree,
                    "-o", fn_output,
                    "-t 2 2>", fn_log)
    system(cmd_astral)
}

# function: calculate sCF and gCF
f_calculate_cf <- function(fn_all_trees, fn_sp_tree, dir_fasta, dir_output, thread, exe_iqtree2) {
    # calculate gCF
    cmd_gcf <- paste(exe_iqtree2,
                     "-t", fn_sp_tree,
                     "--gcf", fn_all_trees,
                     "-T", thread,
                     "--prefix", paste0(dir_output, "/gcf"))
    system(cmd_gcf)

    # calculate sCF
    cmd_scf <- paste(exe_iqtree2,
                     "-te", fn_sp_tree,
                     "-p", dir_fasta,
                     "--scfl 100",
                     "-T", thread,
                     "--prefix", paste0(dir_output, "/scf"))
    system(cmd_scf)
}

# function: extract support values
f_extract_branch_supports <- function(fn_input, ls_header) {
    # open trees
    tre <- ape::read.tree(fn_input)
    tre_mrca <- ape::getMRCA(tre_astral, ls_header)
    tre_supp <- tre$node.label[tre_mrca - length(tre$tip.label)]
    tre_supp <- gsub("'\\[", "", tre_supp)
    tre_supp <- gsub("\\]'", "", tre_supp)
    tre_supp <- strsplit(tre_supp, split=";")[[1]]

    # extract support values
    supp_pairs <- strsplit(tre_supp, "=")
    supp_pairs <- setNames(
        lapply(supp_pairs, function(x) as.numeric(x[2])),
        sapply(supp_pairs, function(x) x[1])
    )

    return(supp_pairs)
}