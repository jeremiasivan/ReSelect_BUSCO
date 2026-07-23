# functions for codes/2_busco_check

# function: check if vector is all NULL or NAs
f_all_null_or_na <- function(vector) {
    all_null_or_na <- all(is.null(vector) | is.na(vector))
    return (all_null_or_na)
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
                     "-redo",
                     "--prefix", paste0(dir_output, "/gcf"))
    system(cmd_gcf)

    # calculate sCF
    cmd_scf <- paste(exe_iqtree2,
                     "-te", fn_sp_tree,
                     "-p", dir_fasta,
                     "--scfl 100",
                     "-T", thread,
                     "-redo",
                     "--prefix", paste0(dir_output, "/scf"))
    system(cmd_scf)
}
