# functions for codes/2_mash_sketch

# function: create Mash sketch for all reference alignments
f_mash_sketch <- function(exe_mash, str_refseq, thread, prefix) {
    cmd <- paste(exe_mash, "sketch",
                "-o", prefix,
                "-p", thread,
                str_refseq)
    system(cmd)
}

# function: run Mash screen on short reads
f_mash_screen <- function(exe_mash, file_msh, file_input, file_output) {
    cmd_mash <- paste(exe_mash, "screen -w", file_msh, file_input, ">", file_output)
    system(cmd_mash)
}

# function: run entrez-direct to extract species from SRA accession number
f_extract_species_from_reads <- function(exe_esearch, exe_efetch, accession) {
    cmd_entrez <- paste(exe_esearch, "-db sra -query", accession,
                        "|", exe_efetch, "-format runinfo")
    metadata <- system(cmd_entrez, intern=T)

    # extract index of the species name
    ls_metadata <- strsplit(metadata, split=",") 
    idx_species <- which(ls_metadata[[1]] == "ScientificName")
    species_name <- ls_metadata[[2]][idx_species]

    return(species_name)
}

# function: run entrez-direct to extract species from assembly accession number
f_extract_species_from_assembly <- function(exe_esearch, exe_efetch, exe_xtract, accession) {
    cmd_entrez <- paste(exe_esearch, "-db assembly -query", accession,
                        "|", exe_efetch, "-format docsum",
                        "|", exe_xtract, "-pattern DocumentSummary -element Organism")
    species_name <- system(cmd_entrez, intern=T)

    return(species_name)
}