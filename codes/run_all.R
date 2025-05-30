#################################
codedir <- "~/ReSelect_BUSCO/codes"
outprefix <- ""
outdir <- "~/busco"
thread <- 5
redo <- TRUE

# data download
file_refseq <- paste0(codedir, "/../data/eucs_refseq.txt")
file_shortreads <- paste0(codedir, "/../data/eucs_shortreads.txt")
file_adapters <- paste0(codedir, "/../data/eucs_adapter.txt")

exe_datasets <- "datasets"
bin_sratoolkit <- "sratoolkit/bin/"

# BUSCO check
exe_adapterremoval <- "AdapterRemoval"
exe_bwamem2 <- "bwa-mem2"
exe_samtools <- "samtools"
exe_bcftools <- "bcftools"
exe_qualimap <- "qualimap"

exe_gffread <- "gffread"
exe_gff2bed <- "bedops/bin/gff2bed"

exe_busco <- "busco"
exe_iqtree2 <- "iqtree2"
exe_mafft <- "mafft"
exe_treeshrink <- "treeshrink"

exe_astral <- "astral.jar"

busco_lineage <- "eudicots_odb10"
busco_mode <- "genome"
type <- "coordinate"

min_busco_depth <- 10

outgroup <- c("Angophora","Corymbia")

# BUSCO tree
file_buscotree <- paste0(codedir, "/../data/eucs_buscotree.txt")
file_genome_treefile <- paste0(codedir, "/../data/eucs.treefile")
busco_tree_mode <- "random"
busco_tree_random_ref <- TRUE
#################################

# extract current time
current_time <- format(Sys.time(), "%y-%m-%d %H:%M")
current_time <- gsub("[: -]", "" , current_time, perl=TRUE)

# set up outdir
outdir_prefix <- paste0(outdir, "/", outprefix, "/")
if (!dir.exists(outdir_prefix)) {
    dir.create(outdir_prefix, recursive=T)
}

# data download
rmarkdown::render(input=paste0(codedir,"/1_data_preparation/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, outprefix, ".download.html"),
                  params=list(prefix=outprefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              file_refseq=file_refseq, file_shortreads=file_shortreads, file_adapters=file_adapters,
                              exe_datasets=exe_datasets, bin_sratoolkit=bin_sratoolkit,
                              exe_adapterremoval=exe_adapterremoval, exe_bwamem2=exe_bwamem2, exe_samtools=exe_samtools, exe_bcftools=exe_bcftools, exe_qualimap=exe_qualimap,
                              current_time=current_time),
                  quiet=TRUE)

# BUSCO check
rmarkdown::render(input=paste0(codedir,"/2_busco_check/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, outprefix, ".check.html"),
                  params=list(prefix=outprefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              file_refseq=file_refseq, file_shortreads=file_shortreads, file_genome_treefile=file_genome_treefile, file_buscotree=file_buscotree,
                              exe_busco=exe_busco, exe_gffread=exe_gffread, exe_samtools=exe_samtools, exe_gff2bed=exe_gff2bed, exe_qualimap=exe_qualimap,
                              exe_iqtree2=exe_iqtree2, exe_mafft=exe_mafft, exe_treeshrink=exe_treeshrink, exe_astral=exe_astral,
                              busco_lineage=busco_lineage, busco_mode=busco_mode, type=type, outgroup=outgroup,
                              min_busco_depth=min_busco_depth, busco_tree_mode=busco_tree_mode, busco_tree_random_ref=busco_tree_random_ref,
                              current_time=current_time),
                  quiet=TRUE)

#################################