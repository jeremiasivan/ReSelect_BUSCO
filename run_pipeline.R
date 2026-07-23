#!/usr/bin/env Rscript

# ============================================================
#  BusIER — 2_phylorbt
#
#  Usage: Rscript run_pipeline.R --config config.yaml
#         Rscript run_pipeline.R --config config.yaml --redo
# ============================================================

# --- Load libraries and function ----------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(yaml))

# create a function to retrieve parameter value with default
f_get_param <- function(value, default) {
  if (is.null(value) || identical(value, "")) {
    default
  } else {
    value
  }
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL,
              help="Path to YAML config file [required]", metavar="FILE"),
  make_option(c("-r", "--redo"), action="store_true", default=FALSE,
              help="Re-run all analyses and override previous results")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if config file is invalid
if (is.null(opt$config)) {
  stop("-c/--config is required.\n
        Usage: Rscript run_pipeline.R --config config.yaml")
}

if (!file.exists(opt$config)) {
  stop(paste("Config file not found:", opt$config))
}

# --- Load and validate config --------------------------------
cfg <- yaml::read_yaml(opt$config)

# set required parameters
required_fields <- c("codedir", "outdir", "file_refseq", "file_shortreads",
                      "file_genome_treefile", "dir_genes_alignment",
                      "exe_iqtree2", "exe_astral")
missing <- setdiff(required_fields, names(cfg))
if (length(missing) > 0) {
  stop(paste("Missing required config fields:", paste(missing, collapse=", ")))
}

# check if input metadata files are invalid
if (is.null(cfg$file_refseq) || cfg$file_refseq == "") {
  stop("file_refseq must be set in the config file.")
}
if (!file.exists(path.expand(cfg$file_refseq))) {
  stop(paste("file_refseq file not found:", cfg$file_refseq))
}

if (is.null(cfg$file_shortreads) || cfg$file_shortreads == "") {
  stop("file_shortreads must be set in the config file.")
}
if (!file.exists(path.expand(cfg$file_shortreads))) {
  stop(paste("file_shortreads file not found:", cfg$file_shortreads))
}

if (is.null(cfg$file_genome_treefile) || cfg$file_genome_treefile == "") {
  stop("file_genome_treefile must be set in the config file.")
}
if (!file.exists(path.expand(cfg$file_genome_treefile))) {
  stop(paste("file_genome_treefile file not found:", cfg$file_genome_treefile))
}

# dir_genes_alignment is produced by 1_data_preparation and must already exist
if (is.null(cfg$dir_genes_alignment) || cfg$dir_genes_alignment == "") {
  stop("dir_genes_alignment must be set in the config file.")
}
if (!dir.exists(path.expand(cfg$dir_genes_alignment))) {
  stop(paste("dir_genes_alignment not found -- run 1_data_preparation first:", cfg$dir_genes_alignment))
}

# exe_iqtree2 / exe_astral have no sensible default -- must be set explicitly
if (is.null(cfg$exe_iqtree2) || cfg$exe_iqtree2 == "") {
  stop("exe_iqtree2 must be set in the config file.")
}
if (is.null(cfg$exe_astral) || cfg$exe_astral == "") {
  stop("exe_astral must be set in the config file.")
}

# check if prefix is set, otherwise derive it from the genome treefile filename
if (is.null(cfg$prefix) || cfg$prefix == "") {
  cfg$prefix <- tools::file_path_sans_ext(basename(cfg$file_genome_treefile))
  message("Note: prefix not set, using \"", cfg$prefix, "\" derived from file_genome_treefile.")
  message("      If you want this stage written into the same folder as 1_data_preparation,")
  message("      set prefix explicitly to match that run's prefix.")
}

# --- Apply CLI overrides -------------------------------------
if (opt$redo) {
  message("Note: --redo flag set via CLI, overriding config.")
  cfg$redo <- TRUE
}

# --- Map config to rmarkdown params --------------------------
render_params <- list(
  codedir              = cfg$codedir,
  prefix               = cfg$prefix,
  outdir               = cfg$outdir,
  thread               = as.integer(f_get_param(cfg$thread, 1)),
  redo                 = as.logical(f_get_param(cfg$redo, FALSE)),

  file_refseq          = cfg$file_refseq,
  file_shortreads      = cfg$file_shortreads,
  file_genome_treefile = cfg$file_genome_treefile,

  dir_genes_alignment  = cfg$dir_genes_alignment,

  exe_iqtree2          = cfg$exe_iqtree2,
  exe_astral           = cfg$exe_astral,

  is_astral_constrained = as.logical(f_get_param(cfg$is_astral_constrained, FALSE)),

  outgroup             = as.character(unlist(f_get_param(cfg$outgroup, list("")))),
  focal_species        = f_get_param(cfg$focal_species, "")
)

# --- Run PhyloRBT --------------------------------------------
rmd_path <- file.path(path.expand(render_params$codedir), "2_phylorbt", "1_main.Rmd")
if (!file.exists(rmd_path)) {
  stop(paste("1_main.Rmd not found:", rmd_path))
}

output_dir <- file.path(path.expand(render_params$outdir), render_params$prefix)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive=TRUE)
}

message("Starting PhyloRBT pipeline...")
message("  Config:          ", opt$config)
message("  Prefix:          ", render_params$prefix)
message("  Output:          ", render_params$outdir)
message("  Gene alignments: ", render_params$dir_genes_alignment)
message("  Genome treefile: ", render_params$file_genome_treefile)
message("  Threads:         ", render_params$thread)

# render the Rmarkdown file
rmarkdown::render(
  input       = rmd_path,
  params      = render_params,
  output_file = paste0(render_params$prefix, "_phylorbt_report.html"),
  output_dir  = output_dir,
  quiet       = FALSE
)

message("Done. Report: ",
        file.path(output_dir, paste0(render_params$prefix, "_phylorbt_report.html")))
