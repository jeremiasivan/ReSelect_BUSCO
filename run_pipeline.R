#!/usr/bin/env Rscript

# ============================================================
#  PhyloRBT
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
required_fields <- c("codedir", "outdir", "dir_genes_alignment", "file_species_treefile")
missing <- setdiff(required_fields, names(cfg))
if (length(missing) > 0) {
  stop(paste("Missing required config fields:", paste(missing, collapse=", ")))
}

# check if input files are invalid
if (is.null(cfg$file_species_treefile) || cfg$file_species_treefile == "") {
  stop("file_species_treefile must be set in the config file.")
}

if (!file.exists(path.expand(cfg$file_species_treefile))) {
  stop(paste("file_species_treefile file not found:", cfg$file_species_treefile))
}

if (is.null(cfg$dir_genes_alignment) || cfg$dir_genes_alignment == "") {
  stop("dir_genes_alignment must be set in the config file.")
}

if (!dir.exists(path.expand(cfg$dir_genes_alignment))) {
  stop(paste("dir_genes_alignment not found:", cfg$dir_genes_alignment))
}

# check if prefix is set
if (is.null(cfg$prefix) || cfg$prefix == "") {
  cfg$prefix <- "PhyloRBT_output"
}

# update the number of thread
cfg$thread <- as.integer(f_get_param(cfg$thread, 1))

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
  thread               = cfg$thread,
  redo                 = as.logical(f_get_param(cfg$redo, FALSE)),

  dir_genes_alignment   = cfg$dir_genes_alignment,
  file_species_treefile = cfg$file_species_treefile,

  exe_iqtree2      = f_get_param(cfg$exe_iqtree2, "iqtree2"),
  exe_astral       = f_get_param(cfg$astral, "astral.jar"),

  is_astral_constrained = as.logical(f_get_param(cfg$is_astral_constrained, FALSE)),

  outgroup             = f_get_param(cfg$outgroup, ""),
  focal_species        = f_get_param(cfg$focal_species, "")
)

# --- Run PhyloRBT --------------------------------------------
rmd_path <- file.path(path.expand(render_params$codedir), "codes", "1_main.Rmd")
if (!file.exists(rmd_path)) {
  stop(paste("1_main.Rmd not found:", rmd_path))
}

message("Starting PhyloRBT pipeline...")
message("  Config:          ", opt$config)
message("  Prefix:          ", render_params$prefix)
message("  Output:          ", render_params$outdir)
message("  Gene alignments: ", render_params$dir_genes_alignment)
message("  Genome treefile: ", render_params$file_species_treefile)
message("  Threads:         ", render_params$thread)

# render the Rmarkdown file
rmarkdown::render(
  input       = rmd_path,
  params      = render_params,
  output_file = paste0(render_params$prefix, "_phylorbt_report.html"),
  output_dir  = file.path(path.expand(render_params$outdir), render_params$prefix),
  quiet       = FALSE
)

message("Done. Report: ",
        file.path(path.expand(render_params$outdir), render_params$prefix, paste0(render_params$prefix, "_phylorbt_report.html")))
