#!/usr/bin/env Rscript

# ============================================================
#  PhyloRBT — 1_data_preparation
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
required_fields <- c("codedir", "outdir", "file_refseq", "file_shortreads", "busco_lineage")
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

# BUSCO lineage must be set explicitly
if (is.null(cfg$busco_lineage) || cfg$busco_lineage == "") {
  stop("busco_lineage must be set in the config file.")
}

# check if prefix is set, otherwise derive it from the reference metadata filename
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
full_params <- list(
  codedir          = cfg$codedir,
  prefix           = cfg$prefix,
  outdir           = cfg$outdir,
  thread           = cfg$thread,
  redo             = as.logical(f_get_param(cfg$redo, FALSE)),

  file_refseq      = cfg$file_refseq,
  file_shortreads  = cfg$file_shortreads,

  exe_bwamem2      = f_get_param(cfg$exe_bwamem2, "bwa-mem2"),
  exe_samtools     = f_get_param(cfg$exe_samtools, "samtools"),
  exe_bcftools     = f_get_param(cfg$exe_bcftools, "bcftools"),
  exe_qualimap     = f_get_param(cfg$exe_qualimap, "qualimap"),

  exe_busco        = f_get_param(cfg$exe_busco, "busco"),
  exe_gff2bed      = f_get_param(cfg$exe_gff2bed, "gff2bed"),
  exe_mafft        = f_get_param(cfg$exe_mafft, "mafft"),
  exe_trimal       = f_get_param(cfg$exe_trimal, "trimal"),
  exe_iqtree2      = f_get_param(cfg$exe_iqtree2, "iqtree2"),
  exe_treeshrink   = f_get_param(cfg$exe_treeshrink, "run_treeshrink.py"),

  thread_bwamem2   = as.integer(f_get_param(cfg$thread_bwamem2, cfg$thread)),
  thread_samtools  = as.integer(f_get_param(cfg$thread_samtools, cfg$thread)),
  thread_qualimap  = as.integer(f_get_param(cfg$thread_qualimap, cfg$thread)),

  busco_lineage    = cfg$busco_lineage,
  busco_mode       = f_get_param(cfg$busco_mode, "genome"),
  thread_busco     = as.integer(f_get_param(cfg$thread_busco, cfg$thread)),
  min_busco_depth  = as.integer(f_get_param(cfg$min_busco_depth, 10))
)

# extract parameters for each analysis
readmap_param_names <- c("codedir", "prefix", "outdir", "thread", "redo",
                         "file_refseq", "file_shortreads",
                         "exe_bwamem2", "exe_samtools", "exe_bcftools", "exe_qualimap",
                         "thread_bwamem2", "thread_samtools", "thread_qualimap")

busco_param_names   <- c("codedir", "prefix", "outdir", "thread", "redo",
                         "file_refseq", "file_shortreads",
                         "exe_busco", "exe_gff2bed", "exe_samtools", "exe_qualimap", "exe_mafft", "exe_trimal", "exe_iqtree2", "exe_treeshrink",
                         "busco_lineage", "busco_mode", "thread_busco", "min_busco_depth")

params_readmap <- full_params[readmap_param_names]
params_busco   <- full_params[busco_param_names]

# --- Run PhyloRBT data preparation ----------------------------
readmap_rmd_path <- file.path(path.expand(full_params$codedir), "data_preparation", "1_readmap", "1_main.Rmd")
if (!file.exists(readmap_rmd_path)) {
  stop(paste("1_main.Rmd not found:", readmap_rmd_path))
}

busco_rmd_path <- file.path(path.expand(full_params$codedir), "data_preparation", "2_extract_busco_loci", "1_main.Rmd")
if (!file.exists(busco_rmd_path)) {
  stop(paste("1_main.Rmd not found:", busco_rmd_path))
}

message("Starting PhyloRBT data preparation pipeline...")
message("  Config:         ", opt$config)
message("  Prefix:         ", full_params$prefix)
message("  Output:         ", full_params$outdir)
message("  Reference meta: ", full_params$file_refseq)
message("  Short-read meta:", full_params$file_shortreads)
message("  Threads:        ", full_params$thread)

# --- Step 1: Read mapping -------------------------------------
message("\n[1/2] Running read mapping (1_readmap)...")
rmarkdown::render(
  input       = readmap_rmd_path,
  params      = params_readmap,
  output_file = paste0(full_params$prefix, "_readmap_report.html"),
  output_dir  = file.path(path.expand(full_params$outdir), full_params$prefix),
  quiet       = FALSE
)

# --- Step 2: BUSCO loci extraction ------------------------------
message("\n[2/2] Running BUSCO loci extraction (2_extract_busco_loci)...")
rmarkdown::render(
  input       = busco_rmd_path,
  params      = params_busco,
  output_file = paste0(full_params$prefix, "_buscoextract_report.html"),
  output_dir  = file.path(path.expand(full_params$outdir), full_params$prefix),
  quiet       = FALSE
)

message("\nDone. Reports:")
message("  ", file.path(file.path(path.expand(full_params$outdir), full_params$prefix), paste0(full_params$prefix, "_readmap_report.html")))
message("  ", file.path(file.path(path.expand(full_params$outdir), full_params$prefix), paste0(full_params$prefix, "_buscoextract_report.html")))
