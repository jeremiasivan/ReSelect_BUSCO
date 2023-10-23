#################################
codedir <- "~/BusIER/codes"
prefix <- ""
outdir <- "~/busco"
thread <- 5
redo <- TRUE

# data download
file_refseq <- ""
file_shortreads <- ""

exe_datasets <- ""
dir_sratoolkit <- ""

# mash sketch
exe_mash <- ""
exe_efetch <- ""

# read mapping
exe_bwamem2 <- ""
#################################

# set up outdir
outdir_prefix <- paste0(outdir, "/", prefix, "/")
if (!dir.exists(outdir_prefix)) {
    dir.create(outdir_prefix, recursive=T)
}

# data download
rmarkdown::render(input=paste0(codedir,"/1_data_download/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, prefix, ".download.html"),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              file_refseq=file_refseq, file_shortreads=file_shortreads,
                              exe_datasets=exe_datasets, dir_sratoolkit=dir_sratoolkit),
                  quiet=TRUE)

# mash sketch
rmarkdown::render(input=paste0(codedir,"/2_mash_sketch/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, prefix, ".mash.html"),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              exe_mash=exe_mash, exe_efetch=exe_efetch),
                  quiet=TRUE)

# read mapping
rmarkdown::render(input=paste0(codedir,"/3_read_mapping/1_main.Rmd"),
                  output_file=paste0(outdir_prefix, prefix, ".readmap.html"),
                  params=list(prefix=prefix, codedir=codedir, outdir=outdir, thread=thread, redo=redo,
                              exe_bwamem2=exe_bwamem2),
                  quiet=TRUE)

#################################