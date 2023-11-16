# functions for codes/4_busco

# function: run BUSCO pipeline
f_run_busco <- function(fn_fasta, lineage, dir_output, mode, thread, exe_busco) {
    cmd_busco <- paste(exe_busco,
                       "-i", fn_fasta,
                       "-l", lineage,
                       "-m", mode,
                       "-o", dir_output,
                       "-c", thread)
    system(cmd_busco)
}