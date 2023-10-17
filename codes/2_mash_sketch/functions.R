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