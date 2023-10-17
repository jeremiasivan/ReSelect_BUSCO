# functions for codes/2_mash_sketch

# function: create Mash Sketch for all reference alignments
f_mash_sketch <- function(exe_mash, str_refseq, thread, prefix) {
    cmd <- paste(exe_mash, "sketch",
                "-o", prefix,
                "-p", thread,
                str_refseq)
    system(cmd)
}
