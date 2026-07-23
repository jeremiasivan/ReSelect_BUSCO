# functions for codes/01_readmap/

# function: readmap fastq to reference sequence using BWA-MEM2
f_read_mapping <- function(refseq, fastq_one, fastq_two, thread, exe_bwamem2, file_sam) {
    # index reference file
    cmd_index <- paste(exe_bwamem2, "index", refseq)
    system(cmd_index)

    # readmap
    cmd_readmap <- paste(exe_bwamem2, "mem -t", thread, refseq)

    if (is.null(fastq_two)) {
        cmd_readmap <- paste(cmd_readmap, fastq_one, ">", file_sam)
    } else {
        cmd_readmap <- paste(cmd_readmap, fastq_one, fastq_two, ">", file_sam)
    }
    system(cmd_readmap)
}

# function: convert SAM to BAM
f_sam_to_bam <- function(prefix, dir_output, thread, exe_samtools) {
    # initiate variables
    fn_sam <- paste0(dir_output, "/", prefix, ".sam")
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")

    # set the number of threads
    nthread <- paste("--threads", thread)

    # run samtools
    cmd_samtools <- paste(exe_samtools, "view", nthread, "-b -u", fn_sam, "|", # convert to uncompressed (-u) BAM (-b)
                          exe_samtools, "collate", nthread, "-O -u - |",       # group reads with the same name together, output as STDOUT (-O)
                          exe_samtools, "fixmate", nthread, "-m -u - - |",     # correct flags used in the file, adding mate score tags (-m)
                          exe_samtools, "sort", nthread, "-u - |",             # sort the reads based on their positions
                          exe_samtools, "markdup -r", nthread, "-", fn_bam)    # remove duplicates based on the mate score tags
    system(cmd_samtools)

    # index BAM file
    cmd_bam_index <- paste(exe_samtools, "index", nthread, fn_bam)
    system(cmd_bam_index)
}

# function: qualiMap
f_qualimap <- function(fn_bam, dir_output, thread, fn_gff, exe_qualimap) {
    cmd_coverage <- paste(exe_qualimap, "bamqc",
                          "-bam", fn_bam,
                          "-outdir", dir_output,
                          "-nt", thread,
                          "--java-mem-size=4G")

    if (!is.null(fn_gff)) {
        cmd_coverage <- paste(cmd_coverage, "-gff", fn_gff)
    }
    
    system(cmd_coverage)
}

# function: variant calling
f_variant_calling <- function(prefix, dir_output, thread, refseq, exe_bcftools) {
    # initiate variables
    fn_bam <- paste0(dir_output, "/", prefix, ".bam")
    fn_vcf <- paste0(dir_output, "/", prefix, ".vcf.gz")
    fn_fas <- paste0(dir_output, "/", prefix, ".fa")

    nthread <- paste("--threads", thread)

    # run bcftools mpileup
    cmd_bcftools <- paste(exe_bcftools, "mpileup", nthread, "-Ou -f", refseq, fn_bam, "|",      # generate genotype likelihoods at each position with coverage
                          exe_bcftools, "call", nthread, "-Ou -mv |",                           # variant calling with default settings (-m) and output only variant sites (-v)
                          exe_bcftools, "view", nthread, "-V indels -i 'QUAL>15 & MQ>30' |",    # filter out variants with low quality score
                          exe_bcftools, "norm", nthread, "-f", refseq, "-Oz -o", fn_vcf)        # normalize variants
    system(cmd_bcftools)

    # index VCF file
    cmd_vcf_index <- paste(exe_bcftools, "index -t", nthread, fn_vcf)
    system(cmd_vcf_index)

    # generate consensus sequence
    cmd_consensus <- paste("cat", refseq, "|", exe_bcftools, "consensus", fn_vcf, ">", fn_fas)
    system(cmd_consensus)
}
