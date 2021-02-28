PARAM_CAGE_XIAOYUN_ADAPT <- '-a "ACACATCTCCGAGCCCACGAGAC;min_overlap=10" -A "CCCTATA;min_overlap=7" '
                                 # ACATCTCCGAGCCCACGAGAC is the one Xiaoyun gave me through slack (I had to reverse complement it to get this sequence)
# this param is adapted from juheon's /scratch/jmaeng/TE_TF_cancer/TE_gene_chimericreads/script/pipe_align_nanoCAGE.sh
# cutadapt suggested that the -a adaptor might not be complete. additional "TAT" might need to be added to the  beginning.

cutadapt.fanc <- function(fastqs, out.files = NULL, main.params = PARAM_CAGE_XIAOYUN_ADAPT, thread,
                          min.read.length, run = T, stdout.file = NULL, other.params = "",
                          cutadapt.path = "/bar/cfan/anaconda2/envs/jupyter/bin/cutadapt") {
  if (is.null(out.files)) {
    out.files <- fastqs %>% sub("READ1", "trimmed_READ1", .) %>% sub("_R1", "_trimmed_R1", .) %>% sub("_1.fastq", "_trimmed_1.fastq", .) %>%
      sub("READ2", "trimmed_READ2", .) %>% sub("_R2", "_trimmed_R2", .) %>% sub("_2.fastq", "_trimmed_2.fastq", .)
  }
  cmd <- paste0(cutadapt.path, " -j ", thread," ", main.params, " -o ", out.files[1], " -p ", out.files[2],
                " -m ", min.read.length,
                " ",paste0(other.params, collapse = " "), " ", paste0(fastqs, collapse =" "))
  utilsFanc::cmd.exec.fanc(cmd = cmd, stdout.file = stdout.file, intern = F, run = run)
  if (sum(!file.exists(out.files)) > 0) {
    stop(paste0("at least one of the files of ", paste0(out.files, collapse = ";"),
                " was not successfully generated"))
  } else
    return(out.files)
}
