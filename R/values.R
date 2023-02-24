BAMCOVERAGE <- "/bar/cfan/anaconda2/envs/jupyter/bin/bamCoverage"
host.name <- Sys.info()["nodename"][1]
if (grepl("ris", host.name)) {
    # ris:
    SAMTOOLS <- "samtools"
    STAR <- "star"
} else {
    SAMTOOLS <- "/bar/cfan/anaconda2/envs/jupyter/bin/samtools"
    STAR <- "/opt/apps/STAR/2.5.4b/STAR"
}
