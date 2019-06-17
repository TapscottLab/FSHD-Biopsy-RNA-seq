#' This script calls up do_tophat_D4Z4.sh and align ../trimmed/*.fastq
#' files to D4Z4 (chr4) regions.
#' 
#'  ml R/3.4.0-foss-2016b-fh1
#'

pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
trimmedDir <- file.path(pkgDir, "trimmed")
fqFiles <- list.files(trimmedDir, pattern=".trim.fastq",
                      full.name=TRUE)
tophatDir <- file.path(pkgDir, "tophat_D4Z4")
setwd(file.path(pkgDir, "scripts"))

#' align filtered/trimmed
for (i in 1:length(fqFiles)) {
    cmd <- sprintf("sbatch -n4 ./do_tophat_D4Z4.sh %s %s",
                   fqFiles[i], tophatDir)
    system(cmd) #54971067
}


#' make sure all the tophat jobs have completed
bamDir <- file.path(pkgDir, "tophat_D4Z4", "bam")
tmp <- list.files(bamDir, pattern=".txt")
length(tmp)
