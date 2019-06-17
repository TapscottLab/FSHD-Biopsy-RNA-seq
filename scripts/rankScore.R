#'
#' R/3.4.0-foss-2016b-fh1
#'
library(DESeq2)
library(ggplot2)

pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
load(file.path(pkgDir, "data", "sanitized.dds.rda"))
load(file.path(pkgDir, "data", "se.rda"))
load(file.path(pkgDir, "data", "mk.rda"))
load(file.path(pkgDir, "data", "fourmk.rda"))


#' get rank scores
library(lattice)
library(pheatmap)
library(gridExtra)
library(xlsx)
#rlg <-  rlog(sanitized.dds)

source(file.path(pkgDir, "scripts", "misc_tools.R"))
rank.fourmk<- getRankScore(dds=sanitized.dds, which.genes=rownames(fourmk),
                           main="FourBiomarkers_AllFirstVisit",  num.cluster=4)
## remove "ENSG00000144188.9": TRIM43CP - only have 21 reads across samples
mk <- mk[!rownames(mk) %in% "ENSG00000144188.9", ]
rank.allmk <- getRankScore(dds=sanitized.dds, which.genes=rownames(mk),
                           main="AllBiomarkers_AllFirstVisit",  num.cluster=5)

desDir <- file.path(pkgDir, "stats")
main <- "FourBiomarkers_AllFirstVisit"
res <- getRankSummary(dds=sanitized.dds, rank=rank.fourmk, main=main,
                      desDir=desDir, sep=FALSE, generate.report=TRUE)
main <- "AllBiomarkers_AllFirstVisit"
res2 <- getRankSummary(dds=sanitized.dds, rank=rank.allmk, main=main,
                       desDir=desDir, sep=TRUE, generate.report=TRUE)

#'
#' figures for manuscript - scatter plot again with better cluster number...
#'

