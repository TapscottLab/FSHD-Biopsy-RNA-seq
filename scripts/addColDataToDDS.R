#'
#' ml R/3.4.0-foss-2016b-fh1
#'
#' This script add column data, RNA.Cluster and other clinical data, to the
#' sanitized.dds instance.
#' 
library(DESeq2)
library(xlsx)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
source(file.path(pkgDir, "scripts", "misc_tools.R"))
load(file.path(pkgDir, "data", "sanitized.dds.rda"))

#'
#' Add FSHD-specific biomarker rank (RNA.Cluster)
#' 
rank4 <- read.xlsx(file=file.path(pkgDir, "stats",
                       "RankSummary_FourBiomarkers_AllFirstVisit.xlsx"),
                   sheetIndex=1)

#' fix cluster label
RNA.Cluster <- rep(1, nrow(rank4))
RNA.Cluster[rank4$cluster==2] <- 1
RNA.Cluster[20:24] <- 2
RNA.Cluster[rank4$cluster==1] <- 3
RNA.Cluster[rank4$cluster %in% c(3,4)] <- 4
RNA.Cluster <-  factor(RNA.Cluster, levels=c(1, 2, 3, 4))
names(RNA.Cluster) <- as.character(rank4[, 1])

if (!any(names(colData(sanitized.dds)) %in% "RNA.Cluster")) {
    sanitized.dds$RNA.Cluster <- RNA.Cluster[colnames(sanitized.dds)]
    save(sanitized.dds, file=file.path(pkgDir, "data", "sanitized.dds.rda"))
}

colData <- addClinicalsToColData(sanitized.dds)
