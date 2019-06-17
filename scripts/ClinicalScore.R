library(DESeq2)
library(xlsx)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
figDir <- file.path(pkgDir, "figures")
statsDir <- file.path(pkgDir, "stats")
source(file.path(pkgDir, "scripts", "misc_tools.R"))
load(file.path(pkgDir, "data", "sanitized.dds.rda"))

#'
#' Zscan4 vs pathology score
#'
rlg <- rlog(sanitized.dds)
clinic <- read.xlsx(file.path(pkgDir, "extdata", "Tawil-Master-4-13-17.xlsx"),
          sheetIndex=1, rowIndex=c(2:38), colIndex=c(1, 19))

rownames(clinic) <- clinic[, 1]
rownames(clinic) <- sub("01--", "", rownames(clinic))
rownames(clinic) <- sub("01-", "", rownames(clinic), fixed=TRUE)

sub <- rlg[, rlg$pheno_type=="FSHD"]
zscan4 <- assay(sub)["ENSG00000180532.10", ]
names(zscan4)[2] <- "0952"
clinic$zscan4 <- zscan4[rownames(clinic)]
rm <- c("2320", "2393", "32-0011", "32-0016")
clinic <- clinic[!rownames(clinic) %in% rm, ] #0.55
