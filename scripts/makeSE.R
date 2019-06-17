#' NOTE: It is important to make the pheno_type as a factor with levels named "FSHD" and "Control".
#' Otherwise the tools in the misc_tools.R would not work.
#' R version: R/3.4.0-foss-2016b-fh1
#' Note: pay attention to 32-0005b and 2567b, size factor and gene counts
# ml R/3.4.0-foss-2016b-fh1
library(GenomicAlignments)

pkgDir <- "~/tapscott/RNA-Seq/hg38.FSHD.biopsyI-IV"
source("~/tapscott/R_package/SeqPipelineTools/R/makeSE.R")
bamDir <- file.path(pkgDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", full.name=TRUE)
names(bamFiles) <- gsub(".bam", "", basename(bamFiles))
workers <- 10L

#'
#' tophat mapping results
#'
tophatFolders <- list.files(file.path(pkgDir, "tophat"), full.name=TRUE)
tophatFolders <- tophatFolders[!grepl("bam", tophatFolders)]
summary_files <- file.path(tophatFolders, "align_summary.txt")
tophat_summary <- lapply(summary_files, function(x) {
    n <- readLines(x)
    sapply(strsplit(n[2:3], ":"), "[[", 2)
})
names(tophat_summary) <- basename(tophatFolders)
tophat_summary <- as.data.frame(do.call(rbind, tophat_summary))
colnames(tophat_summary) <- c("Input", "Mapped")

#'
#' sample information 
#' 
si <- DataFrame(sample_name=sub(".bam", "", basename(bamFiles)),
                file_bam=bamFiles)
load("~/tapscott/RNA-Seq/hg38.FSHD.biopsy4/data/FirstVisit.SE.rda")
df <- colData(FirstVisit.SE)[, c("sample_id", "Batch", "pheno_type")]
si <- append(si, df[si$sample_name, ])
si <- append(si, as(tophat_summary[si$sample_name, ], "DataFrame"))

#'
#' use the comprehensive gencode v24 annotation from Gencode website,
#' hg38.HomoSapiens.Gencode.v24.
#' Don't use the knownGene track from UCSC, TxDb.Hsapiens.UCSC.hg38.GENCODE.v24.
#' It doesn't have the psudogenes we identified as biomarkers.
#' NOTE: This script is written in the way that is reproducible.
#' 
library(hg38.HomoSapiens.Gencode.v24)
library(org.Hs.eg.db)
library(BiocParallel)
library(GenomicAlignments)
mparam <- MulticoreParam(workers = workers, type = "SOCK")
features <- GenomicFeatures::exonsBy(hg38.HomoSapiens.Gencode.v24, by="gene")
se <- GenomicAlignments::summarizeOverlaps(features=features,
                                           reads=si$file_bam,
                                           mode="IntersectionStrict",
                                           inter.feature=TRUE,
                                           singleEnd=TRUE,
                                           ignore.strand=TRUE,
                                           BPPARAM=mparam)
colData(se) <- si
colnames(se) <- si$sample_name
save(se, file=file.path(pkgDir, "data", "se.rda"))

#' 
#' PCA
#' 
library(DESeq2)
library(ggplot2)
pkgDir <- "~/tapscott/RNA-Seq/hg38.FSHD.biopsyI-IV"
load(file.path(pkgDir, "data", "se.rda"))
figDir <- file.path(pkgDir, "figures")
dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > ncol(se), ], design = ~ pheno_type)
rlg <- rlog(dds)
data <- plotPCA(rlg, intgroup=c("pheno_type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
title <- "BatchI-IV"
qplot(PC1, PC2, color=pheno_type, data=data) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    labs(list(title=title))
ggsave(file=file.path(figDir, paste0(title, "_SampleDistancePCAplot.png")))

#'
#' patch the counts for 55 biomarkers (mk) and four major markers (fourmk)
#' 
data(gene.anno)
source(file.path(pkgDir, "scripts",  "misc_tools.R"))
library(xlsx)
sig <- getDUX4Targets()
i <- gene.anno$gene_name %in% as.character(sig$gene.name) ## use gene name

markers <- gene.anno$gene_name[i]
mk <- getMarkerCounts(markers, txdb=hg38.HomoSapiens.Gencode.v24,
                      bamFiles=si$file_bam, anno=gene.anno, txdb.type="Gencode",
                      workers=workers)
fourmk <- getMarkerCounts(markers=c("LEUTX", "PRAMEF2", "TRIM43", "KHDC1L"),
                          txdb=hg38.HomoSapiens.Gencode.v24,
                          bamFiles=si$file_bam, anno=gene.anno, txdb.type="Gencode",
                          workers=workers)

colData(mk) <- colData(fourmk) <- colData(Biopsy5.SE)
mk <- patchAnnotationSE(mk, gene.anno)
fourmk <- patchAnnotationSE(fourmk, gene.anno)
save(mk, file=file.path(pkgDir, "data", "mk.rda"))
save(fourmk, file=file.path(pkgDir, "data", "fourmk.rda"))

#'
#' immune markers counts
#' NOTE: D87024.2 around area of IGL cluster does not exist in Gencode v24
#' Chao-Jen added IGLC1 and IGLC2 cluster genes
#' 
markers <- c("AIRE", "CDKN2A", "DMBT1", "GCNT3", "IGHA1", "IGHA2", "IGHD",
             "IGHV3-30", "IGHV3-33", "IGLC1", "IGLC2", "IGLC3", "IGLV3-10",
             "PDCD1", "RGS1", "TDGF1", "TREM2")
immuneMarkers <- getMarkerCounts(markers, txdb=hg38.HomoSapiens.Gencode.v24,
                      bamFiles=si$file_bame, anno=gene.anno, txdb.type="Gencode",
                                 workers=workers)
immuneMarkers <- patchAnnotationSE(immuneMarkers, gene.anno)
save(immuneMarkers, file=file.path(pkgDir, "data", "immuneMarkers.rda"))

#'
#' Patch the markers to Biopsy4.SE
#' 
assays(se)[["counts"]][rownames(mk),] <- assays(mk)[["counts"]]
assays(se)[["counts"]][rownames(immuneMarkers),] <- assays(immuneMarkers)[["counts"]]
se <- patchAnnotationSE(se, gene.anno)
save(se, file=file.path(pkgDir, "data", "se.rda"))

#'
#' Counts for H3.X and H3.Y
#'
load(file.path(pkgDir, "data", "se.rda"))
gr <- GRanges(seqnames="chr5", IRanges(start=c(17491325, 17655128),
                  end=c(17491768, 17655538)))
mcols(gr)$exon_id <- mcols(gr)$exon_name <- c("H3.X", "H3.Y")
grl <- split(gr, gr$exon_id)
mcols(grl) <- DataFrame(gene_name=c("H3.X", "H3.Y"),
                        gene_type="Histone",
                        gene_status="NOVEL")
se_H3 <- GenomicAlignments::summarizeOverlaps(features=grl,
                                              reads=se$file_bam,
                                              mode="IntersectionStrict",
                                              inter.feature=TRUE,
                                              singleEnd=TRUE,
                                              ignore.strand=TRUE,
                                              BPPARAM=mparam)
colData(se_H3) <- colData(se)
colnames(se_H3) <- colnames(se)
save(se_H3, file=file.path(pkgDir, "data", "se_H3.rda"))

#'
#' sanitized dds: remove 0952, 32-0008-1, 32-0002b and controls 2374-1 and2374-2
#' 

library(DESeq2)
library(ggplot2)
pkgDir <- "~/tapscott/RNA-Seq/hg38.FSHD.biopsyI-IV"
load(file.path(pkgDir, "data", "se.rda"))
load(file.path(pkgDir, "data", "mk.rda"))
load(file.path(pkgDir, "data", "fourmk.rda"))

#' sanitize samples
rm_samples <- c("0952", "32-0008-1", "32-0002b", "2374-1", "2374-2")
se <-  se[, !colnames(se) %in% rm_samples]
sanitized.dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > 20, ],
                              design = ~ 1)
sanitized.dds <- estimateSizeFactors(sanitized.dds)
#' add FPKM to dds
source("~/tapscott/R_package/SeqPipelineTools/R/getScaledCounts.R")
bam_info <- SGSeq::getBamInfo(colData(sanitized.dds))
colData <- colData(sanitized.dds)
colData$read_length <- 100
colData$lib_size <- sapply(strsplit(as.character(colData$Mapped),
                                    " (", fixed=TRUE), "[[", 1)
colData$lib_size <- as.numeric(colData$lib_size)
fpkm <- getScaledCountsPerTx(sanitized.dds, read_length=colData$read_length,
                             lib_size=colData$lib_size)
assays(sanitized.dds)[["FPKM"]] <- fpkm
save(sanitized.dds, file=file.path(pkgDir, "data", "sanitized.dds.rda"))

#' sanitized.rlg
sanitized.rlg <- rlog(sanitized.dds)
save(sanitized.rlg, file=file.path(pkgDir, "data", "sanitized.rlg.rda"))

