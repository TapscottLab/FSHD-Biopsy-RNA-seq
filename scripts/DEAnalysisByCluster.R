#'
#' ml R/3.4.0-foss-2016b-fh1
#'
library(DESeq2)
library(xlsx)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
figDir <- file.path(pkgDir, "figures")
statsDir <- file.path(pkgDir, "stats")
source(file.path(pkgDir, "scripts", "misc_tools.R"))

biomarker.id <- c(PRAMEF2="ENSG00000120952.4",
                  TRIM43="ENSG00000144015.4",
                  KHDC1L="ENSG00000256980.4",
                  LEUTX="ENSG00000213921.6")
other.id <- c(DUXA="ENSG00000258873.2", ZSCAN4="ENSG00000180532.10")
#' NOTE: DUXA=ENSG00000258873.2, zscan4=ENSG00000180532.10

#'
#' loading some data
#' 
load(file.path(pkgDir, "data", "se.rda"))
load(file.path(pkgDir, "data", "sanitized.dds.rda"))
load(file.path(pkgDir, "data", "mk.rda"))
load(file.path(pkgDir, "data", "immuneMarkers.rda"))
historicalDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsy.discovery"
a <- get(load(file.path(historicalDir, "data", "discovery.SE.rda")))
load(file.path(historicalDir, "data_sampling", "discovery.SE.rda"))
#'
#' make sure the RNA.Cluster is added to colData of sanitized.dds
#' 
if (!any(names(colData(sanitized.dds)) %in% "RNA.Cluster")) {
    source(file.path(pkgDir, "scripts", "addColDataToDDS.R"))
    load(file.path(pkgDir, "data", "sanitized.dds.rda"))
}

#'
#' Combine cluster1 with the historical control and make a SummerizedExperiment instance
#' colData=sample_name, file_bam, pheno_type
#' 
all(rownames(se)==rownames(discovery.SE))

#' historical and exclude 2401
getHistoricalControl <- function(discovery.SE, cluster.SE) {
    historical_control <- discovery.SE[, discovery.SE$pheno_type=="Control"]
    historical_control <- historical_control[, -4] # exlude 2401
    colData(historical_control) <- colData(historical_control)[, c("sample_name", "file_bam",
                                                               "pheno_type", "lib_size")]
    historical_control$set_group <- "discovery"
    historical_control <- historical_control[rownames(cluster.SE)]
    historical_control
}
#' combine cluster with historical controls
combClusterWithHistoricalControl <- function(sanitized.dds, RNA.cluster="1", discovery.SE) {
    cluster <- sanitized.dds[, sanitized.dds$RNA.Cluster %in% RNA.cluster]
    lib_size <- sapply(strsplit(as.character(cluster$Mapped), " (", fixed=TRUE), "[[", 1)
    cluster$lib_size <- as.numeric(lib_size)
    colData(cluster) <- colData(cluster)[,  c("sample_name", "file_bam",
                                 "pheno_type", "lib_size")]
    cluster$set_group <- "new"

    historical_control <- getHistoricalControl(discovery.SE, cluster.SE=cluster)
    #' combine
    assays <- cbind(assays(cluster)[[1]], assays(historical_control)[[1]])
    colData <- rbind(colData(cluster), colData(historical_control))

    comb <- SummarizedExperiment(assays=SimpleList(counts=assays),
                                 rowRanges=rowRanges(cluster),
                                 colData=colData)
    pheno_type <- as.character(comb$pheno_type)
    i <- comb$set_group=="discovery"
    pheno_type[i] <- paste0("Discovery_", pheno_type[i])
    #pheno_type <- factor(pheno_type, levels=c("FSHD", "Control", "Discovery_Control"))
    comb$pheno_type <- factor(pheno_type)
    comb
}

comb <- combClusterWithHistoricalControl(sanitized.dds, RNA.cluster="1",
                                         discovery.SE=discovery.SE)

#'
#' PCA for cluster 1, control and historical control
#'
makePCAPlotByPhenoType(se=comb, figDir=figDir, title="PCA_Cluster1AndControls")

#'
#' (1) Cluster 1: FSHD vs Control from new set
#'
sub <- comb[, comb$set_group=="new"]
sub$pheno_type <- factor(sub$pheno_type)
makePCAPlotByPhenoType(se=sub, figDir=figDir, title="PCA_Cluster1FSHDAndControls")
res.1 <- do_deseq(se=sub, mk=mk, title="Cluster1FSHDvsControl",
                  logFC.threshold=2, alpha=0.05,
                  statsDir=statsDir)

## One biomarkers, ZSCAN4, is DE

#'
#' (2) Cluster 1 FSHD vs historical control
#'
idx <- comb$set_group == "discovery" | comb$pheno_type == "FSHD"
sub <- comb[, comb$pheno_type %in% c("FSHD", "Discovery_Control")]
sub$pheno_type <- factor(sub$pheno_type)

#' testing: compare size factor between using comb and comb2
#dds <- DESeqDataSet(sub, design = ~ pheno_type)
#dds <- estimateSizeFactors(dds)
#dds$sizeFactor

#sub2 <- comb2[, comb$pheno_type %in% c("FSHD", "Discovery_Control")]
#sub2$pheno_type <- factor(sub2$pheno_type)
#dds2 <- DESeqDataSet(sub2, design = ~ pheno_type)
#dds2 <- estimateSizeFactors(dds2)
#dds2$sizeFactor
#' end of testing

makePCAPlotByPhenoType(se=sub, figDir=figDir, title="PCA_Cluster1FSHDAndHistControls")
res.2 <- do_deseq(se=sub, mk=mk, title="Cluster1FSHDvsHistControl", logFC.threshold=2,
                  alpha=0.05,
                  statsDir=statsDir)

## PRAMEF15

#'
#' Compare res.1 (Cluster1 vs Control) and res.2 (Cluster1 vs historical control)
#'
sig.1 <- summaryDESeqResults(res.1, padj.threshold=0.05, lfc.threshold=2)
sig.2 <- summaryDESeqResults(res.2, padj.threshold=0.05, lfc.threshold=2)
sum(rownames(sig.1) %in% rownames(sig.2))/nrow(sig.1) #62%

library(VennDiagram)
library(ggplot2)
venn.diagram(x=list(Cluser1VsCont=rownames(sig.1), Cluser1VsHisCont=rownames(sig.2)),
             filename=file.path(pkgDir, "figures",
                 "Venn_Cluster1vsControl_HistControl.png"))
#' scatter plot of logFC of the union of DE
geneid <- union(rownames(sig.1), rownames(sig.2))
df <- data.frame(Cluster1vsCont=res.1[geneid, "log2FoldChange"],
                 Cluster1vsHist=res.2[geneid, "log2FoldChange"],
                 DE="Both", stringsAsFactors=FALSE)
rownames(df) <- geneid
df$gene_name <- rowData(comb[geneid])$gene_name
cont <- setdiff(rownames(sig.1), rownames(sig.2))
hist <- setdiff(rownames(sig.2), rownames(sig.1))
df$DE[rownames(df) %in% cont] <- "vsControl"
df$DE[rownames(df) %in% hist] <- "vsHistoric"
df$moderate_agree <- ifelse(abs(df[, 1]-df[, 2]) > 1, FALSE, TRUE)
df$moderate_agree <- df$moderate_agree | df$DE=="Both"
png(file.path(pkgDir, "figures", "Scatter_DE_Cluster1vsContvsHis1.png"))
ggplot(df, aes(x=Cluster1vsCont, y=Cluster1vsHist)) +
    geom_point(aes(col=DE))
dev.off()

#' scatter for the whole genes list
name <- intersect(rownames(res.1), rownames(res.2))
dfall <- data.frame(Cluster1vsCont=res.1[name, "log2FoldChange"],
                 Cluster1vsHist=res.2[name, "log2FoldChange"],
                    DE="NONE", stringsAsFactors=FALSE)
rownames(dfall) <- name
dfall$DE[name %in% cont] <- "vsControl"
dfall$DE[name %in% hist] <- "vsHistoric"
dfall$DE[name %in% intersect(rownames(sig.1), rownames(sig.2))] <- "Both"
png(file.path(pkgDir, "figures", "Scatter_DE_Cluster1vsContvsHis.png"))
ggplot(dfall, aes(x=Cluster1vsCont, y=Cluster1vsHist)) +
    geom_point(aes(col=DE, alpha=0.3)) +
        geom_smooth(method=lm)
dev.off()
#'
#' (3) Control vs historical control
#' 
sub <- comb[, comb$pheno_type != "FSHD"]
sub$pheno_type <- factor(sub$pheno_type)
makePCAPlotByPhenoType(se=sub, figDir=figDir, title="PCA_Cluster1ControlsAndHistControls")
res.3 <- do_deseq(se=sub, mk=mk, title="Cluster1ControlsvsHistControl", logFC.threshold=2,
                  alpha=0.05,
                  statsDir=statsDir)
res.3[biomarker.id,]
tmp <- summaryDESeqResults(res.3)

#'
#' (3.1) loner cluster 1 FSHD vs controls
#'
sub <- comb[, c(12, 19, 4:11)]
sub$pheno_type <- factor(sub$pheno_type)
tmp <- do_deseq(se=sub, mk=mk, title="Cluster1LonerVsControl", logFC.threshold=2,
                  alpha=0.05,
                statsDir=statsDir)
loner <- summaryDESeqResults(tmp) #' ZSCAN is DE

#'
#' (4) cluster2 vs control and historical controls
#'
comb1 <- combClusterWithHistoricalControl(sanitized.dds, RNA.cluster="1",
                                          discovery.SE=discovery.SE)

comb2 <- combClusterWithHistoricalControl(sanitized.dds, RNA.cluster="2",
                                          discovery.SE=discovery.SE)
#' get control from cluster 1, new set
contr <-  comb1[, comb1$pheno_type=="Control" & comb1$set_group=="new"]
comb <- cbind(comb2, contr)
comb$pheno_type <- factor(comb$pheno_type, levels=c("Control", "Discovery_Control",
                                                    "FSHD"))
makePCAPlotByPhenoType(se=comb, figDir=figDir, title="PCA_Cluster2AndAllControls")

#' cluster2FSHD vs control
sub <- comb[, comb$set_group == "new"]
sub$pheno_type <- factor(sub$pheno_type)
makePCAPlotByPhenoType(se=sub, figDir=figDir, title="PCA_Cluster2FSHDAndControls")
res2.1 <- do_deseq(se=sub, mk=mk, title="Cluster2FSHDvsControl", logFC.threshold=2,
                  alpha=0.05,
                   statsDir=statsDir)
## three biomarkers, ZSCAN4, CCNA1, KHDC1L biomarkers are DE

#' active cluster2FSHD (active disease 2358, 007, 0013) vs control
sub2 <- sub[, -c(1, 3)] # exclude 1614, 2411
res2.3 <- do_deseq(se=sub2, mk=mk, title="Cluster2ActiveFSHDvsControl", logFC.threshold=2,
                  alpha=0.05,
                   statsDir=statsDir)
tmp3 <- summaryDESeqResults(res2.3)
tmp1 <- summaryDESeqResults(res2.1)
sum(rownames(tmp1) %in% rownames(tmp3)) # 254 (out of 306) are in tmp3 (549)
diff <- tmp3[!rownames(tmp3) %in% rownames(tmp1), ]
write.csv(as.data.frame(diff), file=file.path(statsDir, "DESeqResutls_Cluster2_DIFFbetweenActiveAndNonActiveFSHD.csv"))

#' is there any fibrosis genes in the diff list? hows the immune marker genes look like?
sub <- comb[, comb$set_group == "new"]
dds <- DESeqDataSet(sub, design = ~ pheno_type)
dds <- estimateSizeFactors(dds)
rlg <- rlog(dds)
immune <- rlg[rownames(immuneMarkers), ]
rownames(immune) <- mcols(rowRanges(immune))$gene_name
immune$active_status <- c("NO", "Yes", "NO", "Yes", "Yes", rep("Control", 9))
#' pheatmap
library(pheatmap)
annotation_col <- data.frame(active_status=immune$active_status)
rownames(annotation_col) <- colnames(immune)

pheatmap(assay(immune), annotation_col=annotation_col,
         main="Immune genes: cluster 2 FSHD vs Control",
         filename=file.path(pkgDir, "figures", "heatmap_immnue_Cluster2.pdf"))


#'
#' cluster2FSHD vs historial control
#' 
sub <- comb[, comb$pheno_type %in% c("FSHD", "Discovery_Control")]
sub$pheno_type <- factor(sub$pheno_type)
makePCAPlotByPhenoType(se=sub, figDir=figDir, title="PCA_Cluster2FSHDAndHistControls")
res2.2 <- do_deseq(se=sub, mk=mk, title="Cluster2FSHDvsHistControl", logFC.threshold=2,
                  alpha=0.05,
                   statsDir=statsDir)
tmp <- summaryDESeqResults(res2.2)

#' cluster2FSHD (active disease 2358, 007, 0013) vs control



#'
#' (5) cluster 3 vs control
#' 
idx <- sanitized.dds$RNA.Cluster == "3" | sanitized.dds$pheno_type == "Control"
sub <- sanitized.dds[, idx]
makePCAPlotByPhenoType(se=sub, figDir=figDir, title="PCA_Cluster3FSHDAndControls")
res3 <- do_deseq(se=sub, mk=mk, title="Cluster3FSHDvsControl", logFC.threshold=2,
                  alpha=0.05,
                  statsDir=statsDir)
## 47 biomarkers are DE

#'
#' (6) cluster 4 vs control
#' 
idx <- sanitized.dds$RNA.Cluster == "4" | sanitized.dds$pheno_type == "Control"
sub <- sanitized.dds[, idx]
makePCAPlotByPhenoType(se=sub, figDir=figDir, title="PCA_Cluster4FSHDAndControls")
res4 <- do_deseq(se=sub, mk=mk, title="Cluster4FSHDvsControl", logFC.threshold=2,
                  alpha=0.05,
                  statsDir=statsDir)
## 54 biomarkers are DE

#'
#' Visualization of some biomarkers and DE biomarkers
#' 

#' PCA by clusters / pheno_type
tmp <- paste0("Cluster", sanitized.dds$RNA.Cluster, "_", sanitized.dds$pheno_type)
pheno_type <- factor(tmp, levels=c("Cluster1_Control", "Cluster1_FSHD",
                       "Cluster2_FSHD", "Cluster3_FSHD", "Cluster4_FSHD"))
tmp <- sanitized.dds
tmp$pheno_type <- pheno_type
makePCAPlotByPhenoType(se=tmp, figDir=figDir, title="PCA_Clusters")

#' boxplots of four biomarkers and zscan and DuxA?
pdf(file.path(figDir, "BoxPlot_SomeBiomarkers.pdf"))
lapply(c(biomarker.id, other.id), function(x) {
     perGeneBoxplotCount(x, sanitized.dds)

})
dev.off()

#'
#' parallel plots
#' 
library(lattice)
tmp <- paste0("Cluster", sanitized.dds$RNA.Cluster, "_", sanitized.dds$pheno_type)
pheno_type <- factor(tmp, levels=c("Cluster1_Control", "Cluster1_FSHD",
                       "Cluster2_FSHD", "Cluster3_FSHD", "Cluster4_FSHD"))
tmp <- sanitized.dds
tmp$pheno_type <- pheno_type

#'
#' four biomarkers (plus DUXA and ZSCAN4)
#' 
data <- sapply(levels(tmp$pheno_type), function(x) {
    data <- counts(tmp, normalized=TRUE)[c(biomarker.id, other.id), tmp$pheno_type==x]
    rowMeans(data)
})

data <- as.data.frame(data)
data$gene_name <- c(names(biomarker.id), names(other.id))
pdf(file.path(figDir, "ParallelPlot_SomeBiomarkersCounts.pdf"))
parallelplot( ~data[c(1:5)], groups=gene_name,
                 panel=function(...) {
                      panel.parallel(..., common.scale=TRUE)
                 },
                 data=data,lty=1, #col="grey",
             ylab="normalized counts", #alpha=0.5,
             horizontal.axis=FALSE,
             main="Some Biomarkers",
             auto.key=list(space="top", columns=3))
dev.off()
#'
#' 54 biomarkers
#'
id <- intersect(rownames(mk), rownames(tmp))
data <- sapply(levels(tmp$pheno_type), function(x) {
    data <- counts(tmp, normalized=TRUE)[id, tmp$pheno_type==x]
    rowMeans(data)
})

data <- as.data.frame(data)
pdf(file.path(figDir, "ParallelPlot_54BiomarkersCounts.pdf"))
parallelplot( ~data[c(1:5)], 
                 panel=function(...) {
                      panel.parallel(..., common.scale=TRUE)
                 },
                 data=data,lty=1, col="grey",
             ylab="normalized counts", alpha=0.5,
             horizontal.axis=FALSE,
             main="54 Biomarkers")
dev.off()

#'
#' immune markers
#'
id <- intersect(rownames(immuneMarkers), rownames(tmp))
data <- sapply(levels(tmp$pheno_type), function(x) {
    data <- counts(tmp, normalized=TRUE)[id, tmp$pheno_type==x]
    rowMeans(data)
})

data <- as.data.frame(data)
pdf(file.path(figDir, "ParallelPlot_ImmuneMarkerCounts.pdf"))
parallelplot( ~data[c(1:5)], 
                 panel=function(...) {
                      panel.parallel(..., common.scale=TRUE)
                 },
                 data=data,lty=1, col="grey",
             ylab="normalized counts", alpha=0.5,
             horizontal.axis=FALSE,
             main="Immune Markers")
dev.off()

#' pheatmap
dds <- DESeqDataSet(tmp, design = ~ pheno_type)
rlg <- rlog(dds)
anno_col <- data.frame(cluster=rlg$pheno_type)
rownames(anno_col) <- colnames(tmp)
cnt <- rlg[id, ]
rownames(cnt) <- mcols(rowRanges(cnt))$gene_name
pheatmap(assay(cnt), annotation_col=anno_col,
         filename=file.path(pkgDir, "figures", "heatmap_immune_markers.pdf"))
