#'
#' Set up Parameters:
#' ml R/3.4.0-foss-2016b-fh1
#'
pkgDir <- "~/tapscott/RNA-Seq/hg38.FSHD.biopsyI-IV"
source(file.path(pkgDir, "scripts", "viz_tools.R"))
source(file.path(pkgDir, "scripts", "misc_tools.R"))

library(DESeq2)
library(lattice)
library(pheatmap)
library(parallel)
library(gridExtra) #grid.arrange(obj1, obj2, ncol=2)

load(file.path(pkgDir, "data", "sanitized.dds.rda"))
load(file.path(pkgDir, "data", "sanitized.rlg.rda"))
load(file.path(pkgDir, "data", "se_H3.rda"))
load(file.path(pkgDir, "data", "mk.rda"))
load(file.path(pkgDir, "data", "fourmk.rda"))
rownames(sanitized.dds) <- .cleanGeneID(rownames(sanitized.dds))
rownames(sanitized.rlg) <- .cleanGeneID(rownames(sanitized.rlg))
rownames(fourmk) <- .cleanGeneID(rownames(fourmk))
rownames(mk) <-.cleanGeneID(rownames(mk))
biomarkers <- rownames(mk)
#' add clinical data to colData of sanitized.dds and sanitized.dds
colData <- addClinicalsToColData(sanitized.dds)
colData(sanitized.dds) <- colData(sanitized.rlg) <- colData

tmp <- paste0("Cluster", sanitized.dds$RNA.Cluster, "_", sanitized.dds$pheno_type)
fac <- factor(tmp, levels=c("Cluster1_Control", "Cluster1_FSHD",
                       "Cluster2_FSHD", "Cluster3_FSHD", "Cluster4_FSHD"))
sanitized.rlg$fac <- sanitized.dds$fac <- fac

#'
#' H3.X and H3.Y are robust dux4-biomarkers
#'

H3.dds <- DESeqDataSet(se_H3, design = ~1)
H3.dds <- H3.dds[, colnames(sanitized.dds)]
colData(H3.dds) <- colData(sanitized.dds)

add_on <- SummarizedExperiment(assays=SimpleList(counts=rbind(counts(sanitized.dds),
                                                     counts(H3.dds))),
                               colData=colData(sanitized.dds),
                               rowData=rbind(rowData(sanitized.dds),
                                   rowData(H3.dds)))
for (cluster in 1:4) 
    print(counts(H3.dds)[, H3.dds$RNA.Cluster==cluster])

#'
#' GESA: controls vs cluster 4
#' 

getUpSigEachCluster <- function(sanitized.dds, cluster) {
    ## cluster = c("1", "2", "3", "4")
    idx <- sanitized.dds$RNA.Cluster == cluster | sanitized.dds$pheno_type == "Control"
    sub <- sanitized.dds[, idx]
    design(sub) <- ~ pheno_type
    dds <- DESeq(sub)
    res <- results(dds, alpha=0.05, lfcThreshold=2)
    res$gene_name <- rowData(dds)$gene_name
    up_sig <- subset(res, padj < 0.05 & log2FoldChange > 0)
}

#' (1)find robust, differentially expressed genes first
idx <- sanitized.dds$RNA.Cluster == "4" | sanitized.dds$pheno_type == "Control"
sub <- sanitized.dds[, idx]
design(sub) <- ~ pheno_type
dds <- DESeq(sub)
res <- results(dds, alpha=0.05, lfcThreshold=2)
res$gene_name <- rowData(dds)$gene_name
up_sig <- subset(res, padj < 0.05 & log2FoldChange > 0)
#' how do they look like in cluster 1, 2, and 3?

#' Use differentially expressed genes to find the GO terms that are
#' altered.
library(goseq)
library(org.Hs.eg.db)
library(GO.db)
universe <- rownames(res)
isDEGs <- as.integer(universe %in% rownames(up_sig))
names(isDEGs) <- sapply(strsplit(universe, ".", fixed=TRUE), "[[", 1)

library(hg38.HomoSapiens.Gencode.v24)
txsByGene <- transcriptsBy(hg38.HomoSapiens.Gencode.v24, by="gene")
names(txsByGene) <- sapply(strsplit(names(txsByGene), ".", fixed=TRUE), "[[", 1)
lengthData <- median(width(txsByGene))
bias.data <- lengthData[names(isDEGs)]
pwf <- nullp(isDEGs, bias.data=bias.data, plot.fit=FALSE)
GO.BP <- goseq(pwf, "hg38", "ensGene", test.cats=c("GO:BP"))
GO.BP$padj <- p.adjust(GO.BP$over_represented_pvalue, method="BH")
enriched.BP <- subset(GO.BP, padj < 0.01)
write.csv(enriched.BP, file=file.path(pkgDir, "stats", "GOBP_Cluster4_vs_Control.csv"))

#' display
for(go in enriched.BP$category[1:10]){
    print(GOTERM[[go]])
    cat("--------------------------------------\n")
}

#'
#' find out the genes  that's in immune respone, inflammatory respone and defense response
#' 
cat2DEgenes <- function(GOID) {
    gene2cat <- getgo(rownames(pwf), "hg38", "ensGene", fetch.cats = "GO:BP")
    names(gene2cat) <- rownames(pwf)
    cat2gene <- goseq:::reversemapping(gene2cat)
    #' sanity check
    doesIDexist <- GOID %in% names(cat2gene)
    if (!all(doesIDexist)) stop("GOID is not found")
    sig_gene <- rownames(pwf)[pwf$DEgenes==1]
    geneInGOID <- cat2gene[[GOID]]
    sig_gene[sig_gene %in% geneInGOID]
}

immune <- cat2DEgenes("GO:0006955")
inflamm <- cat2DEgenes("GO:0006954")

catGenes <- lapply(enriched.BP$category, cat2DEgenes)
names(catGenes) <- enriched.BP$category
GSEA_results <- list(enriched.BP=enriched.BP, catGenes=catGenes)
save(GSEA_results,
     file=file.path(pkgDir, "data", "GSEA_results.rda"))

#'
#' what genes are in the category
#' ceel differentiation, cell proliferation, multicellular
load(file.path(pkgDir, "data", "GSEA_results.rda"))
catGenes <- GSEA_results[[2]]
all <- unique(c(unlist(catGenes), biomarkers))
in_four_cat <- unique(c(immune=catGenes[["GO:0006955"]],
                         inflamm=catGenes[["GO:0006954"]],
                         cell=catGenes[["GO:0030154"]],
                         biomarkers))

data <- c(inflammatory_response=57,
          immune_process=108,
          cell_differentiation=122,
          robust_dux4_targets=55,
          in_other_BP=length(all)-length(in_four_cat),
          not_in_BP=nrow(up_sig) - length(all))
df <- data.frame(Category=names(data), numInCategory=data)
df
ggplot(data=df, aes(x=Category, y=numInCat)) +
    geom_bar(stat="identity", color="steelblue") +
        coord_flip()

sig <- names(isDEGs[isDEGs==1])
not_in_enriched_cat <- sig[!sig %in% all]
gene.anno[not_in_enriched_cat, "gene_name"]

#' most of them are psudo-genes, TRIM, PRAMEF, IGHG cluster genes

#'
#' most informative genes in immune, inflamm, cell differentiation category
#' 

#'
#' heatmap: use rlog
#' 
source(file.path(pkgDir, "scripts", "viz_tools.R"))
rlg <- rlog(sanitized.dds)
#' immune
filename <- file.path(pkgDir, "figures", "Heatmap_ImmuneDE_Cluster4.pdf")
biomarker.rlg <- rlg[immune, ]
fac <- list(RNA.Cluster=rlg$RNA.Cluster,
            N.R.Inflamm=rlg$Necrosis.Regeneration.Inflammation,
            Disease=rlg$pheno_type)
plotHeatmap(biomarker.rlg, fac=fac, main="Immune DE",
            filename=filename)
#' inflamm
filename <- file.path(pkgDir, "figures", "Heatmap_InflammDE_Cluster4.pdf")
biomarker.rlg <- rlg[inflamm, ]
plotHeatmap(biomarker.rlg, fac=fac, main="Inflamm DE",
            filename=filename)

#' extracellular matrix organization
extracellular <- catGenes[["GO:0030198"]]
biomarker.rlg <- rlg[extracellular, ]
fac <- list(Fibrosis=rlg$Interstitial.Fibrosis,
            Disease=rlg$pheno_type, 
            STIR=rlg$STIR, RNA.Cluster=rlg$RNA.Cluster)
filename <- file.path(pkgDir, "figures", "Heatmap_extracellular_matrix.pdf")
plotHeatmap(biomarker.rlg, fac=fac, main="Extracellular DE",
            filename=filename)

#' cell death - many of them are dux4-targeted genes
cell_death <- catGenes[["GO:0060548"]]
biomarker.rlg <- rlg[cell_death, ]
fac <- list(Disease=rlg$pheno_type, 
            RNA.Cluster=rlg$RNA.Cluster)
filename <- file.path(pkgDir, "figures", "Heatmap_negative_regulation_cell_death.pdf")
plotHeatmap(biomarker.rlg, fac=fac, main="Neg regulation cell death DE",
            filename=filename)

#' Dux4 robust markers
biomarkers <- intersect(biomarkers, rownames(rlg))
biomarker.rlg <- rlg[biomarkers, ]
fac <- list(Disease=rlg$pheno_type, 
            RNA.Cluster=rlg$RNA.Cluster)
filename <- file.path(pkgDir, "figures", "Heatmap_Dux4_markers.pdf")
plotHeatmap(biomarker.rlg, fac=fac, main="54 Dux4 biomarkers",
            filename=filename)

#' four DUX4 biomarkers
biomarker.rlg <- rlg[biomarkers, ]
fac <- list(Disease=rlg$pheno_type, 
            RNA.Cluster=rlg$RNA.Cluster)
filename <- file.path(pkgDir, "figures", "Heatmap_Dux4_markers.pdf")
plotHeatmap(biomarker.rlg, fac=fac, main="54 Dux4 biomarkers",
            filename=filename)

#' four DUX4 biomarkers
biomarker.rlg <- rlg[rownames(fourmk), ]
fac <- list(Disease=rlg$pheno_type, 
            RNA.Cluster=rlg$RNA.Cluster)
filename <- file.path(pkgDir, "figures", "Heatmap_Four_Dux4_markers.pdf")
plotHeatmap(biomarker.rlg, fac=fac, main="Four Dux4 biomarkers",
            filename=filename)
#'
#' parallel: 
#'
load(file.path(pkgDir, "data", "GSEA_results.rda"))
source(file.path(pkgDir, "scripts", "viz_tools.R"))
catGenes <- GSEA_results[[2]]

#' (a) H3.X, H3.Y and biomarkers
biomarker.dds <- add_on.dds[c(rownames(fourmk), "ENSG00000180532", "H3.X", "H3.Y"), ]
pp <- plotParallel(biomarker.dds, fac=fac, main="Some Dux4 Biomarkers",
                   key=TRUE, key.columns=4, FPKM=FALSE)
pdf(file.path(pkgDir, "figures", "Parallel_SomeBiomarkers.pdf"))
plot(pp)
dev.off()
   
#' (b) other enriched go terms: put all together

dux4_mk <- intersect(rownames(mk), rownames(sanitized.dds))
markers <- list(dux4=dux4_mk,
    inflamm=catGenes[["GO:0006954"]],
    immune=catGenes[["GO:0006955"]],
    extracelluler=catGenes[["GO:0030198"]],
    cell_death=catGenes[["GO:0060548"]])

getMedianFPKMEachFactor <- function(biomarkers, fac) {
    biomarker.dds <- sanitized.dds[biomarkers]
    data <- sapply(levels(fac), function(x) {
            data <- assays(biomarker.dds)[["FPKM"]][, fac==x]
            rowMedians(data)})
    data <- as.data.frame(data)
    rownames(data) <- biomarkers
    data
}
med_fpkm <- lapply(markers, getMedianFPKMEachFactor, fac)
med_fpkm <- do.call(rbind, med_fpkm)
med_fpkm$BP <- sapply(strsplit(rownames(med_fpkm), ".", fix=TRUE), "[[", 1)
parallelplot( ~ med_fpkm[1:length(levels(fac))]|BP,
                       panel=function(...) {
                          panel.parallel(..., common.scale=TRUE)
                       },
                       data=med_fpkm, lty=1, col="grey",
                       ylab="Median (FPKM)", alpha=0.5,
                       horizontal.axis=FALSE)

getInformativeDE <- function(biomarkers, top=5) {
    #' find most informative DEs with highest SD across samples
    sub <- sanitized.dds[biomarkers]
    sds <- rowSds(assays(sub)[["FPKM"]])
    sorted <- biomarkers[order(sds, decreasing=TRUE)]
    top_markers <- sorted[1:top]
    print(rowData(sub[top_markers])$gene_name)
    return(top_markers)
}

informative_markers <- unique(unlist(lapply(markers[2:5], getInformativeDE)))
rowData(sanitized.dds[informative_markers])$gene_name
#rm <- which(informative_markers %in% c("ENSG00000164692", "ENSG00000108821" ))
#nformative_markers <- informative_markers[-rm]

#' for each enriched bp, plot the parallel plot for the most informative genes
pdf(file.path(pkgDir, "figures", "Parallel_SomeEnrichedGoTerms.pdf"))
for (i in 2:5) {
    p <- plotParallel(sanitized.dds[markers[[i]]], fac,
                 main=names(markers)[i],
                 key=FALSE, FPKM=TRUE)
    plot(p)
    genes <- getInformativeDE(markers[[i]])
    p <- plotParallel(sanitized.dds[genes], fac,
                 main=paste0(names(markers)[i], " with highest SD"),
                      key=TRUE, FPKM=TRUE)
    plot(p)
}
dev.off()

plotParallel(sanitized.dds[informative_markers], fac=fac,
             main="Most informative within enriched BPs",
             key=TRUE, FPKM=TRUE)

p <- plotHeatmap(sanitized.rlg[c(informative_markers, rownames(fourmk))],
                 fac=list(factor  =sanitized.rlg$fac,
                          fibrosis=sanitized.rlg$Interstitial.Fibrosis),
                 main="Most informative markers", scale="row",
                 filename=file.path(pkgDir, "figures",
                     "Heatmap_InterestingMarkers.pdf"))
plot(p$gtable)

pdf(file.path(pkgDir, "figures", "Boxplot_InformativeMarkers.pdf"))
for (id in c(rownames(fourmk), informative_markers)) {
    p <- perGeneBoxplotCount(id, sanitized.dds)
    plot(p)
}
dev.off()
            
#' COL1A1, COL1A2, COL3A1 are very similar.


#'
#' GSEA for Cluster vs Control
#'
#' NOTE: There are two groups of cluster 1 samples, The 2567 and 0016 are
#' highly expressed in immune response, inflammatory process, and extracellular
#' structure.
#'
#' (1) cluster 1 (except 267, 32-0018 and 2583, nothting come up
idx <- (sanitized.dds$RNA.Cluster==1 |
            sanitized.dds$pheno_type == "Control") &
                !colnames(sanitized.dds) %in% c("2567", "32-0016", "2583")
enriched_bp <- do_goseq(sanitized.dds, idx)
#' (2) cluster 1 2567 and 32-0018: extracellular structures 48/295
#' collagen fibril organization 14/38
idx <- colnames(sanitized.dds) %in% c("2567", "32-0016") |
    sanitized.dds$pheno_type == "Control"
enriched_bp <- do_goseq(sanitized.dds, idx)
write.csv(enriched_bp, file=file.path(pkgDir, "stats", "GOBP_2567And0016_vs_Control.csv"))

#' (3) cluster 2 (active) vs. control: nothing, all cluster 2 vs control: nothing
idx <- colnames(sanitized.dds) %in% c("2358", "32-0013") |
    sanitized.dds$pheno_type == "Control"
enriched_bp <- do_goseq(sanitized.dds, idx)
#' no enrichment for cluster 2 vs. control

#' (4) cluster 3 vs control:negative regulation of apoptotic process,
#' negative regulation of programmed cell death
idx <- (sanitized.dds$RNA.Cluster==3 |
            sanitized.dds$pheno_type == "Control")
enriched_bp <- do_goseq(sanitized.dds, idx)
write.csv(enriched_bp, file=file.path(pkgDir, "stats", "GOBP_Cluster3_vs_Control.csv"))


#'
#' Compare UP_SIG list to Cluster 4
#'
library(VennDiagram)

getUpSigEachCluster <- function(cluster, sanitized.dds, alpha=0.05,
                                lfcThreshold=2) {
    ## cluster = c("1", "2", "3", "4")
    message("Cluster ", cluster)
    idx <- sanitized.dds$RNA.Cluster == cluster | sanitized.dds$pheno_type == "Control"
    sub <- sanitized.dds[, idx]
    design(sub) <- ~ pheno_type
    dds <- DESeq(sub)
    res <- results(dds, alpha=alpha, lfcThreshold=lfcThreshold)
    res$gene_name <- rowData(dds)$gene_name
    #' return up-regulated differentially expressed genes only
    up_sig <- subset(res, padj < alpha & log2FoldChange > 0)
    up_sig <- up_sig[order(up_sig$pvalue), ]
}

filterOnLFC <- function(dds.res, threshold) {
    subset(dds.res, log2FoldChange > threshold)
}

upSigList <- lapply(levels(sanitized.dds$RNA.Cluster), getUpSigEachCluster,
                    sanitized.dds, alpha=0.01, lfcThreshold=0)
upSigList <- lapply(upSigList, filterOnLFC, threshold=2)
names(upSigList) <- levels(sanitized.dds$RNA.Cluster)
elementNROWS(upSigList)

sum(rownames(upSigList[[3]]) %in% rownames(upSigList[[4]]))
sum(rownames(upSigList[[2]]) %in% rownames(upSigList[[4]]))
sum(rownames(upSigList[[1]]) %in% rownames(upSigList[[4]]))

venn <- venn.diagram(list(cluser1=rownames(upSigList[[1]]),
                          cluster4=rownames(upSigList[[4]])), 
             filename=file.path(pkgDir, "figures",
                                "VennDiam_DEGs_EachClusterFSHDvsControl.png"))


bm <- c(RPG4="ENSG00000116690", rownames(fourmk)) 
plotParallel(sanitized.dds[bm], fac=fac,
             main="RPG4",
             key=TRUE, FPKM=TRUE)

bm <- rownames(upSigList[[1]])[1:10]
tmp <- bm[c(1, 5, 7, 4)]
plotParallel(sanitized.dds[tmp], fac=fac,
             main="RPG4",
             key=TRUE, FPKM=TRUE)

## vinndiagram: use alpha=0.05 and cutoff lfc=2, and plot a vendiagram

## parallel: all DE in cluster 1 and some important genes: PRG4, TNMD, ...
plotParallel(sanitized.dds[rownames(upSigList[[1]])], fac=fac,
             main="DE Cluster 1 vs Control", key=FALSE, FPKM=TRUE)
