#' 
#' ml R/3.4.3-foss-2016b-fh2   
#' Supplemental figure 1A and 1B - heatmap of four and 54 DUX4 biomarkers
#' Supplemental Figure 2: boxplot of four biomarkers
#' Supplemental Figure 3: PCA plot of all samples
#' Supplemental Figure 4: genes that show high expression in two group 1 sample
#' 01-0037 and 32-0016
#' Supplemental Figure 5: genes that show high expression in two group 1 sample
#' 01-0037 and 32-0016



#'
#' tools
#' 
makeBioMarkerHeatmap <- function(dds, file_name) {
    #' dds is the biomarkers container
    rownames(dds) <- rowData(dds)$gene_name
    tmp <- rlog(dds, fitType="local")
    lg <- assay(tmp)
    avg <- rowMeans(lg)
    centered.lg <- lg - matrix(avg, nrow=length(avg), ncol=ncol(lg))
    annotation_col <- data.frame(pheno_type=dds$pheno_type)
    rownames(annotation_col) <- colnames(centered.lg)
    fontsize_row <- 10
    fontsize_col <- 10
    if (nrow(centered.lg) > 10) fontsize_row <- 7
    if (ncol(centered.lg) > 30) fontsize_col <- 8

    pheatmap(centered.lg,
             annotation_col=annotation_col, silent=TRUE,
             main="", fontsize_row=fontsize_row,
             fontsize_col=fontsize_col, 
             file=file_name)
}

#'
#' Supplemental Figure 1A and 1B - heapmap of biomarkers
#' 
dds <- sanitized.dds[rownames(fourmk)]
file_name <- file.path(figDir, "SupFig1B_Heatmap_FourBioMarkers.pdf")
makeBioMarkerHeatmap(dds, file_name)

i <- rownames(sanitized.dds) %in% rownames(mk)
dds <- sanitized.dds[i]
file_name <- file.path(figDir, "SupFig1A_Heatmap_AllBioMarkers.pdf")
)makeBioMarkerHeatmap(dds, file_name)

#'
#' Supplemental Figure 2: boxplot of four biomarkers 
#' (count, and RPKM)
#'

#data <- counts(sanitized.dds[rownames(fourmk)], normalized=TRUE)
data <- assays(sanitized.dds[rownames(fourmk)])$FPKM
data <- melt(data)
colnames(data) <- c("id", "sample_name", "value")
data$group <- rep(sanitized.dds$group, each=4)
symbol <- rowData(sanitized.dds[rownames(fourmk)])$gene_name
data$symbol <- rep(symbol, ncol(sanitized.dds))
gg <- ggplot(data, aes(x=group, y=value)) +
    geom_boxplot(outlier.colour=NA) +
        geom_jitter(width=0.2, size=0.6) +
            theme_bw() + 
                facet_wrap(~symbol, ncol=2, scales="free_y") +
                    ylab("RPKM") + xlab("sample group") + 
                        theme(axis.text.x = element_text(angle = 90, hjust = 1),
                              legend.position="none") 

pdf(file.path(figDir, "SupFig2_BoxplotFourBioMarkers.pdf"),
    width=6, height=6)
plot(gg)
dev.off()

#'
#' Supplemental Figure 3: PCA plot
#' 
tmp <- paste0("Group", sanitized.dds$RNA.Cluster, "_", sanitized.dds$pheno_type)
fac <- factor(tmp, levels=c("Group1_Control", "Group1_FSHD",
                       "Group2_FSHD", "Group3_FSHD", "Group4_FSHD"))
sanitized.rlg$fac <- sanitized.dds$fac <- fac
data <- plotPCA(sanitized.rlg, intgroup=c("pheno_type"), returnData=TRUE)
data$RNA.Group <- fac
percentVar <- round(100 * attr(data, "percentVar"))
gg <- ggplot(data, aes(x=PC1, y=PC2, color=RNA.Group, label=rownames(data))) +
    geom_point() +
        theme(legend.position="bottom") +
          geom_text(hjust=0.7, vjust=1.5, show.legend=FALSE) + 
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance"))
pdf(file.path(figDir, "SupFig3_PCAPlot.pdf"))
plot(gg)
dev.off()
                      
#'
#'tools: make panel of genes
#'
make_panel_of_genes <- function(gene_name, 
                                sanitized.dds) {
  i <- which(rowData(sanitized.dds)$gene_name %in% gene_name)
  sub <- sanitized.dds[i]
  data <- log10(assays(sub)[["FPKM"]]+0.5)
  data <- melt(data)
  colnames(data) <- c("id", "sample_name", "value")
  data$group <- rep(sub$group, each=nrow(sub))
  symbol <- rowData(sub)$gene_name
  data$symbol <- rep(symbol, ncol(sub))
  data$status <- "normal"
  data$status[data$sample_name %in% c("01-0037", "32-0016")] <- "outlier"
  gg <- ggplot(data, aes(x=group, y=value)) +
    geom_boxplot(outlier.colour=NA) +
    geom_jitter(width=0.2, size=0.6, alpha=0.6, aes(colour=factor(status))) +
    theme_bw() + 
    facet_wrap(~symbol, ncol=3, scales="free_y") +
    ylab(expression(log[10](RPKM))) + xlab("sample group") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position="none") +
    scale_colour_manual(values=c("black", "red"))
}
#'
#' make panel of nine (1, 2) for stephen (8/10/18)
#' 
gene_name <- c("COMP",   "PLA2G2A", "SAA1",
               "SFRP2",  "CDKN1A",  "COL1A1",
               "COL1A2", "COL3A1",  "CIDEC")
all(gene_name %in% rowData(sanitized.dds)$gene_name)
gg <- make_panel_of_genes(gene_name, sanitized.dds)
pdf(file.path(figDir, "SupFig_panel_of_nine_1.pdf"),
    height=7, width=7)
gg
dev.off()
#' another 9 gene panel
gene_name <- c("TIMP1", "C1QA", "C1QB", 
               "C2", "C3", "C7",
               "CCL18", "CTLA4", "CD44")
all(gene_name %in% rowData(sanitized.dds)$gene_name)
gg <- make_panel_of_genes(gene_name, sanitized.dds)
pdf(file.path(figDir, "SupFig_panel_of_nine_2.pdf"),
    height=7, width=7)
gg
dev.off()

#'
#' Supplemental Figure 5: DE in cluster 4 and related
#' to Immune, cell death, extraceller matrix and 
#' inflammatory
#'
gene_name <- c("C3", "CD14", "COL3A1", "COL1A2", "COL1A1", "C1QA", "TIMP1",
                "CDKN1A", "CD44")
all(gene_name %in% rowData(sanitized.dds)$gene_name)
gg <- make_panel_of_genes(gene_name, sanitized.dds)
pdf("SupFig5_InflamImmuneGenes.pdf", height=6.7, width=7)
plot(gg)
dev.off()

#'
#' Supplemental Figure 4: genes that show high expression in two group 1 sample
#' 01-0037 and 32-0016
#'
gene_name <- c("SAA1", "PLA2G2A", "C7", "COMP", "SFRP2")
all(gene_name %in% rowData(sanitized.dds)$gene_name)
gg <- make_panel_of_genes(gene_name, sanitized.dds)
pdf("SupFig4_SigGenesInGroup1.pdf", width=7, height=5)
plot(gg)
dev.off()