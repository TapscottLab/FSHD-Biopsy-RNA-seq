#' heatmap: column annotation by input argument "fac". row by selected genes
plotHeatmap <- function(biomarker.rlg, fac, print.rownames=TRUE, main=NULL,
                    filename=NULL, ...) {
    require(pheatmap)
    data <- assay(biomarker.rlg)
    rownames(data) <- mcols(rowRanges(biomarker.rlg))$gene_name
    #' make the annotation_col more flexible
    annotation_col <- do.call(cbind, lapply(fac, as.data.frame))
    colnames(annotation_col) <- names(fac)
    rownames(annotation_col) <- colnames(biomarker.rlg)
    #annotation_col <- data.frame(Group=fac)
    #rownames(annotation_col) <- colnames(biomarker.rlg)
    
    pheatmap(data,
             annotation_col=annotation_col,
             main=main,  fontsize_row=6,
             filename=filename, ...)
}

getInformativeDE <- function(biomarkers, top=5) {
    #' find most informative DEs with highest SD across samples
    sub <- sanitized.dds[biomarkers]
    sds <- rowSds(assays(sub)[["FPKM"]])
    sorted <- biomarkers[order(sds, decreasing=TRUE)]
    top_markers <- sorted[1:top]
    #print(rowData(sub[top_markers])$gene_name)
    return(top_markers)
}

.cleanGeneID <- function(genecode.id) {
    sapply(strsplit(genecode.id, ".", fixed=TRUE),
                           "[[", 1)
}

#' parallel: by normalized counts or rlg?
plotParallel <- function(biomarker.dds, fac, main="Some Biomarkers",
                         key=TRUE, key.columns=3, FPKM=FALSE) {
    require(lattice)
    data <- sapply(levels(fac), function(x) {
        if (!FPKM) 
            data <- counts(biomarker.dds, normalized=TRUE)[,
                                              fac==x]
        
        if (FPKM) 
            data <- assays(biomarker.dds)[["FPKM"]][, fac==x]
                                           
        rowMedians(data)
        #rowMeans(data)
    })
    data <- as.data.frame(data)

    ylabel <- ifelse(FPKM, "medians (FPKM)", "medians (normalized counts)")
    ## annotation
    data$gene_name <- rowData(biomarker.dds)$gene_name
    if (key) {
        parallelplot( ~ data[1:length(levels(fac))], groups=gene_name,
                       panel=function(...) {
                          panel.parallel(..., common.scale=TRUE)
                       },
                       data=data, lty=1, #col="grey",
                       ylab=ylabel, alpha=0.5,
                       horizontal.axis=FALSE,
                       main=main,
                     auto.key=list(space="top", columns=key.columns))
    }
    else {
        parallelplot( ~ data[1:length(levels(fac))], 
                       panel=function(...) {
                          panel.parallel(..., common.scale=TRUE)
                       },
                       data=data, lty=1, col="grey",
                       ylab=ylabel, alpha=0.5,
                       horizontal.axis=FALSE,
                     main=main)
    }
}


perGeneBoxplotCount <- function(id, dds) {
    require(ggplot2)
    #' id: only one gene as the time
    #' make boxplot and stripplot
    cnt <- data.frame(count=counts(dds, normalized=TRUE)[id, ],
                      group=paste0("Group", dds$RNA.Cluster, "_", dds$pheno_type))
    gene_name <- rowData(dds[id])$gene_name
    cnt$group <- factor(cnt$group, levels=c("Group1_Control", "Group1_FSHD",
                       "Group2_FSHD", "Group3_FSHD", "Group4_FSHD"))
    ggplot(cnt, aes(x=group, y=count)) + geom_boxplot(outlier.colour=NA) +
        geom_jitter(width=0.2, size=0.8) + labs(title=gene_name, y="normalized counts") +
            theme(axis.text.x = element_text(angle=20, size=8),
                  axis.title.x= element_blank(),
                  plot.title = element_text(size=12)) 
            
}
