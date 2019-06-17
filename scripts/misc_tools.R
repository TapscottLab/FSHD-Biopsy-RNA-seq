renameSampleName <- function(seq_name) {
    require(xlsx)
    df <- read.xlsx(file.path("/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV",
                              "extdata", "WellstoneProj2.Data.Alt.UR.ID.xlsx"),
                    sheetIndex=1, startRow=2, endRow=38, colIndex=1:2,
                    stringsAsFactors=FALSE)
    df[is.na(df$Alt.UR.ID), 1] <- df$record_id[is.na(df$Alt.UR.ID)]
    sample_name <- df$record_id
    sample_name <- sub("--", "-", sample_name)
    sample_name <- sub("01-", "", sample_name)
    alt_id <- df$Alt.UR.ID
    names(alt_id) <- sample_name
    #' add nine control samples
    cntr <- paste0("01-", "00", seq(41,49))
    names(cntr) <- c("2524", "2536", "2538", "2539", "2540",
                     "2550", "2559", "2563", "2578")
    #' 0952-1
    add_on <- "01-0022-1"
    names(add_on) <- "0952-1"
    
    #' combine
    alt_id <- c(alt_id, cntr, add_on)
    
    #' santity check
    if (all(seq_name %in% names(alt_id)))
        return(alt_id[seq_name])
    else {
        i <- which(!seq_name %in% names(alt_id))
        message("Cannot find alternative name to match ", paste(seq_name[i], collapse=","))
    }
}

getLogScale <- function(se) {
    ## simple log scale: linearly adjust the count by size factor than take (log2(count+0.5))
    sizeFactor <- se$sizeFactor
    cnt <- assays(se)[[1]]
    l <- nrow(cnt)
    mf <- matrix(rep(sizeFactor, l), nrow=l, byrow=TRUE)
    cnt <- cnt / mf
    log2(cnt+0.5)
}

getMarkerCounts <- function(markers=c("LEUTX", "PRAMEF2", "TRIM43", "KHDC1L"),
                            txdb, bamFiles, 
                            anno=NULL, txdb.type="Gencode", workers=1) {
    exons <- exonsBy(txdb, by="gene")
    
    if (txdb.type == "UCSC")
        marker_id <- mapIds(org.Hs.eg.db, markers,
                            "ENTREZID", "SYMBOL", multiVals="first")
    if(txdb.type == "Gencode") {
        marker_id <- rownames(anno)[anno$gene_name %in% markers]
    }
    
    features <- exons[marker_id]
    mparam <- MulticoreParam(workers = workers, type = "SOCK")
    se <- summarizeOverlaps(features=features,
                            reads=bamFiles,
                            mode="IntersectionStrict",
                            inter.feature=TRUE,
                            singleEnd=TRUE,
                            ignore.strand=TRUE,
                            BPPARAM=mparam)
    mcols(rowRanges(se)) <- anno[rownames(se), ]
    colnames(se) <- sub(".bam", "", colnames(se))
    se
}

patchAnnotationSE <- function(se, gene.anno) {
    mcols(rowRanges(se))$gene_name <- gene.anno[rownames(se), "gene_name"]
    mcols(rowRanges(se))$gene_type <- gene.anno[rownames(se), "gene_type"]
    mcols(rowRanges(se))$gene_status <- gene.anno[rownames(se), "gene_status"]
    se
}
 
#' DUX4-targeted markers
getDUX4Targets <- function() {
    require(xlsx)
    sig <-  read.xlsx(file=file.path( "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg19.FSHD.biopsy",
                      "inst", "extdata", "ddu251supp_table3.xls"),
                      sheetIndex=1, startRow=2)
}

getRankScore <- function(dds, which.genes, lfcThreshold=2,
                        num.cluster, main=NULL, 
                         simple.log2=FALSE) {
    library(lattice)
    library(pheatmap)
    library(gridExtra)
    #' output: (1) rank by clustering; (2) rank by expressed level
    #' transformation; (3) rank by sum

    #' subset dds
    rank <- vector("list", 3)
    names(rank) <- c("cluster", "log.sum", "calls")

    ## subset the dds, and report which are not in the list
    keep <- intersect(rownames(dds), which.genes)
    low_count_genes <- which.genes[!(which.genes %in% keep)]
    message("biomarkers not keep:")
    print(low_count_genes)
    
    x <- dds[keep, ]
    rownames(x) <- rowData(x)$gene_name
    ## get log transformation: either by simple log2 or rlog
    #if (simple.log2) {
    #    cnt <- counts(x, normalized=TRUE)
    #    lg <- log2(cnt+0.5)
    #}
    #else {
        tmp <- rlog(x, fitType="local")
        lg <- assay(tmp)
    #}
    avg <- rowMeans(lg)
    centered.lg <- lg - matrix(avg, nrow=length(avg), ncol=ncol(lg))
    cont.avg <- rowMeans(lg[, dds$pheno_type == "Control"])
    relative.lg <- lg - matrix(cont.avg, nrow=length(cont.avg), ncol=ncol(lg))
    
    ## (1) cluster, cut tree for num.cluster blanchees
    ##hc <- hclust(dist(t(relative.lg)))
    hc <- hclust(dist(t(lg)))
    rank$cluster <- cutree(hc, k=num.cluster)

    ## (2) calls: how many have relative log greater than the threshold
    ##rank$calls <- factor(colSums(relative.lg > lfcThreshold))
    rank$calls <- colSums(relative.lg > lfcThreshold)
    
    ## (3) sum of relative log expression
    rank$log.sum <- colSums(relative.lg)

    #'
    #' visualization
    #'

    ## (4) heatmap
    annotation_col <- data.frame(pheno_type=dds$pheno_type)
    ##Batch=dds$Batch)
    ##Calls=rank$calls)
    rownames(annotation_col) <- colnames(relative.lg)
    fontsize_row <- 10
    fontsize_col <- 10
    if (nrow(relative.lg) > 10) fontsize_row <- 7
    if (ncol(relative.lg) > 30) fontsize_col <- 8

   pheatmap(centered.lg, annotation_col=annotation_col, silent=TRUE,
            main=main, fontsize_row=fontsize_row,
            fontsize_col=fontsize_col,
            file=file.path(pkgDir, "figures", paste0("Heatmap_", main, ".pdf")))

    ## (5) scatter plot
    sampleID <- paste0(substr(dds$pheno_type, 1, 1), "_", colnames(dds), "_", dds$Batch)
    df <- data.frame(sampleID=sampleID,
                     Batch=dds$Batch,
                     Group=dds$pheno_type,
                     Calls=as.numeric(as.character(rank$calls)),
                     Log.sum=rank$log.sum,
                     Cluster=factor(rank$cluster))
    df$sampleID <- factor(df$sampleID, levels=df$sampleID[order(df$Log.sum)])

    scatterPlot = ggplot(df, aes(x=sampleID, y=Log.sum, color=Cluster)) +
        geom_point(aes(size=Calls+1)) + guides(size=FALSE) +
            theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
                  legend.justification=c(0,1), legend.position=c(0, 1)) +
                      labs(y="Sum of relative expression (rlog)")
    ## scatterPlot = scatterPlot + guides(colour=FALSE)
    barPlot = ggplot(df, aes(x=sampleID, y=Calls)) + geom_bar(stat="identity") +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +
            labs(title=main)
  
    pdf(file.path(pkgDir,"figures", paste0("Scatter_", main, ".pdf")))
    grid.arrange(barPlot, scatterPlot, nrow=2, ncol=1, heights=c(1,3))
    dev.off()

    cnt <- counts(x, normalized=TRUE)  
    return(list(rank=rank, rlog=lg, relative.log=relative.lg, normalized.cnt=cnt))
}

sanitizeSamples <- function(se) {
    ## 0952 from the second batch is better, 32-0002b is the second visit,
    ## and 2374 are dupliced from older cell line. Remove them from downstream analysis. 
    rm_samples <- c("0952", "32-0008-1", "32-0002b", "2374-1", "2374-2")
    se[, !colnames(se) %in% rm_samples]
}
                               
getRankSummary <- function(dds, rank, desDir, main=NULL, generate.report=TRUE,
                           sep=FALSE) {
    ## (1) Print out spreat sheet
    ## (2) boxplot of the markers  expressions

    #' 1-rank
    smry <- as.data.frame(do.call(cbind, rank$rank))
    smry$pheno_type <- dds$pheno_type
    smry$batch <- dds$Batch
    
    smry <- smry[order(smry$log.sum), ]
    smry$rank <- seq(1, nrow(smry))
    smry <- smry[, c("pheno_type", "rank", "cluster", "log.sum", "calls", "batch")]
    smry$score <- smry$log.sum - min(smry$log.sum)
    #' 2-normaized counts, 3-rlog
    rlg <- rank$rlog
    rlg <- rlg[, rownames(smry)]
    cnt <- rank$normalized.cnt
    cnt <- cnt[, rownames(smry)]

    library(xlsx)
    filename <- file.path(desDir, paste0("RankSummary_", main, ".xlsx"))

    if (generate.report) {
        if (!sep) {
            rownames(rlg) <- paste0(rownames(rlg), "_rlog")
            rownames(cnt) <- paste0(rownames(cnt), "_norm.counts")
            tmp <- cbind(smry, t(cnt), t(rlg))
            write.xlsx(tmp, file=filename, sheetName="rank scores", append=FALSE)
        }
        if (sep) {
            write.xlsx(smry, file=filename, sheetName="rank scores", append=FALSE)
            write.xlsx(t(rlg), file=filename, sheetName="rlog", append=TRUE)
            write.xlsx(t(cnt), file=filename, sheetName="normalized counts", append=TRUE)
        }
    }
    ## (2) boxplots
    return(list(smry=smry, rlg=rlg, cnt=cnt))
}



heatmapArrangeByClass <- function(rlg, id) {
    require(pheatmap)
    class <- paste0(dds$class4, "-", dds$pheno_type)
    class <- factor(class, levels=c("cold-Control", "cold-FSHD", "warm-FSHD",
                               "warmer-FSHD", "hot-FSHD"))
    rlg$class <- class
    #' arrange by control, cold, warm, warmer and hot
    data_class <- lapply(levels(rlg$class), function(x) {
        data <- assay(rlg)[id, rlg$class == x]
        data[, order(colMeans(data), decreasing=FALSE)]
    })
    data_class <- do.call(cbind, data_class)
    annotation_col <- data.frame(class=rlg[, colnames(data_class)]$class)
    rownames(annotation_col) <- colnames(data_class)
    pheatmap(data_class,
             cluster_cols=FALSE,
             annotation_col=annotation_col)
}


getClinicalsFromMasterSheet <- function(file=NULL) {
    require(xlsx)
    if (is.null(file)) {
        pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsy4"
        file <- file.path(pkgDir, "inst", "extdata", "Tawil-Master-4-13-17.xlsx")
    }

    
    df <- read.xlsx(file, sheetIndex=1, startRow=2, endRow=38,
                    colIndex=1:29, stringsAsFactors=FALSE)
    rownames(df) <- sub("--", "-", df$record_id)


    colnames(df)[c(28, 29)] <- c("RNA.Score", "RNA.Cluster")
    df
}

addClinicalsToColData <- function(dds) {
    clinic <- getClinicalsFromMasterSheet()
    rownames(clinic) <- sub("01-", "", rownames(clinic), fixed=TRUE)
    clinic <- clinic[, c(1, 19:28)]
    #' matching sanitized.dds and clinic sample names
    which_0952 <- which(rownames(clinic) == "0952")
    rownames(clinic)[which_0952] <- "0952-1"
    #' clean up the column name
    colnames(clinic) <- sapply( strsplit(colnames(clinic), "..", fixed=TRUE),
                               function(x) {
                                   x[length(x)]
                               })
                                        # ADD column by colmn 
    sample_with_clinic <- intersect(rownames(clinic), colnames(sanitized.dds))
    colData <- colData(dds)
    for (column_name in colnames(clinic)) {
        colData[, column_name] <- NA
        colData[sample_with_clinic, column_name] <-
        clinic[sample_with_clinic, column_name]
    }

    #' fix the factors: Variability.in.Fiber.Size, Extent.of.Central.Nucleation,
    #' Necrosis.Regeneration.Inflammation, Interstitial.Fibrosis, STIR and T1
    factorize_column <- c("Variability.in.Fiber.Size",
                          "Extent.of.Central.Nucleation",
                          "Necrosis.Regeneration.Inflammation",
                          "Interstitial.Fibrosis")
    levels <- c("Normal", "Mild", "Moderate", "Severe")
    for (col in factorize_column)
        colData[, col] <- factor(colData[, col], levels=levels)

    colData[, "STIR"] <- factor(colData[, "STIR"],
                                levels=c("0", "<0.5", "0.5", "1", "2"))
    colData[, "T1"] <- factor(colData[, "T1"],
                              levels=c("0", "1", "2", "2 to 3"))
    colData
}

makePCAPlotByPhenoType <- function(se, title="Title", figDir, threshold.count=5) {
    require(DESeq2)
    require(ggplot2)
    dds <- DESeqDataSet(se[rowSums(assays(se)[[1]]) > threshold.count, ],
                        design = ~ pheno_type)
    rlg <- rlog(dds)
    data <- plotPCA(rlg, intgroup=c("pheno_type"), returnData=TRUE)
    percentVar <- round(100 * attr(data, "percentVar"))
    ggplot(data, aes(x=PC1, y=PC2, color=group, label=rownames(data))) +
               geom_point() + 
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                        labs(list(title=title)) +
                            geom_text(hjust=0, vjust=0.5)

    ggsave(file=file.path(figDir, paste0(title, "_SampleDistancePCAplot.png")))
}

do_deseq <- function(se, mk, logFC.threshold=2, alpha=0.05,
                     title, statsDir) {
    require(xlsx)
    require(DESeq2)
    ## mk is a SummarizedExperiment instance containing the counts for biomarkers
    dds <- DESeqDataSet(se, design = ~ pheno_type)
    ## rlg <- rlog(dds, blind=FALSE)
    dds <- DESeq(dds)
    res <- results(dds)
    #res <- results(dds, lfcThreshold = 1)
    res <- res[order(res$padj, decreasing=FALSE), ]
    gene.anno <- mcols(rowRanges(dds))
    rownames(gene.anno) <- rownames(dds)
    res$gene_name <- gene.anno[rownames(res), "gene_name"]
    res$rank <- seq(1, nrow(res))
    res$biomarker <- rownames(res) %in%  rownames(mk)

    ## append normalized count
    cnt <- counts(dds, normalized=TRUE)
    colnames(cnt) <- paste0(dds$pheno_type, "_", colnames(cnt))
    cnt <- cnt[rownames(res), ]
    res <- cbind(res, as(as.data.frame(cnt), "DataFrame"))
    ## 54 biomarkers and four biomarkers
    i <- intersect(rownames(res), rownames(mk))
    tmp <- res[i, ]
    de_idx <- which(tmp$padj < alpha & abs(tmp$log2FoldChange) >= logFC.threshold)
    if (length(de_idx) > 0) {
        message("Number of biomarkers that are DE (padj < 0.05, LFC>2): ", length(de_idx))
        print(as.data.frame(tmp[de_idx, ]))
    }
    id <- c("ENSG00000120952.4", "ENSG00000144015.4", "ENSG00000256980.4",
                "ENSG00000213921.6")
    tmp.fourmk <- res[id, ]
    #' viz? parallel plot? boxplot plot for four biomarkers and zscan4?
    
    #' report
    file <- file.path(statsDir, paste0("DESeqResults_", title, ".xlsx"))
    sub.res <- subset(res, res$padj < 0.05 & abs(res$log2FoldChange) >=logFC.threshold)
    sheetName <- "DE padj<0.05, LFC >=2"
    write.xlsx(as.data.frame(sub.res), file=file, sheetName=sheetName, append=FALSE)
    write.xlsx(as.data.frame(tmp.fourmk), file=file, sheetName="Four biomarkers statistics",
               append=TRUE)
    write.xlsx(as.data.frame(tmp), file=file, sheetName="All biomarkers statistics",
               append=TRUE)
    
    DESeqResults(res)
}

summaryDESeqResults <- function(res, padj.threshold=0.05, lfc.threshold=2) {
    tmp=subset(res, res$padj < padj.threshold & abs(res$log2FoldChange) >= lfc.threshold)
    message("Totle DE:", nrow(tmp))
    tmp
}

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

do_goseq <- function(dds, idx) {
    sub <- dds[, idx]
    design(sub) <- ~ pheno_type
    dds <- DESeq(sub)
    res <- results(dds, alpha=0.05, lfcThreshold=2)
    res$gene_name <- rowData(dds)$gene_name
    up_sig <- subset(res, padj < 0.05 & log2FoldChange > 0)

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
}
