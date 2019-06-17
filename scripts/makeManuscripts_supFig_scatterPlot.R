#' ml R/3.4.3-foss-2016b-fh1
#' supplemental Figures 0: scatter plot of RNA score and DUX4 target expression
library(xlsx)
library(aod)
library(gridExtra)
library(ggplot2) # path score axis: 0, 2, ... 10
library(lattice)
library(pheatmap)

pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
figDir <- file.path(pkgDir, "manuscript", "figures")

rank4 <- read.xlsx(file=file.path(pkgDir, "stats",
                       "RankSummary_FourBiomarkers_AllFirstVisit.xlsx"),
                   sheetIndex=1)

#' fix cluster label
RNA.Group <- rep(1, nrow(rank4))
RNA.Group[rank4$cluster==2] <- 1
RNA.Group[20:24] <- 2
RNA.Group[rank4$cluster==1] <- 3
RNA.Group[rank4$cluster %in% c(3,4)] <- 4

rank4$RNA.Group <- factor(RNA.Group, levels=c(1, 2, 3, 4))

#'
#' rename sampleID using renameSampleID() later
#' 
source(file.path(pkgDir, "scripts", "renameSample.R"))
colnames(rank4)[1] <- "sampleID" #' call renameSampleID() later
rank4$sampleID <- as.character(rank4$sampleID)
rank4$original_name <- rank4$sampleID
rank4$sampleID <- renameSample(rank4$sampleID)

#'
#' append * to FSHD samples
#' 
rank4$sampleID <- as.character(rank4$sampleID)
#' append * to FSHD samples
i <- which(rank4$pheno_type == "FSHD")
rank4$sampleID[i] <- paste0(rank4$sampleID[i], "*")
#rank4$sampleID <- paste0(substr(rank4$pheno_type, 1, 1), "_", rank4$sampleID)
rank4$sampleID <- factor(rank4$sampleID, levels=rank4$sampleID[order(rank4$log.sum)])

scatterPlot = ggplot(rank4, aes(x=sampleID, y=log.sum, color=RNA.Group)) +
        geom_point(aes(size=calls+1)) + guides(size=FALSE) +
            theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
                  legend.justification=c(0,1), legend.position=c(0, 1)) +
                      labs(y="sum of relative expression \n (rlog - average control)",
                           x="sample name")

#' scaled the rlog by the average control values
library(reshape2)
rlg <- rank4[, 13:16]
colnames(rlg) <- sub("_rlog", "", colnames(rlg))
rownames(rlg) <- as.character(rank4$sampleID)
rlg <- t(rlg)
cont.avg <- rowMeans(rlg[, rank4$pheno_type=="Control"])
relative.rlg <- rlg - matrix(cont.avg, nrow=length(cont.avg), ncol=ncol(rlg))
relative.rlg <- relative.rlg[c(1, 2, 4, 3), ]
melted <- melt(relative.rlg)

tilePlot <- ggplot(melted, aes(y=Var1, x=Var2, fill=value)) +
    geom_tile() + theme_minimal() +
        theme(legend.direction="horizontal",
              legend.position=c(0.3, 0.6), legend.justification=c(1, 0),
              #legend.title = element_blank(),
              axis.text.x=element_blank(), axis.title.x=element_blank(),
              axis.title=element_blank(), axis.text.y=element_text(size=8)) +
                  scale_fill_gradient(low="white", high="steelblue", name="relative expression") +
                      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5)) 

        
pdf(file.path(figDir, "SupFig0_Scatter_FourBiomarkers.pdf"))
grid.arrange(tilePlot, scatterPlot, nrow=2, ncol=1, heights=c(1,3))
dev.off()
