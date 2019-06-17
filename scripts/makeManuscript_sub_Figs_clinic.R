#' module load R/3.4.0-foss-2016b-fh1
#' Get clinical data
#' 
library(xlsx)
library(aod)
library(gridExtra)
library(ggplot2) # path score axis: 0, 2, ... 10
library(lattice)
library(pheatmap)

pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
figDir <- file.path(pkgDir, "manuscript", "figures")

file <- file.path(pkgDir, "inst", "extdata", "Tawil-Master-4-13-17.xlsx")
df <- read.xlsx(file, sheetIndex=1, startRow=2, endRow=38, colIndex=1:29, stringsAsFactors=FALSE)
rownames(df) <- sub("--", "-", df$record_id)
#df <- read.xlsx(file, sheetIndex=8, startRow=2, colIndex=1:11)
colnames(df)[c(28, 29)] <- c("RNA.Score", "RNA.Cluster")
df$RNA.Cluster[df$RNA.Cluster==5] <- 4
df$RNA.Cluster <- factor(df$RNA.Cluster, levels=c(1,2,3,4))
df$IV..Interstitial.Fibrosis <- factor(df$IV..Interstitial.Fibrosis,
                                       levels=c("Normal", "Mild", "Moderate", "Severe"))
df$STIR <- factor(as.character(df$STIR), levels=c("0", "<0.5", "0.5", "1", "2"))
df$STIR.num <- as.character(df$STIR)
df$STIR.num[df$STIR.num=="<0.5"] <- "0.25"
df$STIR.num <- as.numeric(df$STIR.num)

#' new MRI data (water and fat fraction) -> add to df
#' there is no 2483 (-> 2583), no 2359 (-> 2358)
new <- read.xlsx(file=file.path(pkgDir, "inst", "extdata", "dixon_biopsy.xlsx"), sheetIndex=1)
new$SUBID <- as.character(new$SUBID)
new$SUBID[17] <- "2583"
new[17:18, ] <- new[18:17, ]
df$Water.Fraction <- new$Water.Fraction
df$Fat.Fraction <- new$Fat.Fraction

#'
#' ggplot default color (hue)
#'
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(4)

library(scales)
cols <- hue_pal()(4)

#'
#' RNA clusters (scores) and pathlogy Score
#'
pdf(file=file.path(figDir, "RNACluster_vs_PathScore.pdf"),
    width=5, height=5)
boxplot(Path.Score ~ RNA.Cluster, data=df,
        outline=FALSE,
        xlab="RNA Cluster",
        ylab="Path Score")
stripchart(Path.Score ~ RNA.Cluster, vertical=TRUE, data=df, method="jitter",
           jitter=0.2, add=TRUE, pch=20, col="grey")
dev.off()

#'
#' RNA scores and pathlogy score
#'
tmp <- df[!is.na(df$RNA.Cluster), ]
pdf(file.path(figDir, "PathScoreVSRNAScore_scatterplot.pdf"),
              height=5, width=5)
xyplot(Path.Score ~ RNA.Score, group=RNA.Cluster, data=df, pch=20, col=cols,
       key=list(corner=c(0.94, 0), text=list(levels(df$RNA.Cluster)), cex=0.8,
           title="RNA Cluster", cex.title=0.9, border=FALSE, size=10,
               points=list(col=cols, pch=20)),
       panel=function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.lmline(x, y, col="grey")
       },
       xlab="RNA Score", ylab="Path Score")
dev.off()

#'
#' RNA clusters (or score) vs MRI-T1, STIR, and Fat Fraction 
#'

pdf(file=file.path(figDir, "MRI_vs_RNAScore.pdf"),
    width=11, height=8)
g1 <- ggplot(tmp, aes(x=Fat.Fraction, y=RNA.Score, color=RNA.Cluster)) +
    geom_point() +
            labs(y="RNA Score", x="Fat Fraction") +
            theme(legend.position = c(1, 0), legend.justification = c(1, 0),
                  legend.background = element_rect(colour = "black"),
                  legend.title=element_text(size=9), 
                  panel.background = element_rect(fill = "white", colour = "black"))

g2 <- ggplot(df, aes(x=T1, y=RNA.Score)) +
    geom_point(position=position_jitter(width=0.2), aes(color=df$RNA.Cluster))+
      geom_boxplot(outlier.colour=NA, fill=NA) +
        theme(panel.background = element_rect(fill = "white", colour = "black")) +
            labs(x="T1", Y="RNA Score") + theme(legend.position="none")
g3 <- ggplot(df, aes(x=STIR, y=RNA.Score)) +
    geom_point(position=position_jitter(width=0.2), aes(color=df$RNA.Cluster))+
      geom_boxplot(outlier.colour=NA, fill=NA) +
        theme(panel.background = element_rect(fill = "white", colour = "black")) +
          labs(x="STIR", y="RNA Score") + theme(legend.position="none")
grid.arrange(g1, g2, g3, nrow=2)
dev.off()

#'
#' PATH score vs Fat Fraction, STIR and RNA Score
#' 
pdf(file.path(pkgDir, "inst", "figures", "FatFraction_vs_PathScore.pdf"))
xyplot(Path.Score ~ Fat.Fraction, df, pch=20, col="black")
dev.off()
pdf(file.path(pkgDir, "inst", "figures", "FatFraction_vs_RNAScore.pdf"))
tmp <- df[!is.na(df$RNA.Score), ]
ggplot(tmp, aes(y=Fat.Fraction, x=RNA.Score, color=RNA.Cluster)) +
    geom_point() + theme_bw() +
        labs(y="Fat Fraction", x="RNA Score") +
            theme(legend.position = c(1, 0), legend.justification = c(1, 0),
                  legend.background = element_rect(colour = "black"))

tmp <- df[!is.na(df$RNA.Cluster), ]
ggplot(tmp, aes(y=Path.Score, x=RNA.Score, color=RNA.Cluster)) +
    geom_point() + theme_bw() +
        labs(y="Path Score", x="RNA Score") +
            theme(legend.position = c(1, 0), legend.justification = c(1, 0),
                  legend.background = element_rect(colour = "black"))

dev.off()

xyplot(Path.Score ~ Fat.Fraction, df, pch=20, col="black")
g1 <- ggplot(df, aes(y=Path.Score , x=Fat.Fraction)) +
    geom_point() +
        #geom_smooth(method="lm", color="steelblue") +
        theme_bw() +
            labs(y="Path Score", x="Fat Fraction")

g2 <- ggplot(df, aes(y=Path.Score , x=STIR)) +
    geom_boxplot() +
        geom_point(color="grey") +
            theme_bw()
            geom_jitter(position=position_jitter(0.2))
            labs(y="Path Score", x="STIR")

tmp <- df[!is.na(df$RNA.Cluster), ]
g3 <- ggplot(tmp, aes(y=Path.Score, x=RNA.Score, color=RNA.Cluster)) +
    geom_point() + theme_bw() +
        labs(y="Path Score", x="RNA Score") +
            theme(legend.position = c(1, 0), legend.justification = c(1, 0),
                  legend.background = element_rect(colour = "black"))

## another example
boxplot(Path.Score ~ STIR, data=df,
        xlab="STIR",
        ylab="Path Score")
stripchart(Path.Score ~ STIR, vertical=TRUE, data=df, method="jitter",
           jitter=0.2, add=TRUE, pch=20, col="grey")

stripplot(STIR ~ Path.Score, data=df, jitter=TRUE)

#'
#' Fat Fraction vs Path subscores
#'
mypanel <- function(x, y, ...) {
    panel.bwplot(x, y, ..., do.out=FALSE, pch="|")
    medy <- by(y, list(x), median)
    xx <- sort(unique(as.numeric(x)))
    panel.segments(xx-.25, medy ,xx+.25, medy, lwd=2.5)
    panel.stripplot(x, y, ..., jitter.data=TRUE, pch=20, col="grey")
}
mypars <- list(box.rectangle = list(col = "black"),
               box.umbrella = list(col="black"))
#' inflamm
df$Inflamm.factor <- factor(df$Inflamm, levels=c(0, 1, 2, 3))
g1 <- bwplot(Fat.Fraction ~ Inflamm.factor, data=df,
       panel=mypanel,
       par.settings=mypars,
             xlab="Subscore - Inflammation", ylab="Fat Fraction")

#' Active
df$Active.factor <- factor(df$Active, levels=c(0, 1, 2, 3))
g2 <- bwplot(Fat.Fraction ~ Active.factor, data=df,
       panel=mypanel,
       par.settings=mypars,
             xlab="Subscore - Active", ylab="Fat Fraction")

#' Variability in Fiber Size
df$I..Variability.in.Fiber.Size <-
    factor(df$I..Variability.in.Fiber.Size,
           levels= c("Normal",  "Mild",  "Moderate",  "Severe"))
g3 <- bwplot(Fat.Fraction ~ I..Variability.in.Fiber.Size , data=df,
       panel=mypanel,
       par.settings=mypars,
             xlab="Subscore - Variability in Fiber Size ",
             ylab="Fat Fraction")

#' Central Nucleation
df$II..Extent.of.Central.Nucleation <-
    factor(df$II..Extent.of.Central.Nucleation,
           levels=c("Normal", "Mild", "Moderate"))
g4 <-  bwplot(Fat.Fraction ~ II..Extent.of.Central.Nucleation , data=df,
       panel=mypanel,
       par.settings=mypars,
             xlab="Subscore - Extent of Central Nucleation ",
             ylab="Fat Fraction")
#' Necrosis/Regenration/Inflammation
df$III..Necrosis.Regeneration.Inflammation <-
    factor(df$III..Necrosis.Regeneration.Inflammation,
           levels= c("Normal",  "Mild",  "Moderate",  "Severe"))
g5 <-  bwplot(Fat.Fraction ~ III..Necrosis.Regeneration.Inflammation , data=df,
       panel=mypanel,
       par.settings=mypars,
             xlab="Subscore - Necrosis Regeneration Inflammation ",
             ylab="Fat Fraction")
#' Interstitial Fribosis
df$IV..Interstitial.Fibrosis <-
    factor(df$IV..Interstitial.Fibrosis,
           levels=c("Normal",  "Mild",  "Moderate",  "Severe"))
g6 <-  bwplot(Fat.Fraction ~ IV..Interstitial.Fibrosis, data=df,
       panel=mypanel,
       par.settings=mypars,
             xlab="Subscore - InterstitialFibrosis",
             ylab="Fat Fraction")
pdf(file.path(pkgDir, "inst", "figures", "FatFraction_vs_PathSubScore.pdf"))
g1
g2
g3
g4
g5
g6
dev.off()

#'
#' scatter plot of RNA score with samples
#'
