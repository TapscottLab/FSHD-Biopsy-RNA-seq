#'
#'   ml R/3.4.3-foss-2016b-fh2
#'
library(xlsx)
library(ggplot2)
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
figDir <- file.path(pkgDir, "manuscript", "figures")

clinics <- read.xlsx(file.path(pkgDir, "extdata", "WellstoneProj2.Data.Alt.UR.ID.xlsx"),
                     sheetIndex=1, startRow=2)

pdf(file.path(figDir, "PathScoreVSRNAScore_scatterplot.pdf"),
              height=3, width=5)
xyplot(Path.Score ~ Score, data=clinics, pch=20,
          panel=function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.lmline(x, y, col="grey")
       },
       xlab="RNA-seq Score", ylab="Muscle Pathology Score")
dev.off()

filename <- file.path(figDir, "PathScore_vs_RNAScore.pdf")
gg <- ggplot(clinics, aes(x=Score, y=Path.Score)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    theme_minimal() +
    labs(y="Muscle Pathology Score",
           x="RNA-seq Score")
ggsave(filename=filename, gg, device="pdf", width=6, height=3)
