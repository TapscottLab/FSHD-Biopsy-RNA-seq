library(DESeq2)
library(lattice)
library(xlsx)

pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
figDir <- file.path(pkgDir, "figures")

#' clinical data
clinic <- getClinicalsFromMasterSheet(file=file.path(pkgDir, "extdata",
                                          "Tawil-Master-4-13-17.xlsx"))
rownames(clinic) <- sub("01-", "", rownames(clinic), fixed=TRUE)

clinic <- clinic[, c("Path.Score", "Inflamm", "Active",
                     "I..Variability.in.Fiber.Size",
                     "IV..Interstitial.Fibrosis",  "STIR",
                     "T1", "RNA.Score", "RNA.Cluster")]
idx <- which(rownames(clinic)=="32-0008")
clinic$RNA.Cluster[idx] <- 3
#' clean up samples without clinical data or RNA data
idx <- which(rownames(clinic) %in% c("2393", "32-0011", "32-0016", "2320"))
clinic <- clinic[-idx, ]

#' correlation: fix fibrosis, STIR and T1
what_class <- sapply(clinic, class)
clinic$IV..Interstitial.Fibrosis <-
    factor(clinic$IV..Interstitial.Fibrosis,
           levels=c("Normal", "Mild", "Moderate", "Severe"))
clinic$I..Variability.in.Fiber.Size <-
    factor(clinic$I..Variability.in.Fiber.Size,
           levels=c("Normal", "Mild", "Moderate", "Severe"))
clinic$STIR <- factor(clinic$STIR,
                      levels=c("0", "<0.5", "0.5", "1", "2"))
clinic$T1 <- factor(clinic$T1,
                    levels=c("0", "1", "2", "2 to 3", "3"))
isFactor <- sapply(clinic, is.factor)
clinic$IV..Interstitial.Fibrosis <- as.numeric(clinic$IV..Interstitial.Fibrosis)
clinic$I..Variability.in.Fiber.Size <- as.numeric(clinic$I..Variability.in.Fiber.Size)
clinic$STIR <- as.numeric(clinic$STIR)
clinic$T1 <- as.numeric(clinic$T1)

#' correlation
isFactor <- sapply(clinic, is.factor)
cor.clinic <- cor(clinic)
ord <- order.dendrogram(as.dendrogram(hclust(dist(cor.clinic))))

panel.corrgram.2 <-
    function(x, y, z, subscripts, at=pretty(z), scale=0.8, ...) {
        require("grid", quietly=TRUE)
        x <- as.numeric(x)[subscripts]
        y <- as.numeric(y)[subscripts]
        z <- as.numeric(z)[subscripts]
        zcal <- level.colors(z, at=at, ...)
        for (i in seq(along=z)) {
            lims <- range(0, z[i])
            tval <- 2 * base::pi *
                seq(from=lims[1], to=lims[2], by=0.01)
            grid.polygon(x=x[i] + 0.5 * scale * c(0, sin(tval)),
                         y=y[i] + 0.5 * scale * c(0, cos(tval)),
                         default.units="native",
                         gp=gpar(fill=zcal[i]))
            grid.circle(x=x[i], y=y[i], r=0.5*scale,
                        default.units="native")
        }
    }

corgram <- levelplot(cor.clinic[ord, ord], xlab=NULL, ylab=NULL,
          at=do.breaks(c(-.01, 1.01), 101),
          panel=panel.corrgram.2,
          scales=list(x=list(rot=90)),
          corlorkey=list(space="top"),
          col.regions=colorRampPalette(c("white", "blue")))

pdf(file.path(pkgDir, "figures", "Corrgram_ClinicalScores.pdf"))
print(corgram)
dev.off()
