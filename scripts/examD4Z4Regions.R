#' 
#'  ml R/3.4.0-foss-2016b-fh1
#'
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
bamDir <- file.path(pkgDir, "tophat_D4Z4", "bam")
#'make sure all the tophat jobs have completed
tmp <- list.files(bamDir, pattern=".txt")
length(tmp)
names(tmp) <- sapply(strsplit(tmp, ".", fixed=TRUE), "[[", 1)

bamFiles <- list.files(bamDir, pattern="\\.bam$",
                       full.name=TRUE)
si <- data.frame(sample_name = sub(".bam", "", basename(bamFiles)),
                 file_bam = bamFiles,
                 stringsAsFactors=FALSE)
rownames(si) <- si$sample_name
#si <- SGSeq::getBamInfo(si)

#'
#' count exon1b exon2b and exon3
#'
library(rtracklayer)
library(GenomicAlignments)

features <- import("/fh/fast/tapscott_s/CompBio/D4Z4/sequence/D4Z4_chr4Extended.BED")
exons <- features[3:5]
exons <- c(exons, exons[1])
end(exons[4]) <- end(exons[2])
exons$name[4] <- "exon1b-2b"
names(exons) <- c("exon1", "exon2", "exon3", "exon1_intron_2")
se_D4Z4 <- summarizeOverlaps(reads=si$file_bam,
                             features=exons,
                             mode="Union",
                             inter.feature=FALSE,
                             ignore.strand=TRUE)
colData(se_D4Z4) <- as(si, "DataFrame")
colnames(se_D4Z4) <- si$sample_name

#'
#' find cluster arragement
#' (sanitized.dds must have the information)
#'
library(DESeq2)
load(file.path(pkgDir, "data", "sanitized.dds.rda"))
se <- se_D4Z4[, colnames(sanitized.dds)]
se$pheno_type <- sanitized.dds$pheno_type
se$RNA.Cluster <- sanitized.dds$RNA.Cluster
se$sizeFactor <- sanitized.dds$sizeFactor
se_D4Z4 <- se
save(se_D4Z4, file=file.path(pkgDir, "data", "se_D4Z4.rda"))

#'
#' Report and Viz
#' (1) count and normalized counts
#' 
se <- get(load(file.path(pkgDir, "data", "se_D4Z4.rda")))
sf <- matrix(rep(se$sizeFactor, each=4), nrow=4)
cnt <- assays(se)[[1]]
cnt_norm <- cnt * sf
cnt <- as.data.frame(cnt)
cnt_norm <- as.data.frame(cnt_norm)
cnt[5, ] <- colSums(cnt[1:3, ])
cnt_norm[5, ] <- colSums(cnt_norm[1:3, ])
rownames(cnt)[5] <- rownames(cnt_norm)[5] <- "total_exon1_exon2_exon3"
rownames(cnt_norm) <- paste0("norm_", rownames(cnt_norm))

#'
#' Viz: boxplot with dots
#'
library(ggplot2)
si <- as.data.frame(colData(se))
si <- cbind(si, as.data.frame(t(cnt)))
si <- cbind(si, as.data.frame(t(cnt_norm)))

write.csv(si, file=file.path(pkgDir, "stats", "D4Z4_re_counts.csv"))

mycolor=c("green", "blue")
fill = mycolor[si$pheno_type]

g1 <- ggplot(si, aes(x=RNA.Cluster, y=exon1)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon1",
                             y="Number of reads aligned to exon1")
g2 <- ggplot(si, aes(x=RNA.Cluster, y=exon2)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon2",
                             y="Number of reads aligned to exon2")
g3 <- ggplot(si, aes(x=RNA.Cluster, y=exon3)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon3",
                             y="Number of reads aligned to exon3")
g4 <- ggplot(si, aes(x=RNA.Cluster, y=exon1_intron_2)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon1-intron-exon2",
                             y="Number of reads aligned to exon1, intron and exon2")

g5 <- ggplot(si, aes(x=RNA.Cluster, y=total_exon1_exon2_exon3)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon1 to 3",
                             y="Number of reads aligned to exon1 to 3 \n (no intron)")


si[si$pheno_type=="Control", ] ## control 2538 has 11 counts on exon1 to 3

pdf(file.path(pkgDir, "figures", "D4Z4_counts.pdf"))
plot(g1)
plot(g2)
plot(g3)
plot(g4)
plot(g5)
dev.off()

#' norm counts
g1 <- ggplot(si, aes(x=RNA.Cluster, y=norm_exon1)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon1",
                             y="Number of reads aligned to exon1")
g2 <- ggplot(si, aes(x=RNA.Cluster, y=norm_exon2)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon2",
                             y="Number of reads aligned to exon2")
g3 <- ggplot(si, aes(x=RNA.Cluster, y=norm_exon3)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon3",
                             y="Number of reads aligned to exon3")
g4 <- ggplot(si, aes(x=RNA.Cluster, y=norm_exon1_intron_2)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon1-intron-exon2",
                             y="Number of reads aligned to exon1, intron and exon2")

g5 <- ggplot(si, aes(x=RNA.Cluster, y=norm_total_exon1_exon2_exon3)) +
    geom_boxplot(outlier.size=0) +
        geom_jitter(shape=16, position=position_jitter(width=0.2, height=0.05),
                    color=fill) +
                        labs(title="exon1 to 3",
                             y="Number of reads aligned to exon1 to 3 \n (no intron)")


si[si$pheno_type=="Control", ] ## control 2538 has 11 counts on exon1 to 3

pdf(file.path(pkgDir, "figures", "D4Z4_norm_counts.pdf"))
plot(g1)
plot(g2)
plot(g3)
plot(g4)
plot(g5)
dev.off()

