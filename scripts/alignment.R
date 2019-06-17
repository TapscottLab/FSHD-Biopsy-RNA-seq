#'  ml R/3.4.0-foss-2016b-fh1
#'
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
trimmedDir <- file.path(pkgDir, "trimmed")
ngsDir <- file.path("/shared/ngs/illumina/lsnider",
                    c("151002_SN367_0557_AHC5HLBCXX",
                      "160208_SN367_0613_BHJFJ3BCXX",
                      "160908_D00300_0314_AC9NKGANXX",
                      "170119_D00300_0372_BH77YFBCXY"))
sampleDir <- list.files(file.path(ngsDir, "Unaligned/Project_lsnider"),
                        full.name=TRUE)
#' 
sampleName <- sub("Sample_", "", basename(sampleDir))
sampleName[c(8, 14, 44)] <- paste0(sampleName[c(8, 14, 44)], "-1")
sampleName[22] <- paste0(sampleName[22], "-2")

#' move TruSeq3-SE.fa to trimmedDir or customize a fa?

#' sbatch work
for (i in 2:length(sampleDir)) {
    cmt <- sprintf("sbatch -n6 ./do_tophat.sh %s %s %s",
                   sampleDir[i], sampleName[i], pkgDir)
    system(cmt) #53579191
}

#' do 32-0008 from II again and trim the deteriorated tail
i <- 12
cmt <- sprintf("sbatch -n6 ./do_tophat_trimTail.sh %s %s %s",
               sampleDir[i], sampleName[i], pkgDir)
system(cmt) #53579239

#' make sure all the tophat jobs have completed
tophatDir <- file.path(pkgDir, "tophat", "bam")
tmp <- list.files(tophatDir, pattern=".txt")
length(tmp)

###################
#'
#' testing
#'
#' ################
library(ShortRead)
f1="~/tapscott/RNA-Seq/hg38.FSHD.biopsyI-IV/trimmed/2374.fastq"
trimFile="~/tapscott/RNA-Seq/hg38.FSHD.biopsyI-IV/trimmed/2374.trim.fastq"
files <- list.files(sampleDir[1], pattern=".gz", full.name=TRUE)

#'
#' ok, something is wrong with filtering process
#' 
r1=readFastq(files[12]) #2374_ATCACG_L002_R1_001.fastq.gz; filtered ones working
f1.raw <- file.path(ngsDir[1], "Unaligned/Project_lsnider/Sample_2374",
                    "2374_ATCACG_L002_R1_001.fastq.gz")
r1.raw=readFastq(fq.raw)
# problematic file: 2374.fastq (filtered fastq) error: at line 112932427, expected a line starting with '+'

tmp="2374_ATCACG_L002_R1_001.fastq.gz" lines: 126146


#' let's just filter the 2374 samples
pkgDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg38.FSHD.biopsyI-IV"
library(ShortRead)
filteredDir <- file.path(pkgDir, "filtered")
f2 <- file.path(filteredDir, "2374_ATCACG_L002_R1_001.fastq.gz") # no problem
r2=readFastq(f2)
f3 <- file.path(filteredDir, "2374.fastq")
r3=readFastq(f3)
f4 <- file.path(filteredDir, "2374.cutadapt.fastq")
r4=readFastq(f4)
