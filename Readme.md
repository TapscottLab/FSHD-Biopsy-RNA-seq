This package contains both datasets and scripts (R) supporting the manuscript of MRI-informed muscle biopsies correlate MRI with pathology and DUX4 target gene expression in FSHD.  The shell script `/scripts/do_tophat.sh` performed the preprocessing and `/scripts/*.R` the analysis for the manuscript. See more detail in `/scripts/ReadMe.txt`.

## Folders
List below are the major datasets and scripts that were used for the analysis. Note that the scripts for data exploratory are also included in the *scripts* folder, and so they might be a bit messy.

<pre>
|- data    
  |- se.rda: `RangedSummarizedExperiment` instance as a container containing gene counts, gene annotation and 
     metadata. 
  |- sanitized.dds.rda: sanitied/filtered `DESeqDataSet` of `se`, containing dispersion estimation and size 
     factor, used for downstream analysis  
  |- sanitized.rlg.rda: sanitized/filtered rlog-transformed dataset of sanitized.dds, used to visualize 
     sample space    
  
|- scripts      
  |- *.Rmd:  R markdown document and notebook  
  |- makeSE.R: get gene count and make a `RangedSummarizedExperiment` instance  
  |- rankScore.R: get DUX4 rank and scores
  |- makeManuscripts_scatterPlot.R: make Figure 1
</pre>

### Main analysis

makeSE.R -> rankScore.R -> makeManuscripts_scatterPlot.R -> Figure 1

## system requirement
- R (3.5): pheatmap, xlsx, ggplot2, gridExtra
- Bioconductor (3.7)
  - DESeqs (differential analysis)
  - goseq (GO analysis)
  - GenomicAlignments (gene counts)
  - org.Hs.eg.db (annotation)
  - GO.db (GO annotation)
- Optional (Customized pakcages for gene counts)
  - hg38.HomoSapiens.Gencode.v24 (standard TxDb package built from gencode GFF v24)
