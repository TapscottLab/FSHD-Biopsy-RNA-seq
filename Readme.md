This package contains both datasets and scripts (R) supporting the manuscript of MRI-informed muscle biopsies correlate MRI with pathlogy and DUX4 target gene expression in FSHD.  The shell script `/scripts/do_tophat.sh` performed the preprocessing and `/scripts/*.R` the analysis for the manuscript. 

## Folders
List below are the major datasets and scripts that were used for the analysis. Note that the scripts for data exploratory are also included in the *scripts* folder, and so they might be a bit messy.

<pre>
\- data    
  |\- se.rda: `RangedSummariedExperiment` instance as a container containing gene counts, gene annotation and metadata. Unfiltered.  
  |\- sanitized.dds.rda:  sanitied/filtered `DESeqDataSet` used for downstream analysis  
  |\- sanitized.rlg.rda: sanitized/filtered rlog-transformed dataset used to visualize sample space    
\- scripts      
  |\- \*.Rmd:  R markdown document and notebook  
  |\- makeSE.R: get gene count and make a `RangedSummarizedExperiment` instance  
  |\- \*.R: scripts for exploratory analysis and for manuscript  
</pre>
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
