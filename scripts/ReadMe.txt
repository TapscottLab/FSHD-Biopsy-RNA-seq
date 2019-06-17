## Raw Data

I:   `/shared/ngs/illumina/lsnider/151002_SN367_0557_AHC5HLBCXX` 

II:  `/shared/ngs/illumina/lsnider/160208_SN367_0613_BHJFJ3BCXX`

III: `/shared/ngs/illumina/lsnider/160908_D00300_0314_AC9NKGANXX`

IV: `/shared/ngs/illumina/lsnider/170119_D00300_0372_BH77YFBCXY`

Discovery: `/shared/ngs/illumina/apfong/130419_SN367_0277_AC2192ACXX` and
`/shared/ngs/illumina/apfong/130426_SN367_0278_AD26C4ACXX`.


## Work Flow

The preprocess pipeline begins with removing low quality reads that
did not pass the CHASTITY filtering step and then trimming the TruSeq
adapter by using Trimmomatic. The trimmed reads are aligned to hg38 by
Tophat2/BowTie2. 

1.  Remove reads that did not pass the CHASTITY filtering step
2.  Adapter trimming: trimmomatic-0.32
java -jar /home/solexa/apps/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 6 $fqFile $trimFq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 MINLEN:36
3.  Trim deteriorated tails
4.  Fastqc again
5. Alignment: Tophat2/Bowtie2

(step 2 - 5: do_tophat.sh)

6. Checking fastqc by eyes
7. makeSE: sanitize samples, gene counts, input, mapping rate, lib_size and sizeFactor
8. Rank !
9. Produce tables and figures

makeSE.R -> rankScore.R -> makeManuscripts_scatterPlot.R

## The forth and Fifth batch for second visit
The raw reads are here: /fh/fast/tapscott_s/SR/ngs/illumina/acampbel/171025_D00300_0490_ACBY30ANXX/Unaligned/Project_acampbel/
~                                       
