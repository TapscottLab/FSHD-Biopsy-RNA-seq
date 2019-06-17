#!/bin/bash
#./do_tophat_D4Z4.sh
#SBATCH -n6 -t 4-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

fq1=$1
tophatOutDir=$2

base=`basename ${fq1}`
sampleName=`echo ${base//.trim.fastq}`
echo $sampleName

#use tophat and bowtie2
bowtieVersion=/home/solexa/apps/bowtie/bowtie2-2.2.6
tophatVersion=/home/solexa/apps/tophat/tophat-2.1.0.Linux_x86_64
genomeBuild=/fh/fast/tapscott_s/CompBio/D4Z4/Chr4Extended/Bowtie2Index/sequence
maxIntronLength=500000
innerDist=80
segmentLength=25
libraryType=fr-unstranded

#umask 0002
export PATH=/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/opt/moab/bin:$bowtieVersion:$tophatVersion:/home/solexa/apps/samtools/samtools-0.1.19:/home/solexa/apps/FastQC
#remove /app/bin from the PATH, as it will default to running bowtie 0.12.8
export PATH=${PATH/\/app\/bin:}

#define some directories
mkdir -p $tophatOutDir
bamDir=$tophatOutDir/bam
mkdir -p $bamDir
tmpDir=$tophatOutDir/$sampleName
mkdir -p $tmpDir


tophat --mate-inner-dist $innerDist --num-threads 6 --read-mismatches 4 --read-edit-dist 4 --library-type $libraryType -I $maxIntronLength --segment-length $segmentLength -o $tmpDir $genomeBuild $fq1

mv $tmpDir/accepted_hits.bam $bamDir/$sampleName.bam

# sort bam
cd $bamDir
samtools sort -@ 4 $sampleName.bam $sampleName.bam.sorted
mv $sampleName.bam.sorted.bam $sampleName.bam
samtools index $sampleName.bam
touch $sampleName.tophatDone.txt
exit 0



	 
