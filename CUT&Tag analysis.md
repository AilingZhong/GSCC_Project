# H3K27me3 CUT&Tag analysis

```shell
*********************H3K27me3_in_TP_and_TPE************************
*********************H3K27me3_in_TP_and_TPE************************
*********************H3K27me3_in_TP_and_TPE************************
*********************H3K27me3_in_TP_and_TPE************************

/mnt/data/sequencedata/CUTTAG/CUT_tag_6_ZMS_Ezh2_stomach_4samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J001/2.cleandata

fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPS1_H3K27me3_FKDL210000787-1a_1.clean.fq.gz
fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPS1_H3K27me3_FKDL210000787-1a_2.clean.fq.gz
fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPS2_H3K27me3_FKDL210000788-1a_1.clean.fq.gz
fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPS2_H3K27me3_FKDL210000788-1a_2.clean.fq.gz
fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPE1_H3K27me3_FKDL210000789-1a_1.clean.fq.gz
fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPE1_H3K27me3_FKDL210000789-1a_2.clean.fq.gz
fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPE2_H3K27me3_FKDL210000790-1a_1.clean.fq.gz
fastqc -t 15 -o /mnt/data/user_data/zhaolei/project/zms_EZH2/chip_H3K27me3/fastqFileQC \
       -f fastq ZMS_TPE2_H3K27me3_FKDL210000790-1a_2.clean.fq.gz


vim All_H3K27me3_samplelist.csv
TP_H3K27me3_rep1 /mnt/data2/userdata/abao/8922/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/fastqFileQC/ZMS_TPS1_H3K27me3_FKDL210000787-1a
TP_H3K27me3_rep2 /mnt/data2/userdata/abao/8922/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/fastqFileQC/ZMS_TPS2_H3K27me3_FKDL210000788-1a
TPE_H3K27me3_rep1 /mnt/data2/userdata/abao/8922/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/fastqFileQC/ZMS_TPE1_H3K27me3_FKDL210000789-1a
TPE_H3K27me3_rep2 /mnt/data2/userdata/abao/8922/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/fastqFileQC/ZMS_TPE2_H3K27me3_FKDL210000790-1a

cat All_H3K27me3_samplelist.csv | while read id ; do
arr=($id)
fq2=${arr[1]}'_2.clean.fq.gz'
fq1=${arr[1]}'_1.clean.fq.gz'
sample=${arr[0]}
echo $fq2
echo $fq1
echo $sample 
bowtie2 -p 20 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
-x /mnt/data/user_data/abao/2_reference/mm10/mm10 \
-1 $fq1 -2 $fq2 | samtools sort -O bam -@ 20 -o - > $sample.bowtie2.sorted.bam ; 
done


mkdir -p filter_bam/filter_MT_bam
mkdir bed
mkdir peak_file

cat All_H3K27me3_samplelist.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample 
/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk --java-options "-Xmx30G -Djava.io.tmpdir=./" MarkDuplicates \
-I ${sample}.bowtie2.sorted.bam \
-O ./filter_bam/${sample}.filtdup.bam \
-M ./filter_bam/${sample}.dups.txt \
-REMOVE_DUPLICATES=true ;
done


cat All_H3K27me3_samplelist.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $fq2
echo $fq1
echo $sample 
samtools index  ./filter_bam/$sample.filtdup.bam
mtReads=$(samtools idxstats  ./filter_bam/$sample.filtdup.bam | grep 'chrM' | cut -f 3)
totalReads=$(samtools idxstats  ./filter_bam/$sample.filtdup.bam | awk '{SUM += $3} END {print SUM}')
echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'
samtools flagstat  ./filter_bam/$sample.filtdup.bam > $sample.rmdup.stat
samtools view  -h  -f 2 -q 30 ./filter_bam/$sample.filtdup.bam |grep -v chrM > ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam
samtools sort ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam  -@ 8 -o - > ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam
rm -r ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam
samtools index ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam 
samtools flagstat  ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam > $sample.filter_dupli_chrM_last.stat 
bedtools bamtobed -i ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam  > ./bed/$sample.bed ;
done 


ls *filter_dupli_chrM_last.bam |while read id;do
bamCoverage -p 20 -bs=1 --normalizeUsing BPM -b $id -o ./${id%%.*}.bs1.last.bw ;
done



******************************Merge组内重复************************************************
******************************Merge组内重复************************************************
******************************Merge组内重复************************************************
******************************Merge组内重复************************************************

cd /mnt/data2/userdata/abao/8922/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/
samtools merge -@ 20 Merge_TPE_H3K27me3.sorted.bam renew_TPE2_H3K27me3.filter_dup_chrM.bam renew_TPE1_H3K27me3.filter_dup_chrM.bam
samtools merge -@ 20 Merge_TP_H3K27me3.sorted.bam renew_TP2_H3K27me3.filter_dup_chrM.bam renew_TP1_H3K27me3.filter_dup_chrM.bam

samtools sort -@ 20 Merge_TPE_H3K27me3.sorted.bam -o Merge_TPE_H3K27me3.sorted.bam
samtools index Merge_TPE_H3K27me3.sorted.bam
bamCoverage -p 25 -bs=1 --normalizeUsing BPM -b Merge_TPE_H3K27me3.sorted.bam -o Merge_TPE_H3K27me3.bs1.bw

samtools sort -@ 20 Merge_TP_H3K27me3.sorted.bam -o Merge_TP_H3K27me3.sorted.bam
samtools index Merge_TP_H3K27me3.sorted.bam
bamCoverage -p 25 -bs=1 --normalizeUsing BPM -b Merge_TP_H3K27me3.sorted.bam -o Merge_TP_H3K27me3.bs1.bw


```

```shell

*****************************H3K27me3_peak_calling************************************************
*****************************H3K27me3_peak_calling************************************************
*****************************H3K27me3_peak_calling************************************************
*****************************H3K27me3_peak_calling************************************************

macs2 callpeak -t renew_TPE2_H3K27me3.filter_dup_chrM.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_TPE_H3K27me3_rep2 \
-n TPE_H3K27me3_rep2

macs2 callpeak -t renew_TPE1_H3K27me3.filter_dup_chrM.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_TPE_H3K27me3_rep1 \
-n TPE_H3K27me3_rep1

macs2 callpeak -t renew_TP2_H3K27me3.filter_dup_chrM.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_TP_H3K27me3_rep2 \
-n TP_H3K27me3_rep2

macs2 callpeak -t renew_TP1_H3K27me3.filter_dup_chrM.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_TP_H3K27me3_rep1 \
-n TP_H3K27me3_rep1
```

```R
&&&&&&&&&&&&&&&&&&&&&Annotation&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&Annotation&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&Annotation&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicFeatures)
library(ReactomePA)
library(AnnotationDbi)
library(DOSE)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

library(future.apply)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

files <- list.files(path = "/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/q001_H3K27me3_macs2_peaks", pattern = "_peaks.narrowPeak$", full.names = TRUE)
name <- gsub("_peaks.narrowPeak","",basename(files))

for (i in c(1:length(files))){
  print(i)
  peakAnno <- annotatePeak(files[[i]], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db") 
  peakAnno <- as.data.frame(peakAnno)
  aa <- name[i]
  hh <- paste(aa,"annotatePeak.csv",sep="_")
  bb <- "/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/q001_H3K27me3_macs2_peaks"
  cc <- paste(bb,hh,sep="/")
  write.csv(peakAnno, file=cc)
}

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE)

pdf("K27me3_plotAnnoBar.pdf")
plotAnnoBar(peakAnnoList)
dev.off()


TPE_VS_TP_H3K27me3 <- read.csv("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/1_TPE_VS_TP_H3K27me3_all_anno_peaks.csv")
TPE_VS_TP_H3K27me3 <- TPE_VS_TP_H3K27me3[,c("seqnames","start","end")]
write.table(TPE_VS_TP_H3K27me3,"TPE_VS_TP_H3K27me3_all.txt",row.names =FALSE, col.names =FALSE,sep="\t",quote =FALSE)

cd /mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/
bed=/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/TPE_VS_TP_H3K27me3_all.bed
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
-R $bed \
-S Merge_TP_H3K27me3.bs1.bw Merge_TPE_H3K27me3.bs1.bw \
--numberOfProcessors 30 --skipZeros -o All_peaks_H3K27me3_levels.mat.gz 

plotHeatmap -m All_peaks_H3K27me3_levels.mat.gz -out All_peaks_H3K27me3_levels.png \
--colorList  "#edf8fb,#b3cde3,#8c96c6,#88419d" "#edf8fb,#b3cde3,#8c96c6,#88419d" \
--samplesLabel "TP" "TPE" \
--plotFileFormat png


```

```R

*********************H3K27me3_PCA_analysis**********************************
*********************H3K27me3_PCA_analysis**********************************
*********************H3K27me3_PCA_analysis**********************************
*********************H3K27me3_PCA_analysis**********************************

library("DESeq2")
Counts_H3K27me3 <- read.csv("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/1_TPE_VS_TP_H3K27me3_all_anno_peaks.csv")
#这是抽取他们merge的所有common的peaks开始定量计算每一个region的counts
library(dplyr)
library(limma)
library(trqwe
library(DESeq2)
library(ggplot2)
consensusToCount <- readRDS("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/TPE_VS_TP_H3K27me3_consensusToCount.rds")
counts <- readRDS("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/TPE_VS_TP_H3K27me3_countsFromATAC.rds")

Group <- factor(c("TPE_H3K27me3","TPE_H3K27me3","TP_H3K27me3","TP_H3K27me3"))
metaData <- data.frame(Group, row.names = colnames(counts))
metaData$type <- rownames(metaData)
DDS <- DESeqDataSetFromMatrix(counts, metaData, ~Group, rowRanges = consensusToCount)
rld <- rlog(DDS, blind=FALSE)
ff <- plotPCA(rld, intgroup=c("Group", "type"))
ggsave(ff,file="H3K27me3_rld_PCA_plot.svg",height=6,width=6)
```

```R
****************Differential_modification_region_analysis************************
****************Differential_modification_region_analysis************************
****************Differential_modification_region_analysis************************
****************Differential_modification_region_analysis************************

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicFeatures)
library(ReactomePA)
library(AnnotationDbi)
library(DOSE)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

TPE_VS_TP_H3K27me3 <- read.csv("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/1_TPE_VS_TP_H3K27me3_all_anno_peaks.csv")
TP_hi <- subset(TPE_VS_TP_H3K27me3,pvalue < 0.05 & log2FoldChange < 0)
TPE_hi <- subset(TPE_VS_TP_H3K27me3,pvalue < 0.05 & log2FoldChange > 0)


TP_hi <- TP_hi[,c("seqnames","start","end")]
write.table(TP_hi,"TP_hi_all.txt",row.names =FALSE, col.names =FALSE,sep="\t",quote =FALSE)

TPE_hi <- TPE_hi[,c("seqnames","start","end")]
write.table(TPE_hi,"TPE_hi_all.txt",row.names =FALSE, col.names =FALSE,sep="\t",quote =FALSE)

/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/TPE_hi_all_log0.bed
/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/TP_hi_all_log0.bed

TPE_hi <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/TPE_hi_all_log0.bed")
TP_hi <- readPeakFile("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/TP_hi_all_log0.bed")
TPE_atacDDS_results <- annotatePeak(TPE_hi, tssRegion=c(-3000, 3000),TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
pdf("plotAnnoPie_TPE_hi.pdf",width=7,height=7)
plotAnnoPie(TPE_atacDDS_results)
dev.off()

TP_atacDDS_results <- annotatePeak(TP_hi, tssRegion=c(-3000, 3000),TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")
pdf("plotAnnoPie_TP_hi.pdf",width=7,height=7)
plotAnnoPie(TP_atacDDS_results)
dev.off()



TPE_VS_TP_H3K27me3 <- read.csv("/mnt/data/userdata/abao/project/6_CUTTAG/4_ZMS_EZH2_project/H3K27me3_filter_MT_bam/filter_bam/filter_chrM_bam/1_TPE_VS_TP_H3K27me3_all_anno_peaks.csv")
TP_hi <- subset(TPE_VS_TP_H3K27me3,pvalue < 0.05 & log2FoldChange < 0)
TPE_hi <- subset(TPE_VS_TP_H3K27me3,pvalue < 0.05 & log2FoldChange > 0)

TP_Distal <- subset(TP_hi,annotation=="Distal Intergenic")
TP_Promoter <- subset(TP_hi,annotation=="Promoter (<=1kb)" | annotation=="Promoter (1-2kb)" | annotation=="Promoter (2-3kb)")
TP_Intron <- subset(TP_hi, grepl("^Intron", TP_hi$annotation))
TP_Exon <- subset(TP_hi, grepl("^Exon", TP_hi$annotation))
TP_3UTR <- subset(TP_hi, grepl("^3", TP_hi$annotation))
TP_5UTR <- subset(TP_hi, grepl("^5", TP_hi$annotation))
TP_DN <- subset(TP_hi, grepl("^Down", TP_hi$annotation))

TPE_Distal <- subset(TPE_hi,annotation=="Distal Intergenic")
TPE_Promoter <- subset(TPE_hi,annotation=="Promoter (<=1kb)" | annotation=="Promoter (1-2kb)" | annotation=="Promoter (2-3kb)")
TPE_Intron<- subset(TPE_hi, grepl("^Intron", TPE_hi$annotation))
TPE_Exon <- subset(TPE_hi, grepl("^Exon", TPE_hi$annotation))
TPE_3UTR <- subset(TPE_hi, grepl("^3", TPE_hi$annotation))
TPE_5UTR <- subset(TPE_hi, grepl("^5", TPE_hi$annotation))
TPE_DN <- subset(TPE_hi, grepl("^Down", TPE_hi$annotation))

TP_1 <- data.frame(c(length(TP_Distal$region),length(TP_Promoter$region),length(TP_Intron$region),length(TP_Exon$region),length(TP_3UTR$region),length(TP_5UTR$region),length(TP_DN$region)))
names(TP_1) <- "Num"
TP_1$Type <-c("Distal","Promoter","Intron","Exon","3UTR","5UTR","DN")
TP_1$Group <-c("TP")


TPE_1 <- data.frame(c(length(TPE_Distal$region),length(TPE_Promoter$region),length(TPE_Intron$region),length(TPE_Exon$region),length(TPE_3UTR$region),length(TPE_5UTR$region),length(TPE_DN$region)))
names(TPE_1) <- "Num"
TPE_1$Type <-c("Distal","Promoter","Intron","Exon","3UTR","5UTR","DN")
TPE_1$Group <-c("TPE")

tmp <- rbind(TP_1,TPE_1)

# 根据性别对num赋值
tmp$Group <- factor(tmp$Group,levels = c("TP","TPE"))

tmp$Type <- factor(tmp$Type,levels = c("Promoter","Distal","Intron","Exon","3UTR","5UTR","DN"))

ff <- ggplot(tmp,aes(x=Type, y=Num,fill=Group))+
  geom_bar(stat='identity',position=position_dodge())+
  scale_fill_manual(values=c("TP" = "#cb181d","TPE" = "#2b8cbe"))
```



# EZH2 and TFAP2C CUT&Tag analysis

```shell
******************fastp_QC**EZH2**TFAP2C*CUTtag************************
******************fastp_QC**EZH2**TFAP2C*CUTtag************************
******************fastp_QC**EZH2**TFAP2C*CUTtag************************
******************fastp_QC**EZH2**TFAP2C*CUTtag************************
******************fastp_QC**EZH2**TFAP2C*CUTtag************************

/mnt/data/user_data/abao/4_packages/fastp -i /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/S-E-3/S-E-3_2.fq.gz \
-o /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep3_2.fq.gz \
-I /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/S-E-3/S-E-3_1.fq.gz \
-O /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep3_1.fq.gz \
--thread=15 -c

/mnt/data/user_data/abao/4_packages/fastp -i /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/S-E-2/S-E-2_2.fq.gz \
-o /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep2_2.fq.gz \
-I /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/S-E-2/S-E-2_1.fq.gz \
-O /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep2_1.fq.gz \
--thread=15 -c

/mnt/data/user_data/abao/4_packages/fastp -i /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/S-E-1/S-E-1_2.fq.gz \
-o /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep1_2.fq.gz \
-I /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/S-E-1/S-E-1_1.fq.gz \
-O /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep1_1.fq.gz \
--thread=15 -c

/mnt/data/user_data/abao/4_packages/fastp -i /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/E-2C-3/E-2C-3_1.fq.gz \
-o /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep3_2.fq.gz \
-I /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/E-2C-3/E-2C-3_2.fq.gz \
-O /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep3_1.fq.gz \
--thread=15 -c


/mnt/data/user_data/abao/4_packages/fastp -i /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/E-2C-2/E-2C-2_1.fq.gz \
-o /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep2_2.fq.gz \
-I /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/E-2C-2/E-2C-2_2.fq.gz \
-O /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep2_1.fq.gz \
--thread=15 -c

/mnt/data/user_data/abao/4_packages/fastp -i /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/E-2C-1/E-2C-1_1.fq.gz \
-o /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep1_2.fq.gz \
-I /mnt/data/sequencedata/CUTTAG/CUT_tag_23_ZMS_EZH2_20221026_12samples/CP2020120100066/H101SC20122270/RSCS0500/X101SC20122270-Z01/X101SC20122270-Z01-J066/01.RawData/E-2C-1/E-2C-1_2.fq.gz \
-O /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep1_1.fq.gz \
--thread=15 -c



*********************mm10_bowtie2比对并输出结果文件****************************************************
*********************mm10_bowtie2比对并输出结果文件****************************************************
*********************mm10_bowtie2比对并输出结果文件****************************************************
*********************mm10_bowtie2比对并输出结果文件****************************************************
*********************mm10_bowtie2比对并输出结果文件****************************************************

vim All_samplelist.csv
TP_antiEZH2_rep1 /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep1
TP_antiEZH2_rep2 /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep2
TP_antiEZH2_rep3 /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/TP_antiEZH2_rep3
TPE_antiTFAP2C_rep1 /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep1
TPE_antiTFAP2C_rep2 /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep2
TPE_antiTFAP2C_rep3 /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/TPE_antiTFAP2C_rep3

cat All_samplelist.csv | while read id ; do
arr=($id)
fq2=${arr[1]}'_2.fq.gz'
fq1=${arr[1]}'_1.fq.gz'
sample=${arr[0]}
echo $fq2
echo $fq1
echo $sample 
bowtie2 -p 20 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
-x /mnt/data/user_data/abao/2_reference/mm10/mm10 \
-1 $fq1 -2 $fq2 | samtools sort -O bam -@ 20 -o - > $sample.bowtie2.sorted.bam ; 
done


mkdir -p filter_bam/filter_MT_bam
mkdir bed
mkdir peak_file

cat All_samplelist.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample 
/mnt/data/user_data/xiangyu/programme/gatk-4.1.3.0/gatk --java-options "-Xmx30G -Djava.io.tmpdir=./" MarkDuplicates \
-I ${sample}.bowtie2.sorted.bam \
-O ./filter_bam/${sample}.filtdup.bam \
-M ./filter_bam/${sample}.dups.txt \
-REMOVE_DUPLICATES=true ;
done


cat All_samplelist.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $fq2
echo $fq1
echo $sample 
samtools index  ./filter_bam/$sample.filtdup.bam
mtReads=$(samtools idxstats  ./filter_bam/$sample.filtdup.bam | grep 'chrM' | cut -f 3)
totalReads=$(samtools idxstats  ./filter_bam/$sample.filtdup.bam | awk '{SUM += $3} END {print SUM}')
echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'
samtools flagstat  ./filter_bam/$sample.filtdup.bam > $sample.rmdup.stat
samtools view  -h  -f 2 -q 30 ./filter_bam/$sample.filtdup.bam |grep -v chrM > ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam
samtools sort ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam  -@ 8 -o - > ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam
rm -r ./filter_bam/filter_MT_bam/$sample.filter_dupli.tmp.bam
samtools index ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam 
samtools flagstat  ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam > $sample.filter_dupli_chrM_last.stat 
bedtools bamtobed -i ./filter_bam/filter_MT_bam/$sample.filter_dupli_chrM_last.bam  > ./bed/$sample.bed ;
done 


ls *filter_dupli_chrM_last.bam |while read id;do
bamCoverage -p 20 -bs=1 --normalizeUsing BPM -b $id -o ./${id%%.*}.bs1.last.bw ;
done

******************************Merge组内重复************************************************
******************************Merge组内重复************************************************
******************************Merge组内重复************************************************
******************************Merge组内重复************************************************

cd /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/filter_bam/filter_MT_bam/TFAP2C_bam/
samtools merge -@ 20 Merge_TPE_antiTFAP2C.bam TPE_antiTFAP2C_rep3.filter_dupli_chrM_last.bam TPE_antiTFAP2C_rep2.filter_dupli_chrM_last.bam TPE_antiTFAP2C_rep1.filter_dupli_chrM_last.bam
samtools sort -@ 20 Merge_TPE_antiTFAP2C.bam -o Merge_TPE_antiTFAP2C.sorted.bam
samtools index Merge_TPE_antiTFAP2C.sorted.bam
bamCoverage -p 25 -bs=1 --normalizeUsing BPM -b Merge_TPE_antiTFAP2C.sorted.bam -o Merge_TPE_antiTFAP2C.bs1.bw

cd /mnt/data/user_data/abao/1_project/Bulk_project/5_CUTTAG/5_ZMS_GC/ZMS_EZH2_TFAP2C/filter_bam/filter_MT_bam/EZH2_bam/
samtools merge -@ 20 Merge_TP_antiEZH2.bam TP_antiEZH2_rep3.filter_dupli_chrM_last.bam TP_antiEZH2_rep2.filter_dupli_chrM_last.bam TP_antiEZH2_rep1.filter_dupli_chrM_last.bam
samtools sort -@ 20 Merge_TP_antiEZH2.bam -o Merge_TP_antiEZH2.sorted.bam
samtools index Merge_TP_antiEZH2.sorted.bam
bamCoverage -p 25 -bs=1 --normalizeUsing BPM -b Merge_TP_antiEZH2.sorted.bam -o Merge_TP_antiEZH2.bs1.bw

```

```shell
**************************MACS2_call_peaks****************************************************
**************************MACS2_call_peaks****************************************************
**************************MACS2_call_peaks****************************************************
**************************MACS2_call_peaks****************************************************

macs2 callpeak -t TP_antiEZH2_rep1.filter_dupli_chrM_last.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_EZH2_macs2_peaks \
-n TP_antiEZH2_rep1

macs2 callpeak -t TP_antiEZH2_rep2.filter_dupli_chrM_last.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_EZH2_macs2_peaks \
-n TP_antiEZH2_rep2

macs2 callpeak -t TP_antiEZH2_rep3.filter_dupli_chrM_last.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_EZH2_macs2_peaks \
-n TP_antiEZH2_rep3


macs2 callpeak -t TPE_antiTFAP2C_rep1.filter_dupli_chrM_last.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_TPE_antiTFAP2C_rep1 \
-n TPE_antiTFAP2C_rep1

macs2 callpeak -t TPE_antiTFAP2C_rep2.filter_dupli_chrM_last.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_TPE_antiTFAP2C_rep2 \
-n TPE_antiTFAP2C_rep2

macs2 callpeak -t TPE_antiTFAP2C_rep3.filter_dupli_chrM_last.bam \
-q 0.01 \
-f BAMPE \
-g mm \
--keep-dup all \
--outdir q001_TPE_antiTFAP2C_rep3 \
-n TPE_antiTFAP2C_rep3
```

