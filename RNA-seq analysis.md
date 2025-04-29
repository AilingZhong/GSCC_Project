# RNA-seq analysis

```R

%%%%%%%%%%%%%%%%%%%%%%%%sgScr_vs_sgEZH2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%sgScr_vs_sgEZH2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%sgScr_vs_sgEZH2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%sgScr_vs_sgEZH2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vim sample_fq.csv
sgScr_1 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Scr_1_FRAS210211291-1r/Cas9_Scr_1_FRAS210211291-1r
sgScr_2 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Scr_2_FRAS210211292-1r/Cas9_Scr_2_FRAS210211292-1r
sgScr_3 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Scr_3_FRAS210211293-1r/Cas9_Scr_3_FRAS210211293-1r
sgEzh2_g1_1 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Ezh2-1__1_FRAS210211294-1r/Cas9_Ezh2-1__1_FRAS210211294-1r
sgEzh2_g1_2 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Ezh2-1__2_FRAS210211295-1r/Cas9_Ezh2-1__2_FRAS210211295-1r
sgEzh2_g1_3 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Ezh2-1__3_FRAS210211296-1r/Cas9_Ezh2-1__3_FRAS210211296-1r
sgEzh2_g4_1 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Ezh2-4__1_FRAS210211297-1r/Cas9_Ezh2-4__1_FRAS210211297-1r
sgEzh2_g4_2 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Ezh2-4__2_FRAL210211298-1a/Cas9_Ezh2-4__2_FRAL210211298-1a
sgEzh2_g4_3 /mnt/data/sequencedata/RNAseq/RNAseq_93_ZMS_EZH2_20211108_18samples/CP2020120100066/H101SC21030878/RSCS3300/X101SC21030878-Z03/X101SC21030878-Z03-F046/2.cleandata/Cas9_Ezh2-4__3_FRAS210211299-1r/Cas9_Ezh2-4__3_FRAS210211299-1r

cat sample_fq.csv |while read id;
do echo $id
arr=($id)
fq2=${arr[1]}'_2.clean.fq.gz'
fq1=${arr[1]}'_1.clean.fq.gz'
sample=${arr[0]}
echo $fq2
echo $fq1
echo $sample
STAR  --readFilesCommand zcat  \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 20 \
--readFilesIn $fq1 $fq2 \
--genomeDir /mnt/data/public_data/genome_index/genome_index/mm10_genome_index_150 \
--outFileNamePrefix /mnt/data/user_data/abao/1_project/Bulk_project/1_bulk_RNA/ZhangMengsha_EZH2/$sample. ;
done


samtools index sgEzh2_g1_1.Aligned.sortedByCoord.out.bam
samtools index sgEzh2_g1_2.Aligned.sortedByCoord.out.bam
samtools index sgEzh2_g1_3.Aligned.sortedByCoord.out.bam
samtools index sgEzh2_g4_1.Aligned.sortedByCoord.out.bam
samtools index sgEzh2_g4_2.Aligned.sortedByCoord.out.bam
samtools index sgEzh2_g4_3.Aligned.sortedByCoord.out.bam
samtools index sgScr_1.Aligned.sortedByCoord.out.bam
samtools index sgScr_2.Aligned.sortedByCoord.out.bam
samtools index sgScr_3.Aligned.sortedByCoord.out.bam



project <- c("Normal")
deal_design <- c("sgEzh2","sgScr")
significant_cutoff <- c(1)
organism <- "mouse"

file_path <- paste(getwd(),"/Normal_sgEzh2_VS_sgScr_workfile",sep="")
sample_sampletable.path <- getwd()
dir.create(file_path)
load(file="/mnt/data/userdata/xiangyu/workshop/WORKFLOW_RNAseq/mm10_anno/ebg_mm10.RData")
load(file="/mnt/data/userdata/xiangyu/workshop/WORKFLOW_RNAseq/mm10_anno/txdb_mm10.RData")
deal_sgEzh2 <- deal_design[1]
deal_sgScr <- deal_design[2]
deal_names <- paste(deal_sgEzh2,deal_sgScr,sep="_VS_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"0_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"1_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"1_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"1_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"2_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"3",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"4",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"5",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"6",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"7",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"8",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"9",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})



indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sampletables.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",
singleEnd=FALSE,
ignore.strand=TRUE,
fragments=FALSE)
setwd(file_path)
save(se,file=se.save_RData)

colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL



pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()


pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()


colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","sgEzh2","sgScr"))
res_1 <- cbind(normalized_DEseq,DEseq_res)


if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
}

res_1$ENSEMBL <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$ENTREZID <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
all_summry <- cbind(countdata,res_1)
names(all_summry) <- c("sgEzh2_1","sgEzh2_2","sgEzh2_3","sgScr_1","sgScr_2","sgScr_3","DESeq2_sgEzh2_1","DESeq2_sgEzh2_2","DESeq2_sgEzh2_3","DESeq2_sgScr_1","DESeq2_sgScr_2","DESeq2_sgScr_3","baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "ENSEMBL","ENTREZID","GENENAME")
write.csv(all_summry, "Normal_DEseq2normalized_sgEzh2_VS_sgScr_allsummry.csv")


***********************Tumor***TPE_VS_TP*************************************************
***********************Tumor***TPE_VS_TP*************************************************
***********************Tumor***TPE_VS_TP*************************************************
***********************Tumor***TPE_VS_TP*************************************************
***********************Tumor***TPE_VS_TP*************************************************

STAR --runThreadN 20 \
--genomeDir /mnt/data/public_data/reference/Mus_musculus_UCSC/UCSC/mm10/star_index/mm10_150_index \
--readFilesCommand zcat \
--readFilesIn /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/TCP-GC-R1-0715_FRAS202410078-1r/TCP-GC-R1-0715_FRAS202410078-1r_2.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/TCP-GC-R1-0715_FRAS202410078-1r/TCP-GC-R1-0715_FRAS202410078-1r_1.clean.fq.gz \
--outFileNamePrefix /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/TuOrg_TPS_1_ \
--outSAMtype BAM SortedByCoordinate

STAR --runThreadN 20 \
--genomeDir /mnt/data/public_data/reference/Mus_musculus_UCSC/UCSC/mm10/star_index/mm10_150_index \
--readFilesCommand zcat \
--readFilesIn /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/R1-NR-TPV2TC_FRAS202410077-1r/R1-NR-TPV2TC_FRAS202410077-1r_2.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/R1-NR-TPV2TC_FRAS202410077-1r/R1-NR-TPV2TC_FRAS202410077-1r_1.clean.fq.gz \
--outFileNamePrefix /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/TuOrg_TPS_2_ \
--outSAMtype BAM SortedByCoordinate


STAR --runThreadN 20 \
--genomeDir /mnt/data/public_data/reference/Mus_musculus_UCSC/UCSC/mm10/star_index/mm10_150_index \
--readFilesCommand zcat \
--readFilesIn /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/R3-LRR-TPE-0921_FRAS202410079-1r/R3-LRR-TPE-0921_FRAS202410079-1r_2.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/R3-LRR-TPE-0921_FRAS202410079-1r/R3-LRR-TPE-0921_FRAS202410079-1r_1.clean.fq.gz \
--outFileNamePrefix /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/ZMS_TuOrg_TPEzh2_1_ \
--outSAMtype BAM SortedByCoordinate

STAR --runThreadN 20 \
--genomeDir /mnt/data/public_data/reference/Mus_musculus_UCSC/UCSC/mm10/star_index/mm10_150_index \
--readFilesCommand zcat \
--readFilesIn /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/L2-LRR-TPE-0921_FRAS202410080-1r/L2-LRR-TPE-0921_FRAS202410080-1r_2.clean.fq.gz /mnt/data/sequencedata/RNAseq/RNAseq_67_LHY_gastric_cancer_20201101_7samples/CP2018091214925/H101SC20030543/RSCS7000/X101SC20030543-Z01/X101SC20030543-Z01-J051/2.cleandata/L2-LRR-TPE-0921_FRAS202410080-1r/L2-LRR-TPE-0921_FRAS202410080-1r_1.clean.fq.gz \
--outFileNamePrefix /mnt/data/userdata/abao/project/1_bulk_sequence/11_Liuhongyu/ZMS_TuOrg_TPEzh2_2_ \
--outSAMtype BAM SortedByCoordinate



project <- c("Tumor")
deal_design <- c("TPE","TP")
significant_cutoff <- c(1)
organism <- "mouse"

file_path <- paste(getwd(),"/Tumor_TPE_VS_TP_workfile",sep="")
sample_sampletable.path <- getwd()
dir.create(file_path)
load(file="/mnt/data/userdata/xiangyu/workshop/WORKFLOW_RNAseq/mm10_anno/ebg_mm10.RData")
load(file="/mnt/data/userdata/xiangyu/workshop/WORKFLOW_RNAseq/mm10_anno/txdb_mm10.RData")
deal_TPE <- deal_design[1]
deal_TP <- deal_design[2]
deal_names <- paste(deal_TPE,deal_TP,sep="_VS_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"0_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"1_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"1_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"1_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"2_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"3",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"4",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"5",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"6",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"7",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"8",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"9",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})



indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sampletables.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",
singleEnd=FALSE,
ignore.strand=TRUE,
fragments=FALSE)
setwd(file_path)
save(se,file=se.save_RData)

colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL



pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()


pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()


colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
DEseq_res <- results(DESeq_counts,contrast=c("deal","TPE","TP"))
res_1 <- cbind(normalized_DEseq,DEseq_res)


if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
}

res_1$ENSEMBL <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
res_1$ENTREZID <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(res_1)
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
all_summry <- cbind(countdata,res_1)
names(all_summry) <- c("TPE_1","TPE_2","TP_1","TP_2","DESeq2_TPE_1","DESeq2_TPE_2","DESeq2_TP_1","DESeq2_TP_2","baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj", "ENSEMBL","ENTREZID","GENENAME")
write.csv(all_summry, "Tumor_DEseq2normalized_TPE_VS_TP_allsummry.csv")



%%%%%%%%%%%%%%%%%%%%%%%%TPE_sgScr_vs_TPE_sgTfap2c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%TPE_sgScr_vs_TPE_sgTfap2c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%TPE_sgScr_vs_TPE_sgTfap2c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%TPE_sgScr_vs_TPE_sgTfap2c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vim fq_list.csv
TPE_Scr_rep1 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_Scr_rep1/E-Scr__1_
TPE_Scr_rep2 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_Scr_rep2/E-Scr__2_
TPE_Scr_rep3 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_Scr_rep3/E-Scr__3_
TPE_TFAP2C_sgRNA1_rep1 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_TFAP2C_sgRNA1_rep1/2C-1__1_
TPE_TFAP2C_sgRNA1_rep2 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_TFAP2C_sgRNA1_rep2/2C-1__2_
TPE_TFAP2C_sgRNA1_rep3 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_TFAP2C_sgRNA1_rep3/2C-1__3_
TPE_TFAP2C_sgRNA2_rep1 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_TFAP2C_sgRNA2_rep1/2C-2__1_
TPE_TFAP2C_sgRNA2_rep2 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_TFAP2C_sgRNA2_rep2/2C-2__2_
TPE_TFAP2C_sgRNA2_rep3 /mnt/data/sequencedata/RNAseq/RNAseq_104_ZMS_EZH2_20220530_9samples/CP2020120100066/H101SC21030878/RSHR00104/X101SC21030878-Z03/X101SC21030878-Z03-J110/00.CleanData/TPE_TFAP2C_sgRNA2_rep3/2C-2__3_

cat fq_list.csv | while read id;
do echo $id
arr=($id)
fq2=${arr[1]}'_2.clean.fq.gz'
fq1=${arr[1]}'_1.clean.fq.gz'
sample=${arr[0]}
echo $fq1
echo $fq2
echo $sample
STAR --runThreadN 20 --outSAMtype BAM SortedByCoordinate \
--genomeDir /mnt/data/public_data/genome_index/genome_index/mm10_genome_index_150/ \
--readFilesCommand zcat \
--readFilesIn $fq1 $fq2 \
--outFileNamePrefix /mnt/data/user_data/abao/1_project/Bulk_project/1_bulk_RNA/ZhangMengsha_EZH2/$sample. ;
done



samtools index -@ 20 TPE_TFAP2C_sgRNA1_rep1.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_TFAP2C_sgRNA1_rep2.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_TFAP2C_sgRNA1_rep3.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_TFAP2C_sgRNA2_rep1.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_TFAP2C_sgRNA2_rep2.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_TFAP2C_sgRNA2_rep3.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_Scr_rep1.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_Scr_rep2.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TPE_Scr_rep3.Aligned.sortedByCoord.out.bam
samtools index -@ 20 TuOrg_TPS_1_Aligned.sortedByCoord.out.bam
samtools index -@ 20 TuOrg_TPS_2_Aligned.sortedByCoord.out.bam
samtools index -@ 20 ZMS_TuOrg_TPEzh2_1_Aligned.sortedByCoord.out.bam
samtools index -@ 20 ZMS_TuOrg_TPEzh2_2_Aligned.sortedByCoord.out.bam


.libPaths( "/mnt/data/user_data/xiangyu/programme/R_PACKAGES/R/x86_64-pc-linux-gnu-library/3.6")

project <- c("EZH2_Project")
deal_design <- c("sgTfap2c","sgScr")
significant_cutoff <- c(1)
organism <- "mouse"
file_path <- paste(getwd(),"/GSCC_sgTfap2c_workfile",sep="")
sample_sampletable.path <- getwd()
dir.create(file_path)
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/txdb_mm10.RData")
load(file="/mnt/data/user_data/xiangyu/workshop/WORKFLOW_RNAseq/ebg_mm10.RData")
deal_sgTfap2c	 <- deal_design[1]
deal_sgScr <- deal_design[2]
deal_names <- paste(deal_sgTfap2c,deal_sgScr,sep="_VS_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"0_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"1_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"1_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"1_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"2_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"3",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"4",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"5",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"6",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"7",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"8",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"9",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})



indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sampletable.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",
singleEnd=FALSE,
ignore.strand=TRUE,
fragments=FALSE)
setwd(file_path)
save(se,file=se.save_RData)



colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)
DESeq_counts <- DESeq(dds)
normalized_DEseq <- counts(DESeq_counts, normalized=TRUE)
rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL

pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()


pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()

normalized_DEseq <- data.frame(normalized_DEseq)
if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
} else {
	anno_data=org.Hs.eg.db
	print("organism is mouse")
}
normalized_DEseq$ENSEMBL <- mapIds(x = anno_data,
                        keys = rownames(normalized_DEseq),
						keytype ="SYMBOL",
						column ="ENSEMBL",
						multiVals="first")
normalized_DEseq$entrez <- mapIds(x = anno_data,
                        keys = rownames(normalized_DEseq),
						keytype ="SYMBOL",
						column ="ENTREZID",
						multiVals="first")
AA <- rownames(normalized_DEseq)
AA <- as.character(AA)
normalized_DEseq$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
all_summry <- cbind(countdata,normalized_DEseq)
names(all_summry) <- c("TFAP2C_sgRNA1_rep1","TFAP2C_sgRNA1_rep2","TFAP2C_sgRNA1_rep3","TFAP2C_sgRNA2_rep1","TFAP2C_sgRNA2_rep2","TFAP2C_sgRNA2_rep3","Scr_1","Scr_2","Scr_3",
	"DESeq2_TFAP2C_sg1_rep1","DESeq2_TFAP2C_sg1_rep2","DESeq2_TFAP2C_sg1_rep3",	"DESeq2_TFAP2C_sg2_rep1","DESeq2_TFAP2C_sg2_rep2","DESeq2_TFAP2C_sg2_rep3",
	"DESeq2_Scr_1","DESeq2_Scr_2","DESeq2_Scr_3",
	"ENSEMBL","entrez","GENENAME")
write.csv(all_summry, "All_DEseq2normalized_9samples_of_GSCC_Project.csv")
```

