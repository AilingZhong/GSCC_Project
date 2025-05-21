# Data visualization of Figure 5

```R
EZH2_Binding_in_TP <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/CUTTag_data/1_EZH2_Binding_in_TP_all_anno_peaks.csv")
EZH2_bed <- EZH2_Binding_in_TP[,c(2:4)]
TPE_VS_TP_H3K27me3 <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/CUTTag_data/1_TPE_VS_TP_H3K27me3_all_anno_peaks.csv")
H3K27me3_decrease <- subset(TPE_VS_TP_H3K27me3,pvalue < 0.05 & log2FoldChange < 0)
common_genes <- data.frame(intersect(unique(na.omit(EZH2_Binding_in_TP$SYMBOL)),unique(na.omit(H3K27me3_decrease$SYMBOL))))
names(common_genes) <- "SYMBOL"
Tumor_RNAseq  <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/RNAseq_data/Tumor_DEseq2normalized_TPE_VS_TP_allsummry.csv")
CUTtag_in_RNA <- subset(Tumor_RNAseq, Tumor_RNAseq$X %in% common_genes$SYMBOL )
CUTtag_in_RNA <- CUTtag_in_RNA[order(-CUTtag_in_RNA$log2FoldChange),]
rownames(CUTtag_in_RNA) <- CUTtag_in_RNA$X


TPE_VS_TP_RNA <- CUTtag_in_RNA[,c("X","baseMean","log2FoldChange")]
names(TPE_VS_TP_RNA) <- c("SYMBOL","RNA_basemean","TPE_VS_TP_log2FoldChange")

EZH2_CUTtag_signal <- EZH2_Binding_in_TP[,c("SYMBOL","SignalValue")]
library(dplyr)
EZH2_CUTtag_signal <- EZH2_CUTtag_signal %>% dplyr::arrange(SYMBOL)
EZH2_CUTtag_signal_symbol_sum <- aggregate(x = EZH2_CUTtag_signal[,2],
                            by = list(EZH2_CUTtag_signal$SYMBOL),
                            FUN = sum)
head(EZH2_CUTtag_signal_symbol_sum)
names(EZH2_CUTtag_signal_symbol_sum) <- c("SYMBOL","EZH2_Binding_Signal")


H3K27me3_CUTtag_signal <- TPE_VS_TP_H3K27me3[,c("SYMBOL","TPE_H3K27me3_rep1","TPE_H3K27me3_rep2","TP_H3K27me3_rep1","TP_H3K27me3_rep2")]
library(dplyr)
H3K27me3_CUTtag_signal <- H3K27me3_CUTtag_signal %>% dplyr::arrange(SYMBOL)
H3K27me3_CUTtag_signal_symbol_sum <- aggregate(x = H3K27me3_CUTtag_signal[,2:5],
                            by = list(H3K27me3_CUTtag_signal$SYMBOL),
                            FUN = sum)
head(H3K27me3_CUTtag_signal_symbol_sum)
TP_VS_TPE <- (H3K27me3_CUTtag_signal_symbol_sum$TP_H3K27me3_rep1+ H3K27me3_CUTtag_signal_symbol_sum$TP_H3K27me3_rep2)/(H3K27me3_CUTtag_signal_symbol_sum$TPE_H3K27me3_rep1+ H3K27me3_CUTtag_signal_symbol_sum$TPE_H3K27me3_rep2)
log2FC <- log(TP_VS_TPE,2)
Basemean <- (H3K27me3_CUTtag_signal_symbol_sum$TPE_H3K27me3_rep1+ H3K27me3_CUTtag_signal_symbol_sum$TPE_H3K27me3_rep2+
H3K27me3_CUTtag_signal_symbol_sum$TP_H3K27me3_rep1+ H3K27me3_CUTtag_signal_symbol_sum$TP_H3K27me3_rep2)/4
H3K27me3_CUTtag_signal_final <- data.frame(SYMBOL= H3K27me3_CUTtag_signal_symbol_sum$Group.1, H3K27me3_Basemean= Basemean, TP_VS_TPE_log2FC=log2FC )


tmp <- merge(TPE_VS_TP_RNA, EZH2_CUTtag_signal_symbol_sum,by="SYMBOL")
tmp1 <- merge(tmp, H3K27me3_CUTtag_signal_final,by="SYMBOL")
tmp1[is.na(tmp1)] <- 0

tmp1$final_score <- tmp1$RNA_basemean * tmp1$TPE_VS_TP_log2FoldChange * tmp1$EZH2_Binding_Signal * tmp1$H3K27me3_Basemean * tmp1$TP_VS_TPE_log2FC

down <- subset(tmp1, final_score < 0 )
down$log2_score <- log2_score <- log(-1*down$final_score,2)*-1
up <- subset(tmp1, final_score > 0 )
up$log2_score <- log2_score <- log(up$final_score,2)

tmp_all <- rbind(down,up)
tmp_all <- tmp_all[order(-tmp_all$log2_score),]
tmp_all$label <- ifelse(tmp_all$SYMBOL=="Tfap2c", "TFAP2C","Others")
library(ggpubr)
ff <- ggdotchart(tmp_all, x = "SYMBOL", y = "log2_score",
           color = "label",                                # Color by groups
           palette = c("#00AFBB", "#E7B800"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = FALSE,                                # Rotate vertically
           dot.size = 1,                                 # Large dot size
           label = round(tmp_all$log2_score),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 0,
                             vjust = 0.3),               # Adjust label parameters
           ggtheme = theme_pubr())                 # ggplot2 theme
ggsave(ff,file="TFAP2C_rank_535.png",height=10,width=10)
```

![TFAP2C_rank_535](C:\Users\Irene\Desktop\NC_markdown_整理脚本\assets\TFAP2C_rank_535.png)

```R
STAD_FPKM_tumor <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/TCGA_data/STAD_FPKM_tumor.csv")
TCGA_STAD_CNV <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/TCGA_data/TCGA-STAD_CNV.csv")
TCGA_STAD_Mutation <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/TCGA_data/TCGA_STAD_Mutation.csv")

Tumor_TFAP2C <- subset(STAD_FPKM_tumor,SYMBOL=="TFAP2C")
rownames(Tumor_TFAP2C) <- Tumor_TFAP2C$SYMBOL
Tumor_TFAP2C <- Tumor_TFAP2C[,c(-1,-2,-378)]
data1 <- data.frame(t(Tumor_TFAP2C))
data1$Group <- "Tumor"
data1$Tumor_Sample_Barcode <- substring(rownames(data1),1,12)


EZH2_mut <- subset(TCGA_STAD_Mutation,Hugo_Symbol == "EZH2")


EZH2_mut <- data.frame(EZH2_mut$Tumor_Sample_Barcode)
names(EZH2_mut) <- "ID"
EZH2_mut$Type <- "Loss"

chr7 <- subset(TCGA_STAD_CNV,Chromosome=="7")
chr7_1 <- subset(chr7,Start > 148504475 & End < 148581383)
chr7_2 <- subset(chr7,Start < 148504475 & End > 148504475  & End < 148581383)
chr7_3 <- subset(chr7,Start > 148504475 & Start < 148581383 & End > 148581383)
chr7_4 <- subset(chr7,Start < 148504475 & End > 148581383)

EZH2_variation <- rbind(chr7_1,chr7_2,chr7_3,chr7_4)
EZH2_del <- subset(EZH2_variation,Segment_Mean <= -0.2)


EZH2_del_CNV <- substring(EZH2_del$Sample,1,12)
EZH2_del_CNV <- data.frame(EZH2_del_CNV)
names(EZH2_del_CNV) <- "ID"
EZH2_del_CNV$Type <- "Loss"
EZH2_var <- rbind(EZH2_mut,EZH2_del_CNV)
EZH2_var$Tumor_Sample_Barcode <- gsub("-", "\\.", EZH2_var$ID) 


EZH2_non_var <- data.frame(setdiff(data1$Tumor_Sample_Barcode,EZH2_var$Tumor_Sample_Barcode))
names(EZH2_non_var) <- "ID"
EZH2_non_var$Type <- "Intact"

EZH2_var <- EZH2_var[,c(3,2)]
names(EZH2_var) <- c("ID","Type")
EZH2_info <- rbind(EZH2_var,EZH2_non_var)

tmp_all <- merge(EZH2_info,data1,by.x="ID",by.y="Tumor_Sample_Barcode")
tmp_all$log2_TFAP2C <- log(tmp_all$TFAP2C+1,2)

library(ggpubr)
library(ggplot2)
my_comparisons <- list(c("Intact", "Loss"))
ff <- ggboxplot(tmp_all, x = "Type", y = "log2_TFAP2C",add = "jitter",
               color = "Type",palette = c("#08519c","#e31a1c")) +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="TFAP2C_exp_in_EZH2_variant_or_WT_genome.png",width=4,height=6.4)
```

![TFAP2C_exp_in_EZH2_variant_or_WT_genome](C:\Users\Irene\Desktop\NC_markdown_整理脚本\assets\TFAP2C_exp_in_EZH2_variant_or_WT_genome.png)

```R
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

TCGA_ESCA <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/TCGA_data/TCGA-ESCA_RNA_FPKM_expr.csv")
TCGA_ESCA$SYMBOL <- mapIds(x = org.Hs.eg.db,
                        keys = as.character(TCGA_ESCA$X),
            keytype ="ENSEMBL",
            column ="SYMBOL",
            multiVals="first")  
TCGA_ESCA <- na.omit(TCGA_ESCA)
TCGA_ESCA <- TCGA_ESCA[!duplicated(TCGA_ESCA$SYMBOL),]
rownames(TCGA_ESCA) <- TCGA_ESCA$SYMBOL
TCGA_ESCA_matrix <- TCGA_ESCA[,c(-1,-175)]
TCGA_ESCA_matrix <- data.frame(t(TCGA_ESCA_matrix))

TCGA_ESCA_matrix$malignant <- substring(rownames(TCGA_ESCA_matrix),14,15)
TCGA_ESCA_matrix_1 <- subset(TCGA_ESCA_matrix, malignant < 10)
TCGA_ESCA_matrix_1 <- TCGA_ESCA_matrix_1[,-162]
TCGA_ESCA_matrix_1$ID <- rownames(TCGA_ESCA_matrix_1)
TCGA_ESCA_matrix_1$ID <- substring(TCGA_ESCA_matrix_1$ID,1,12)
TCGA_ESCA_matrix_1$ID <- gsub("\\.","-",TCGA_ESCA_matrix_1$ID)

clinical_ESCA <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/TCGA_data/TCGA-ESCA_clinical.csv")
clinical_ESCC <- subset(clinical_ESCA, primary_diagnosis=="Basaloid squamous cell carcinoma" | primary_diagnosis=="Squamous cell carcinoma, keratinizing, NOS" | primary_diagnosis=="Squamous cell carcinoma, NOS")
clinical_ESCC_1 <- clinical_ESCC[,c("bcr_patient_barcode","primary_diagnosis")]
all <- merge(TCGA_ESCA_matrix_1,clinical_ESCC_1,by.x="ID",by.y="bcr_patient_barcode")
all <- all[!duplicated(all$ID),]
rownames(all) <- all$ID
ESCC <- all[,c(-1,-35083,-35084)]

clinical_ESCA <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/TCGA_data/TCGA-ESCA_clinical.csv")
clinical_Adno <- subset(clinical_ESCA, primary_diagnosis=="Adenocarcinoma, NOS" | primary_diagnosis=="Tubular adenocarcinoma" | primary_diagnosis=="Mucinous adenocarcinoma")
clinical_Adno_1 <- clinical_Adno[,c("bcr_patient_barcode","primary_diagnosis")]
all <- merge(TCGA_ESCA_matrix_1,clinical_Adno_1,by.x="ID",by.y="bcr_patient_barcode")
all <- all[!duplicated(all$ID),]
rownames(all) <- all$ID
ESCA <- all[,c(-1,-25538,-25539)]

genes <- na.omit(intersect(colnames(ESCC),colnames(ESCA)))
ESCA <- ESCA[,genes]
ESCC <- ESCC[,genes]
TFAP2C_ESCA <- data.frame(ESCA[,c("TFAP2C")])
TFAP2C_ESCC <- data.frame(ESCC[,c("TFAP2C")])
rownames(TFAP2C_ESCA) <- rownames(ESCA)
rownames(TFAP2C_ESCC) <- rownames(ESCC)
names(TFAP2C_ESCA) <- "TFAP2C"
names(TFAP2C_ESCC) <- "TFAP2C"
TFAP2C_ESCA$group <- "ESCA"
TFAP2C_ESCC$group <- "ESCC"
tmp_all <- rbind(TFAP2C_ESCA,TFAP2C_ESCC)
tmp_all$log_TFAP2C <- log(tmp_all$TFAP2C+1,2)

library(ggpubr)
library(ggplot2)
my_comparisons <- list(c("ESCA", "ESCC"))
ff <- ggboxplot(tmp_all, x = "group", y = "log_TFAP2C",add = "jitter",
               color = "group",palette = c("#08519c","#e31a1c")) +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="TFAP2C_IN_ESCC_AND_ESCA.png",width=3,height=6.4)

```

![TFAP2C_IN_ESCC_AND_ESCA](C:\Users\Irene\Desktop\NC_markdown_整理脚本\assets\TFAP2C_IN_ESCC_AND_ESCA.png)

```shell
java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185hi_vs_186low.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185hi_vs_186low.cls#TFAP2C_hi_versus_TFAP2C_low \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185c5 -gui false


java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185hi_vs_186low.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185hi_vs_186low.cls#TFAP2C_hi_versus_TFAP2C_low \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185c2 -gui false


java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185hi_vs_186low.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185hi_vs_186low.cls#TFAP2C_hi_versus_TFAP2C_low \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/TFAP2C_Clinical/TFAP2C_185h -gui false
```

