# TCGA-STAD analysis

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

STAD_FPKM_tumor <- read.csv("/mnt/data/user_data/abao/1_project/0_GSCC_NatureCommunications/TCGA_STAD/STAD_FPKM_tumor.csv")
ESCA_squa_sig  <- read.csv("/mnt/data/user_data/abao/1_project/0_GSCC_NatureCommunications/TCGA_STAD/1_TCGA_ESCA_suqa_vs_other_res_pj_mouse.csv")
ESCA_sig <- subset(ESCA_squa_sig,padj < 0.05 & log2FoldChange > 1)
ESCA_sig <- ESCA_sig[order(-ESCA_sig$log2FoldChange),]

squamous_sig <- data.frame(ESCA_sig$SYMBOL)
names(squamous_sig) <- "SYMBOL"
squamous_Exp <- merge(STAD_FPKM_tumor,squamous_sig,by="SYMBOL")
squamous_Exp <- squamous_Exp[!duplicated(squamous_Exp$SYMBOL),]
rownames(squamous_Exp) <- squamous_Exp$SYMBOL
squamous_Exp <- squamous_Exp[,c(-1,-2,-3)]
squamous_Exp_1 <- data.frame(apply(squamous_Exp,2,mean)) 
squamous_Exp <- data.frame(squamous_Exp_1)
names(squamous_Exp) <- "squamous_exp"

EZH2_Exp <- subset(STAD_FPKM_tumor,SYMBOL=="EZH2")
rownames(EZH2_Exp) <- EZH2_Exp$SYMBOL
EZH2_Exp_1 <- EZH2_Exp[,c(-1,-2)]
EZH2_Exp <- data.frame(t(EZH2_Exp_1[,c(-376)]))

tmp <- cbind(squamous_Exp,EZH2_Exp)
tmp_log2 <- log(tmp+1,2)


STAD_clinical <- read.csv("/mnt/data/user_data/abao/1_project/0_GSCC_NatureCommunications/TCGA_STAD/TCGA-STAD_clinical.csv")
tmp_log2$Tumor_Sample_Barcode <- rownames(tmp_log2) 
tmp_log2$Tumor_Sample_Barcode <- substring(tmp_log2$Tumor_Sample_Barcode,1,12)
tmp_log2_1 <- merge(tmp_log2,STAD_clinical,by="Tumor_Sample_Barcode")
meta <- tmp_log2_1
meta[is.na(meta$days_to_last_follow_up),]$days_to_last_follow_up <- "HHH"
tmp <- subset(meta,days_to_last_follow_up=="HHH")
tmp$days_to_last_follow_up <- tmp$days_to_death
no_na <- meta[setdiff(rownames(meta),rownames(tmp)),]
all_merge <- rbind(tmp,no_na)
all_merge <- subset(all_merge,days_to_last_follow_up != "HHH")

#all_merge <- meta
all_merge$vital_status <- as.character(all_merge$vital_status)
all_merge$status <- ifelse(all_merge$vital_status=="Alive",0,1)
all_merge$days_to_last_follow_up <- as.numeric(all_merge$days_to_last_follow_up)
all_merge <- all_merge[order(-all_merge$squamous_exp),]

library("survival")
library("survminer")
all_merge.cut <- surv_cutpoint(
   all_merge,
   time = "days_to_last_follow_up",
   event = "status",
   variables = c("squamous_exp"),
   progressbar=TRUE,
   minprop=0.05)
summary(all_merge.cut)
all_merge.cut.cat <- surv_categorize(all_merge.cut) 
library(survival)
fit <- survfit(Surv(days_to_last_follow_up, status) ~ squamous_exp, data = all_merge.cut.cat)
pdf("Fig_1f_Survival_curve_final.pdf")
ggsurvplot(fit, data = all_merge.cut.cat,
surv.median.line = "hv",
pval = TRUE,
ggtheme = theme_bw(),
risk.table=TRUE)
dev.off()
```

![squamous_exp_Survival_curve_final_页面_2](assets/squamous_exp_Survival_curve_final_页面_2.png)

```R
ff2 <- ggplot(all_merge, aes(y=EZH2, x= squamous_exp))+geom_point() + stat_smooth(method=lm)+stat_cor(data=all_merge, method = "pearson")
ggsave(ff2,file="0426_EZH2_squamous_correlation_final.png",width=6,height=6.4)
```

![Fig_1g_EZH2_squamous_correlation_final](assets/Fig_1g_EZH2_squamous_correlation_final.png)

```R

Patient_Group <- tmp_1[,c("Tumor_Sample_Barcode","Group")]
STAD_FPKM_tumor <- read.csv("/mnt/data/user_data/abao/1_project/0_GSCC_NatureCommunications/TCGA_STAD/STAD_FPKM_tumor.csv")
rownames(STAD_FPKM_tumor) <- STAD_FPKM_tumor$ENSG
STAD_FPKM_tumor <- STAD_FPKM_tumor[,c(-1,-2)]
STAD_FPKM_tumor <- STAD_FPKM_tumor[,c(-376)]
STAD_FPKM <- data.frame(t(STAD_FPKM_tumor))

STAD_FPKM$Tumor_Sample_Barcode <- substring(rownames(STAD_FPKM),1,12)

STAD_FPKM_1 <- merge(STAD_FPKM,Patient_Group,by="Tumor_Sample_Barcode")

Squamous_high <- subset(STAD_FPKM_1,Group=="Squamous_hi")
rownames(Squamous_high) <- Squamous_high$Tumor_Sample_Barcode
Squamous_high <- Squamous_high[,c(-1,-56514)]
Squamous_high_1 <- data.frame(t(Squamous_high))

Squamous_low <- subset(STAD_FPKM_1,Group=="Squamous_low")
rownames(Squamous_low) <- Squamous_low$Tumor_Sample_Barcode
Squamous_low <- Squamous_low[,c(-1,-56514)]
Squamous_low_1 <- data.frame(t(Squamous_low))


# calculate Squa hi vs low DEGs and logFoldChanges

logFC <- log2(rowMeans(as.matrix(Squamous_high_1)+1) / rowMeans(as.matrix(Squamous_low_1)+1))

library("future.apply")
p_values <- future_lapply(seq(1,nrow(Squamous_high_1)), function(x){
  res <- wilcox.test(x = t(Squamous_high_1[x,])[,1], y = t(Squamous_low_1[x,])[,1] )
  res$p.value
})
p <- unlist(p_values)
p.adj <- p.adjust(p, method = "fdr")
genelist <- as.data.frame(logFC)
genelist$p_values <- p
genelist$p_adj <- p.adj
genelist$gene <- rownames(genelist)

aa <- cbind(Squamous_high_1,Squamous_low_1)
basemean <- apply(aa,1,mean)
all <- cbind(basemean,genelist)

all$SYMBOL <- mapIds(x = org.Hs.eg.db,
                        keys = rownames(all),
                  keytype ="ENSEMBL",
                  column ="SYMBOL",
                  multiVals="first")
write.csv(all,"TCGA_STAD_Squamous_58hi_vs_313low_Wilcox_test.csv")



********************************************************************************************
********************************volcano_plot************************************************
********************************************************************************************
********************************************************************************************

library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")

STAD_tpm <- read.csv("/mnt/data/user_data/abao/1_project/0_GSCC_NatureCommunications/0_zhenglijiaoben/TCGA_STAD_Squamous_58hi_vs_313low_Wilcox_test.csv")
genes_plant <- c("MUC17","MUC13","LGR5","TFF1","ASCL2","MUC3A","MUC1","MKI67","ERBB2","EZH2","KRT5","KRT6A","KRT6C","KRT14","KRT6B","KRT24")
STAD_tpm <- na.omit(STAD_tpm)
library(ggplot2)

logFC <-STAD_tpm$logFC
pvalue <- STAD_tpm$p_values
data <- data.frame(logFC=logFC,pvalue=pvalue)
data$sig[(data$pvalue > 0.05|data$pvalue=="NA" | -1 < data$pvalue |  data$pvalue < -1 )] <- "no"

data$sig[data$pvalue <= 0.05 & data$logFC >= 0.5] <- "up"
data$sig[data$pvalue <= 0.05& data$logFC <= -0.5] <- "down"
data$symbol <- STAD_tpm$SYMBOL
data$log10_pvalue <- -log10(data$pvalue)

range(data$log10_pvalue)
data$log10_pvalue[data$log10_pvalue > 22] <- 22

data$logFC[data$logFC > 5] <- 5
data$logFC[data$logFC < -5] <- -5

sel_genes <- function(gene){
  tmp <- subset(data,symbol==gene)
  return(tmp)
}
aa <- future_lapply(as.list(as.character(genes_plant)),sel_genes)
all_res_plant <- do.call(rbind,aa)

ff <- ggplot(data, aes(x = logFC, y = log10_pvalue, color = sig)) +
geom_point(alpha = 0.6, size = 1) +

scale_colour_manual(values  = c('#FF3333','#0066CC','grey'), limits = c('up', 'down', 'no')) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.3) +
geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
xlim(-5, 5) + ylim(0, 22) +
labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)\n', color = '', title = '\n')+ theme(legend.position = 'right') +
geom_text_repel(data = all_res_plant, aes(x = logFC, y = -log10(pvalue), label = symbol),
        size = 5,fontface="bold", 
    color="grey50", box.padding=unit(0.35, "lines"), 
    point.padding=unit(0.5, "lines"), segment.colour = "grey50")
ggsave(ff,file="FigS1e_volcono.png",width=6,height=5.7)


```

![FigS1e_volcono](assets/FigS1e_volcono.png)

```R
library(ggpubr)
library(ggplot2)
rownames(tmp_2) <- tmp_2$Tumor_Sample_Barcode
my_comparisons <- list(c("squamous_hi", "squamous_low"))
tmp_3 <- tmp_2[,c("Group_squamous","EZH2")]
ff1 <- ggboxplot(tmp_3, x ="Group_squamous", y = "EZH2",add = "jitter",
               color ="Group_squamous",palette = c("#08519c","#e31a1c")) +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")

ggsave(ff1,file="Fig_1d_EZH2_squamous_group.png",width=5,height=6.4)
```

![Fig_1d_EZH2_squamous_group](assets/Fig_1d_EZH2_squamous_group.png)

```shell
java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low.cls#Squamous_hi_versus_Squamous_low \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low_c5 -gui false


java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low.cls#Squamous_hi_versus_Squamous_low \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low_c2 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low.cls#Squamous_hi_versus_Squamous_low \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_human_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Clinical_data/Squamous_Clinical/Squamous_58hi_vs_313low_h -gui false
```

