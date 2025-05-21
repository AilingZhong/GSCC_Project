# Data visualization of Figure 2

```R
library(Seurat)
library(ggplot2)
library(pheatmap)
require(RColorBrewer)
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")

sgEZH2_VS_sgScr <- read.csv("/mnt/data2/userdata/abao/Published_Project_code/GSCC_code/RNAseq_data/sgEZH2_VS_sgScr_DEseq_Results.csv")
sig_Normal_RNA <- subset(sgEZH2_VS_sgScr,pvalue <=0.05 & abs(log2FoldChange) > 0.5 )

aa <- subset(sig_Normal_RNA, log2FoldChange > 0.5)
bb <- subset(sig_Normal_RNA, log2FoldChange < -0.5)

genes_sel <- c("Tfap2c","Krt5","Top2a","Krt6b","Krt6a","Krt14","Krt12","Krt76","Twist1","Twist2","Hoxb4","Hoxc6","Hoxc8","Hoxc5","Hoxb2")
rownames(sig_Normal_RNA) <- sig_Normal_RNA$X
sig_Normal_RNA <- sig_Normal_RNA[order(-sig_Normal_RNA$log2FoldChange),]
sig_Normal_RNA <- sig_Normal_RNA[,c(2:10)]
sig_Normal_RNA <- na.omit(sig_Normal_RNA)

tmp_data <- log(sig_Normal_RNA+1,2)
chonglai_zscore_1 <- t(apply(tmp_data, 1, function(x) (x-mean(x))/sd(x)))
range(chonglai_zscore_1)
chonglai_zscore_1[chonglai_zscore_1 > 1] <- 1
chonglai_zscore_1[chonglai_zscore_1 < -1] <- -1
SeuratObject <- CreateSeuratObject(counts = chonglai_zscore_1, project = "TM")
gene <- rownames(SeuratObject)
SeuratObject@meta.data$group <- rownames(SeuratObject@meta.data)
pdf("Fig2d_Normal_heatmap_label.pdf",width=8,height=9)
adjusted_heatmap(seurat_obj=SeuratObject,group="orig.ident",gene = gene,all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),
  min_and_max_cut=2,show_row_names=FALSE,mark_gene=genes_sel,label_size=0,scale = FALSE)
dev.off()
```

![Fig2d_Normal_heatmap_label-01](C:\Users\Irene\Desktop\NC_markdown_整理脚本\assets\Fig2d_Normal_heatmap_label-01.png)

```shell
java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal.cls#sgEzh2_versus_sgScr \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal_c5 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal.cls#sgEzh2_versus_sgScr \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal_c2 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea-3.0.jar -Xmx10240m xtools.gsea.Gsea \
-res /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal.gct \
-cls /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal.cls#sgEzh2_versus_sgScr \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 200 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out /mnt/data/userdata/abao/project/1_bulk_sequence/13_ZhangMengsha/Normal_Organoids/bam_bai_files/Useful_data/Normal_Orgnoids_workfile/Normal_h -gui false
```

