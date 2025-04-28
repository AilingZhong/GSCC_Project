# Targeted bisulfite sequencing analysis

```shell
/mnt/data/user_data/abao/4_packages/fastp -i P1_Tumor_R1.fastq.gz -I P1_Tumor_R2.fastq.gz -o P1_Tumor_trimmed_R1.fastq -O P1_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P4_Tumor_R1.fastq.gz -I P4_Tumor_R2.fastq.gz -o P4_Tumor_trimmed_R1.fastq -O P4_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P4_Normal_R1.fastq.gz -I P4_Normal_R2.fastq.gz -o P4_Normal_trimmed_R1.fastq -O P4_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P5_Normal_R1.fastq.gz -I P5_Normal_R2.fastq.gz -o P5_Normal_trimmed_R1.fastq -O P5_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P5_Tumor_R1.fastq.gz -I P5_Tumor_R2.fastq.gz -o P5_Tumor_trimmed_R1.fastq -O P5_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P8_Normal_R1.fastq.gz -I P8_Normal_R2.fastq.gz -o P8_Normal_trimmed_R1.fastq -O P8_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P8_Tumor_R1.fastq.gz -I P8_Tumor_R2.fastq.gz -o P8_Tumor_trimmed_R1.fastq -O P8_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P9_Tumor_R1.fastq.gz -I P9_Tumor_R2.fastq.gz -o P9_Tumor_trimmed_R1.fastq -O P9_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P13_Tumor_R1.fastq.gz -I P13_Tumor_R2.fastq.gz -o P13_Tumor_trimmed_R1.fastq -O P13_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P13_Normal_R1.fastq.gz -I P13_Normal_R2.fastq.gz -o P13_Normal_trimmed_R1.fastq -O P13_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P14_Tumor_R1.fastq.gz -I P14_Tumor_R2.fastq.gz -o P14_Tumor_trimmed_R1.fastq -O P14_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P14_Normal_R1.fastq.gz -I P14_Normal_R2.fastq.gz -o P14_Normal_trimmed_R1.fastq -O P14_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P15_Normal_R1.fastq.gz -I P15_Normal_R2.fastq.gz -o P15_Normal_trimmed_R1.fastq -O P15_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P15_Tumor_R1.fastq.gz -I P15_Tumor_R2.fastq.gz -o P15_Tumor_trimmed_R1.fastq -O P15_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P17_Normal_R1.fastq.gz -I P17_Normal_R2.fastq.gz -o P17_Normal_trimmed_R1.fastq -O P17_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P18_Normal_R1.fastq.gz -I P18_Normal_R2.fastq.gz -o P18_Normal_trimmed_R1.fastq -O P18_Normal_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P18_Tumor_R1.fastq.gz -I P18_Tumor_R2.fastq.gz -o P18_Tumor_trimmed_R1.fastq -O P18_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P19_Tumor_R1.fastq.gz -I P19_Tumor_R2.fastq.gz -o P19_Tumor_trimmed_R1.fastq -O P19_Tumor_trimmed_R2.fastq
/mnt/data/user_data/abao/4_packages/fastp -i P20_Tumor_R1.fastq.gz -I P20_Tumor_R2.fastq.gz -o P20_Tumor_trimmed_R1.fastq -O P20_Tumor_trimmed_R2.fastq


#EZH2
Chromosome 7: 148,883,700-148,885,800

samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa 7:148883700-148885800 > EZH2_promoter_region.fa


#build index
/mnt/data2/userdata/abao/programm/Bismark-0.24.2/bismark_genome_preparation \
/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/EZH2_region/



#align
vim TBS_sample.csv
P20_Tumor
P19_Tumor
P18_Tumor
P18_Normal
P17_Normal
P15_Tumor
P15_Normal
P14_Tumor
P14_Normal
P13_Tumor
P13_Normal
P9_Tumor
P8_Tumor
P8_Normal
P5_Tumor
P5_Normal
P4_Tumor
P4_Normal
P1_Tumor


cat TBS_sample.csv | while read id ; do
arr=($id)
sample=${arr[0]}
fq2=${arr[0]}'_trimmed_R2.fastq'
fq1=${arr[0]}'_trimmed_R1.fastq'
echo ${fq2}
echo ${fq1}
echo $sample
/mnt/data2/userdata/abao/programm/Bismark-0.24.2/bismark \
/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/EZH2_region/ \
-1 ./trimmed_fq/$fq1 -2 ./trimmed_fq/$fq2 \
-o ./align_results/${sample}.results;
done


samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P1_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P1_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P4_Normal_trimmed_R1_bismark_bt2_pe.bam -o P4_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P4_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P4_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P5_Normal_trimmed_R1_bismark_bt2_pe.bam -o P5_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P5_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P5_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P8_Normal_trimmed_R1_bismark_bt2_pe.bam -o P8_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P8_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P8_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P9_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P9_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P13_Normal_trimmed_R1_bismark_bt2_pe.bam -o P13_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P13_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P13_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P14_Normal_trimmed_R1_bismark_bt2_pe.bam -o P14_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P14_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P14_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P15_Normal_trimmed_R1_bismark_bt2_pe.bam -o P15_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P15_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P15_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P17_Normal_trimmed_R1_bismark_bt2_pe.bam -o P17_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P18_Normal_trimmed_R1_bismark_bt2_pe.bam -o P18_Normal.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P18_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P18_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P19_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P19_Tumor.sorted.bam
samtools sort /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/P20_Tumor_trimmed_R1_bismark_bt2_pe.bam -o P20_Tumor.sorted.bam


cat TBS_sample.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
cd /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/methy_extractor/
/mnt/data2/userdata/abao/programm/Bismark-0.24.2/bismark_methylation_extractor --comprehensive --no_overlap --bedGraph \
--genome_folder /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/EZH2_region/ --cytosine_report \
-p /mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/fastp/align_results/${sample}_trimmed_R1_bismark_bt2_pe.bam ;
done
```

```R

P18_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P18_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P18_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P18_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")

P15_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P15_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P15_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P15_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")

P14_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P14_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P14_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P14_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")

P13_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P13_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P13_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P13_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")

P8_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P8_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P8_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P8_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")

P5_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P5_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P5_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P5_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")

P4_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P4_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P4_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P4_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")



library(rtracklayer)
library(data.table)
library(patchwork)
P13_Normal_EZH2 <- subset(P13_Normal,V1=="chr:7")
P13_Normal_EZH2_freq <- data.frame(position=148883700+P13_Normal_EZH2$V2,CpG_methy_ratio=P13_Normal_EZH2$V4/(P13_Normal_EZH2$V4+P13_Normal_EZH2$V5))
P13_Normal_EZH2_freq[is.na(P13_Normal_EZH2_freq)] <- 0
P13_Normal_EZH2_freq <- P13_Normal_EZH2_freq[order(P13_Normal_EZH2_freq$position),] 
P13_Tumor_EZH2 <- subset(P13_Tumor,V1=="chr:7")
P13_Tumor_EZH2_freq <- data.frame(position=148883700+P13_Tumor_EZH2$V2,CpG_methy_ratio=P13_Tumor_EZH2$V4/(P13_Tumor_EZH2$V4+P13_Tumor_EZH2$V5))
P13_Tumor_EZH2_freq[is.na(P13_Tumor_EZH2_freq)] <- 0
P13_Tumor_EZH2_freq <- P13_Tumor_EZH2_freq[order(P13_Tumor_EZH2_freq$position),] 
P13_data <- data.frame(position=P13_Normal_EZH2_freq$position,N_ratio=P13_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P13_Tumor_EZH2_freq$CpG_methy_ratio, P13_T_vs_N=P13_Tumor_EZH2_freq$CpG_methy_ratio-P13_Normal_EZH2_freq$CpG_methy_ratio)
library(dplyr)
P13_data <- P13_data %>%
  mutate(
    Normal_me_state = ifelse(N_ratio == 0 & T_ratio == 0, 0, ifelse(P13_T_vs_N < 0, 1, 0)),
    Tumor_me_state = ifelse(N_ratio == 0 & T_ratio == 0, 0, ifelse(P13_T_vs_N > 0, 1, 0))
  )

peak_1 <- subset(P13_data, position >=  148883786  & position <= 148883935)
peak_1 <- peak_1 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_1$Normal_me_state <- ifelse(peak_1$N_ratio < 0.025, "0","1")
peak_1$Tumor_me_state <- ifelse(peak_1$T_ratio < 0.025, "0","1")
P13_Normal_peak_1 <- ggplot(peak_1, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P13_Tumor_peak_1 <- ggplot(peak_1, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P13_Normal_peak_1+P13_Tumor_peak_1
ggsave(file="peak_1_P13_Normal_and_Tumor.pdf", ff, width = 9, height = 4)


peak_2 <- subset(P13_data, position >=  148884063  & position <= 148884256)
peak_2 <- peak_2 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_2$Normal_me_state <- ifelse(peak_2$N_ratio < 0.025, "0","1")
peak_2$Tumor_me_state <- ifelse(peak_2$T_ratio < 0.025, "0","1")
P13_Normal_peak_2 <- ggplot(peak_2, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P13_Tumor_peak_2 <- ggplot(peak_2, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P13_Normal_peak_2+P13_Tumor_peak_2
ggsave(file="peak_2_P13_Normal_and_Tumor.pdf", ff, width = 6, height = 3)



peak_3 <- subset(P13_data, position >=  148884749  & position <= 148884918)
peak_3 <- peak_3 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_3$Normal_me_state <- ifelse(peak_3$N_ratio < 0.025, "0","1")
peak_3$Tumor_me_state <- ifelse(peak_3$T_ratio < 0.025, "0","1")
P13_Normal_peak_3 <- ggplot(peak_3, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P13_Tumor_peak_3 <- ggplot(peak_3, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P13_Normal_peak_3+P13_Tumor_peak_3
ggsave("peak_3_P13_Normal_and_Tumor.pdf", ff, width = 8, height = 4)


peak_4 <- subset(P13_data, position >=  148885601  & position <= 148885782)
peak_4 <- peak_4 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_4$Normal_me_state <- ifelse(peak_4$N_ratio < 0.025, "0","1")
peak_4$Tumor_me_state <- ifelse(peak_4$T_ratio < 0.025, "0","1")
P13_Normal_peak_4 <- ggplot(peak_4, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P13_Tumor_peak_4 <- ggplot(peak_4, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  # 调整行列间距
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P13_Normal_peak_4+P13_Tumor_peak_4
ggsave("peak_4_P13_Normal_and_Tumor_test.pdf", ff, width = 8, height = 4)



P15_Normal_EZH2 <- subset(P15_Normal,V1=="chr:7")
P15_Normal_EZH2_freq <- data.frame(position=148883700+P15_Normal_EZH2$V2,CpG_methy_ratio=P15_Normal_EZH2$V4/(P15_Normal_EZH2$V4+P15_Normal_EZH2$V5))
P15_Normal_EZH2_freq[is.na(P15_Normal_EZH2_freq)] <- 0
P15_Normal_EZH2_freq <- P15_Normal_EZH2_freq[order(P15_Normal_EZH2_freq$position),] 
P15_Tumor_EZH2 <- subset(P15_Tumor,V1=="chr:7")
P15_Tumor_EZH2_freq <- data.frame(position=148883700+P15_Tumor_EZH2$V2,CpG_methy_ratio=P15_Tumor_EZH2$V4/(P15_Tumor_EZH2$V4+P15_Tumor_EZH2$V5))
P15_Tumor_EZH2_freq[is.na(P15_Tumor_EZH2_freq)] <- 0
P15_Tumor_EZH2_freq <- P15_Tumor_EZH2_freq[order(P15_Tumor_EZH2_freq$position),] 
P15_data <- data.frame(position=P15_Normal_EZH2_freq$position,N_ratio=P15_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P15_Tumor_EZH2_freq$CpG_methy_ratio, P15_T_vs_N=P15_Tumor_EZH2_freq$CpG_methy_ratio-P15_Normal_EZH2_freq$CpG_methy_ratio)
library(dplyr)
P15_data <- P15_data %>%
  mutate(
    Normal_me_state = ifelse(N_ratio == 0 & T_ratio == 0, 0, ifelse(P15_T_vs_N < 0, 1, 0)),
    Tumor_me_state = ifelse(N_ratio == 0 & T_ratio == 0, 0, ifelse(P15_T_vs_N > 0, 1, 0))
  )

peak_1 <- subset(P15_data, position >=  148883786  & position <= 148883935)
peak_1 <- peak_1 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_1$Normal_me_state <- ifelse(peak_1$N_ratio < 0.025, "0","1")
peak_1$Tumor_me_state <- ifelse(peak_1$T_ratio < 0.025, "0","1")
P15_Normal_peak_1 <- ggplot(peak_1, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P15_Tumor_peak_1 <- ggplot(peak_1, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P15_Normal_peak_1+P15_Tumor_peak_1
ggsave(file="peak_1_P15_Normal_and_Tumor.pdf", ff, width = 9, height = 4)


peak_2 <- subset(P15_data, position >=  148884063  & position <= 148884256)
peak_2 <- peak_2 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_2$Normal_me_state <- ifelse(peak_2$N_ratio < 0.025, "0","1")
peak_2$Tumor_me_state <- ifelse(peak_2$T_ratio < 0.025, "0","1")
P15_Normal_peak_2 <- ggplot(peak_2, aes(x = col_pos * 0.7, y = -row_group * 0.8)) +  
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P15_Tumor_peak_2 <- ggplot(peak_2, aes(x = col_pos * 0.7, y = -row_group * 0.8)) + 
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P15_Normal_peak_2+P15_Tumor_peak_2
ggsave(file="peak_2_P15_Normal_and_Tumor.pdf", ff, width = 6, height = 3)



peak_3 <- subset(P15_data, position >=  148884749  & position <= 148884918)
peak_3 <- peak_3 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_3$Normal_me_state <- ifelse(peak_3$N_ratio < 0.025, "0","1")
peak_3$Tumor_me_state <- ifelse(peak_3$T_ratio < 0.025, "0","1")
P15_Normal_peak_3 <- ggplot(peak_3, aes(x = col_pos * 0.7, y = -row_group * 0.8)) + 
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P15_Tumor_peak_3 <- ggplot(peak_3, aes(x = col_pos * 0.7, y = -row_group * 0.8)) + 
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P15_Normal_peak_3+P15_Tumor_peak_3
ggsave("peak_3_P15_Normal_and_Tumor.pdf", ff, width = 8, height = 4)



peak_4 <- subset(P15_data, position >=  148885601  & position <= 148885782)
peak_4 <- peak_4 %>%
  mutate(row_group = (row_number() - 1) %/% 8,   
         col_pos = (row_number() - 1) %% 8)      
peak_4$Normal_me_state <- ifelse(peak_4$N_ratio < 0.025, "0","1")
peak_4$Tumor_me_state <- ifelse(peak_4$T_ratio < 0.025, "0","1")
P15_Normal_peak_4 <- ggplot(peak_4, aes(x = col_pos * 0.7, y = -row_group * 0.8)) + 
  geom_point(aes(fill = factor(Normal_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
P15_Tumor_peak_4 <- ggplot(peak_4, aes(x = col_pos * 0.7, y = -row_group * 0.8)) + 
  geom_point(aes(fill = factor(Tumor_me_state)), shape = 21, size = 3, color = "black") +  
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +  
  theme_minimal() +
  labs(x = "Position in Group", y = "Read", fill = "Methylation Status") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1.0)  # 适当调整比例，保持点之间间距
ff <- P15_Normal_peak_4+P15_Tumor_peak_4
ggsave("peak_4_P15_Normal_and_Tumor_test.pdf", ff, width = 8, height = 4)

```

```R
#EZH2: Chromosome 7: 148,883,700-148,885,800


P18_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P18_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P18_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P18_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P15_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P15_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P15_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P15_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P14_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P14_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P14_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P14_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P13_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P13_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P13_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P13_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P8_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P8_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P8_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P8_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P5_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P5_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P5_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P5_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P4_Tumor <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P4_Tumor_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")
P4_Normal <- read.table("/mnt/data2/userdata/abao/ZhangMengsha/Target_Methylation_GSCC/methylation_extractor/P4_Normal_trimmed_R1_bismark_bt2_pe.CpG_report.txt",head=FALSE,sep="\t")



P18_Normal_EZH2 <- subset(P18_Normal,V1=="chr:7")
P18_Tumor_EZH2 <- subset(P18_Tumor,V1=="chr:7")
library(rtracklayer)
library(data.table)
P18_Normal_EZH2_freq <- data.frame(position=148883700+P18_Normal_EZH2$V2,CpG_methy_ratio=P18_Normal_EZH2$V4/(P18_Normal_EZH2$V4+P18_Normal_EZH2$V5))
P18_Normal_EZH2_freq[is.na(P18_Normal_EZH2_freq)] <- 0
P18_Normal_EZH2_freq <- P18_Normal_EZH2_freq[order(P18_Normal_EZH2_freq$position),] 
P18_Tumor_EZH2_freq <- data.frame(position=148883700+P18_Tumor_EZH2$V2,CpG_methy_ratio=P18_Tumor_EZH2$V4/(P18_Tumor_EZH2$V4+P18_Tumor_EZH2$V5))
P18_Tumor_EZH2_freq[is.na(P18_Tumor_EZH2_freq)] <- 0
P18_Tumor_EZH2_freq <- P18_Tumor_EZH2_freq[order(P18_Tumor_EZH2_freq$position),] 
P18_data <- data.frame(position=P18_Normal_EZH2_freq$position,N_ratio=P18_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P18_Tumor_EZH2_freq$CpG_methy_ratio, P18_T_vs_N=P18_Tumor_EZH2_freq$CpG_methy_ratio-P18_Normal_EZH2_freq$CpG_methy_ratio)
P18_tmp <- subset(P18_data,P18_T_vs_N > 0 | P18_T_vs_N )
P18_tmp <- P18_tmp[order(P18_tmp$position),] 
list <- data.frame(chr7_pos=c(148883700:148885800),order=1:2101)
P18_tmp <- merge(list,P18_tmp,by.x="chr7_pos",by.y="position",all.x=TRUE)
P18_tmp[is.na(P18_tmp)] <- 0
P18_Normal_EZH2 <- data.frame(P18_tmp[,c("N_ratio")])
P18_Tumor_EZH2 <- data.frame(P18_tmp[,c("T_ratio")])
names(P18_Normal_EZH2) <- "Methy"
names(P18_Tumor_EZH2) <- "Methy"
P18_Normal_EZH2$Group <- "Normal"
P18_Tumor_EZH2$Group <- "Tumor"
rownames(P18_Normal_EZH2) <- P18_tmp$chr7_pos
rownames(P18_Tumor_EZH2) <- P18_tmp$chr7_pos
P18_Normal_EZH2$SampleID <- P18_tmp$chr7_pos
P18_Tumor_EZH2$SampleID <- P18_tmp$chr7_pos
P18_Normal_EZH2 <- P18_Normal_EZH2[order(P18_Normal_EZH2$SampleID),]
P18_Tumor_EZH2 <- P18_Tumor_EZH2[order(P18_Tumor_EZH2$SampleID),]
ff_P18_Normal <- ggplot(P18_Normal_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 
ff_P18_Tumor <- ggplot(P18_Tumor_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 


P15_Normal_EZH2 <- subset(P15_Normal,V1=="chr:7")
P15_Tumor_EZH2 <- subset(P15_Tumor,V1=="chr:7")
library(rtracklayer)
library(data.table)
P15_Normal_EZH2_freq <- data.frame(position=148883700+P15_Normal_EZH2$V2,CpG_methy_ratio=P15_Normal_EZH2$V4/(P15_Normal_EZH2$V4+P15_Normal_EZH2$V5))
P15_Normal_EZH2_freq[is.na(P15_Normal_EZH2_freq)] <- 0
P15_Normal_EZH2_freq <- P15_Normal_EZH2_freq[order(P15_Normal_EZH2_freq$position),] 
P15_Tumor_EZH2_freq <- data.frame(position=148883700+P15_Tumor_EZH2$V2,CpG_methy_ratio=P15_Tumor_EZH2$V4/(P15_Tumor_EZH2$V4+P15_Tumor_EZH2$V5))
P15_Tumor_EZH2_freq[is.na(P15_Tumor_EZH2_freq)] <- 0
P15_Tumor_EZH2_freq <- P15_Tumor_EZH2_freq[order(P15_Tumor_EZH2_freq$position),] 
P15_data <- data.frame(position=P15_Normal_EZH2_freq$position,N_ratio=P15_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P15_Tumor_EZH2_freq$CpG_methy_ratio, P15_T_vs_N=P15_Tumor_EZH2_freq$CpG_methy_ratio-P15_Normal_EZH2_freq$CpG_methy_ratio)
P15_tmp <- subset(P15_data,P15_T_vs_N > 0 | P15_T_vs_N )
P15_tmp <- P15_tmp[order(P15_tmp$position),] 
list <- data.frame(chr7_pos=c(148883700:148885800),order=1:2101)
P15_tmp <- merge(list,P15_tmp,by.x="chr7_pos",by.y="position",all.x=TRUE)
P15_tmp[is.na(P15_tmp)] <- 0
P15_Normal_EZH2 <- data.frame(P15_tmp[,c("N_ratio")])
P15_Tumor_EZH2 <- data.frame(P15_tmp[,c("T_ratio")])
names(P15_Normal_EZH2) <- "Methy"
names(P15_Tumor_EZH2) <- "Methy"
P15_Normal_EZH2$Group <- "Normal"
P15_Tumor_EZH2$Group <- "Tumor"
rownames(P15_Normal_EZH2) <- P15_tmp$chr7_pos
rownames(P15_Tumor_EZH2) <- P15_tmp$chr7_pos
P15_Normal_EZH2$SampleID <- P15_tmp$chr7_pos
P15_Tumor_EZH2$SampleID <- P15_tmp$chr7_pos
P15_Normal_EZH2 <- P15_Normal_EZH2[order(P15_Normal_EZH2$SampleID),]
P15_Tumor_EZH2 <- P15_Tumor_EZH2[order(P15_Tumor_EZH2$SampleID),]

ff_P15_Normal <- ggplot(P15_Normal_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 
ff_P15_Tumor <- ggplot(P15_Tumor_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 



P14_Normal_EZH2 <- subset(P14_Normal,V1=="chr:7")
P14_Tumor_EZH2 <- subset(P14_Tumor,V1=="chr:7")
library(rtracklayer)
library(data.table)
P14_Normal_EZH2_freq <- data.frame(position=148883700+P14_Normal_EZH2$V2,CpG_methy_ratio=P14_Normal_EZH2$V4/(P14_Normal_EZH2$V4+P14_Normal_EZH2$V5))
P14_Normal_EZH2_freq[is.na(P14_Normal_EZH2_freq)] <- 0
P14_Normal_EZH2_freq <- P14_Normal_EZH2_freq[order(P14_Normal_EZH2_freq$position),] 
P14_Tumor_EZH2_freq <- data.frame(position=148883700+P14_Tumor_EZH2$V2,CpG_methy_ratio=P14_Tumor_EZH2$V4/(P14_Tumor_EZH2$V4+P14_Tumor_EZH2$V5))
P14_Tumor_EZH2_freq[is.na(P14_Tumor_EZH2_freq)] <- 0
P14_Tumor_EZH2_freq <- P14_Tumor_EZH2_freq[order(P14_Tumor_EZH2_freq$position),] 
P14_data <- data.frame(position=P14_Normal_EZH2_freq$position,N_ratio=P14_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P14_Tumor_EZH2_freq$CpG_methy_ratio, P14_T_vs_N=P14_Tumor_EZH2_freq$CpG_methy_ratio-P14_Normal_EZH2_freq$CpG_methy_ratio)
P14_tmp <- subset(P14_data,P14_T_vs_N > 0 | P14_T_vs_N )
P14_tmp <- P14_tmp[order(P14_tmp$position),] 
list <- data.frame(chr7_pos=c(148883700:148885800),order=1:2101)
P14_tmp <- merge(list,P14_tmp,by.x="chr7_pos",by.y="position",all.x=TRUE)
P14_tmp[is.na(P14_tmp)] <- 0
P14_Normal_EZH2 <- data.frame(P14_tmp[,c("N_ratio")])
P14_Tumor_EZH2 <- data.frame(P14_tmp[,c("T_ratio")])
names(P14_Normal_EZH2) <- "Methy"
names(P14_Tumor_EZH2) <- "Methy"
P14_Normal_EZH2$Group <- "Normal"
P14_Tumor_EZH2$Group <- "Tumor"
rownames(P14_Normal_EZH2) <- P14_tmp$chr7_pos
rownames(P14_Tumor_EZH2) <- P14_tmp$chr7_pos
P14_Normal_EZH2$SampleID <- P14_tmp$chr7_pos
P14_Tumor_EZH2$SampleID <- P14_tmp$chr7_pos
P14_Normal_EZH2 <- P14_Normal_EZH2[order(P14_Normal_EZH2$SampleID),]
P14_Tumor_EZH2 <- P14_Tumor_EZH2[order(P14_Tumor_EZH2$SampleID),]
ff_P14_Normal <- ggplot(P14_Normal_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 
ff_P14_Tumor <- ggplot(P14_Tumor_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 






P13_Normal_EZH2 <- subset(P13_Normal,V1=="chr:7")
P13_Tumor_EZH2 <- subset(P13_Tumor,V1=="chr:7")
library(rtracklayer)
library(data.table)
P13_Normal_EZH2_freq <- data.frame(position=148883700+P13_Normal_EZH2$V2,CpG_methy_ratio=P13_Normal_EZH2$V4/(P13_Normal_EZH2$V4+P13_Normal_EZH2$V5))
P13_Normal_EZH2_freq[is.na(P13_Normal_EZH2_freq)] <- 0
P13_Normal_EZH2_freq <- P13_Normal_EZH2_freq[order(P13_Normal_EZH2_freq$position),] 
P13_Tumor_EZH2_freq <- data.frame(position=148883700+P13_Tumor_EZH2$V2,CpG_methy_ratio=P13_Tumor_EZH2$V4/(P13_Tumor_EZH2$V4+P13_Tumor_EZH2$V5))
P13_Tumor_EZH2_freq[is.na(P13_Tumor_EZH2_freq)] <- 0
P13_Tumor_EZH2_freq <- P13_Tumor_EZH2_freq[order(P13_Tumor_EZH2_freq$position),] 
P13_data <- data.frame(position=P13_Normal_EZH2_freq$position,N_ratio=P13_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P13_Tumor_EZH2_freq$CpG_methy_ratio, P13_T_vs_N=P13_Tumor_EZH2_freq$CpG_methy_ratio-P13_Normal_EZH2_freq$CpG_methy_ratio)
P13_tmp <- subset(P13_data,P13_T_vs_N > 0 | P13_T_vs_N )
P13_tmp <- P13_tmp[order(P13_tmp$position),] 
list <- data.frame(chr7_pos=c(148883700:148885800),order=1:2101)
P13_tmp <- merge(list,P13_tmp,by.x="chr7_pos",by.y="position",all.x=TRUE)
P13_tmp[is.na(P13_tmp)] <- 0
P13_Normal_EZH2 <- data.frame(P13_tmp[,c("N_ratio")])
P13_Tumor_EZH2 <- data.frame(P13_tmp[,c("T_ratio")])
names(P13_Normal_EZH2) <- "Methy"
names(P13_Tumor_EZH2) <- "Methy"
P13_Normal_EZH2$Group <- "Normal"
P13_Tumor_EZH2$Group <- "Tumor"
rownames(P13_Normal_EZH2) <- P13_tmp$chr7_pos
rownames(P13_Tumor_EZH2) <- P13_tmp$chr7_pos
P13_Normal_EZH2$SampleID <- P13_tmp$chr7_pos
P13_Tumor_EZH2$SampleID <- P13_tmp$chr7_pos
P13_Normal_EZH2 <- P13_Normal_EZH2[order(P13_Normal_EZH2$SampleID),]
P13_Tumor_EZH2 <- P13_Tumor_EZH2[order(P13_Tumor_EZH2$SampleID),]
ff_P13_Normal <- ggplot(P13_Normal_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 
ff_P13_Tumor <- ggplot(P13_Tumor_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 


P8_Normal_EZH2 <- subset(P8_Normal,V1=="chr:7")
P8_Tumor_EZH2 <- subset(P8_Tumor,V1=="chr:7")
library(rtracklayer)
library(data.table)
P8_Normal_EZH2_freq <- data.frame(position=148883700+P8_Normal_EZH2$V2,CpG_methy_ratio=P8_Normal_EZH2$V4/(P8_Normal_EZH2$V4+P8_Normal_EZH2$V5))
P8_Normal_EZH2_freq[is.na(P8_Normal_EZH2_freq)] <- 0
P8_Normal_EZH2_freq <- P8_Normal_EZH2_freq[order(P8_Normal_EZH2_freq$position),] 
P8_Tumor_EZH2_freq <- data.frame(position=148883700+P8_Tumor_EZH2$V2,CpG_methy_ratio=P8_Tumor_EZH2$V4/(P8_Tumor_EZH2$V4+P8_Tumor_EZH2$V5))
P8_Tumor_EZH2_freq[is.na(P8_Tumor_EZH2_freq)] <- 0
P8_Tumor_EZH2_freq <- P8_Tumor_EZH2_freq[order(P8_Tumor_EZH2_freq$position),] 
P8_data <- data.frame(position=P8_Normal_EZH2_freq$position,N_ratio=P8_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P8_Tumor_EZH2_freq$CpG_methy_ratio, P8_T_vs_N=P8_Tumor_EZH2_freq$CpG_methy_ratio-P8_Normal_EZH2_freq$CpG_methy_ratio)
P8_tmp <- subset(P8_data,P8_T_vs_N > 0 | P8_T_vs_N )
P8_tmp <- P8_tmp[order(P8_tmp$position),] 
list <- data.frame(chr7_pos=c(148883700:148885800),order=1:2101)
P8_tmp <- merge(list,P8_tmp,by.x="chr7_pos",by.y="position",all.x=TRUE)
P8_tmp[is.na(P8_tmp)] <- 0
P8_Normal_EZH2 <- data.frame(P8_tmp[,c("N_ratio")])
P8_Tumor_EZH2 <- data.frame(P8_tmp[,c("T_ratio")])
names(P8_Normal_EZH2) <- "Methy"
names(P8_Tumor_EZH2) <- "Methy"
P8_Normal_EZH2$Group <- "Normal"
P8_Tumor_EZH2$Group <- "Tumor"
rownames(P8_Normal_EZH2) <- P8_tmp$chr7_pos
rownames(P8_Tumor_EZH2) <- P8_tmp$chr7_pos
P8_Normal_EZH2$SampleID <- P8_tmp$chr7_pos
P8_Tumor_EZH2$SampleID <- P8_tmp$chr7_pos
P8_Normal_EZH2 <- P8_Normal_EZH2[order(P8_Normal_EZH2$SampleID),]
P8_Tumor_EZH2 <- P8_Tumor_EZH2[order(P8_Tumor_EZH2$SampleID),]
ff_P8_Normal <- ggplot(P8_Normal_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 
ff_P8_Tumor <- ggplot(P8_Tumor_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 


P5_Normal_EZH2 <- subset(P5_Normal,V1=="chr:7")
P5_Tumor_EZH2 <- subset(P5_Tumor,V1=="chr:7")
library(rtracklayer)
library(data.table)
P5_Normal_EZH2_freq <- data.frame(position=148883700+P5_Normal_EZH2$V2,CpG_methy_ratio=P5_Normal_EZH2$V4/(P5_Normal_EZH2$V4+P5_Normal_EZH2$V5))
P5_Normal_EZH2_freq[is.na(P5_Normal_EZH2_freq)] <- 0
P5_Normal_EZH2_freq <- P5_Normal_EZH2_freq[order(P5_Normal_EZH2_freq$position),] 
P5_Tumor_EZH2_freq <- data.frame(position=148883700+P5_Tumor_EZH2$V2,CpG_methy_ratio=P5_Tumor_EZH2$V4/(P5_Tumor_EZH2$V4+P5_Tumor_EZH2$V5))
P5_Tumor_EZH2_freq[is.na(P5_Tumor_EZH2_freq)] <- 0
P5_Tumor_EZH2_freq <- P5_Tumor_EZH2_freq[order(P5_Tumor_EZH2_freq$position),] 
P5_data <- data.frame(position=P5_Normal_EZH2_freq$position,N_ratio=P5_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P5_Tumor_EZH2_freq$CpG_methy_ratio, P5_T_vs_N=P5_Tumor_EZH2_freq$CpG_methy_ratio-P5_Normal_EZH2_freq$CpG_methy_ratio)
P5_tmp <- subset(P5_data,P5_T_vs_N > 0 | P5_T_vs_N )
P5_tmp <- P5_tmp[order(P5_tmp$position),] 
list <- data.frame(chr7_pos=c(148883700:148885800),order=1:2101)
P5_tmp <- merge(list,P5_tmp,by.x="chr7_pos",by.y="position",all.x=TRUE)
P5_tmp[is.na(P5_tmp)] <- 0
P5_Normal_EZH2 <- data.frame(P5_tmp[,c("N_ratio")])
P5_Tumor_EZH2 <- data.frame(P5_tmp[,c("T_ratio")])
names(P5_Normal_EZH2) <- "Methy"
names(P5_Tumor_EZH2) <- "Methy"
P5_Normal_EZH2$Group <- "Normal"
P5_Tumor_EZH2$Group <- "Tumor"
rownames(P5_Normal_EZH2) <- P5_tmp$chr7_pos
rownames(P5_Tumor_EZH2) <- P5_tmp$chr7_pos
P5_Normal_EZH2$SampleID <- P5_tmp$chr7_pos
P5_Tumor_EZH2$SampleID <- P5_tmp$chr7_pos
P5_Normal_EZH2 <- P5_Normal_EZH2[order(P5_Normal_EZH2$SampleID),]
P5_Tumor_EZH2 <- P5_Tumor_EZH2[order(P5_Tumor_EZH2$SampleID),]
ff_P5_Normal <- ggplot(P5_Normal_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 
ff_P5_Tumor <- ggplot(P5_Tumor_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 




P4_Normal_EZH2 <- subset(P4_Normal,V1=="chr:7")
P4_Tumor_EZH2 <- subset(P4_Tumor,V1=="chr:7")
library(rtracklayer)
library(data.table)
P4_Normal_EZH2_freq <- data.frame(position=148883700+P4_Normal_EZH2$V2,CpG_methy_ratio=P4_Normal_EZH2$V4/(P4_Normal_EZH2$V4+P4_Normal_EZH2$V5))
P4_Normal_EZH2_freq[is.na(P4_Normal_EZH2_freq)] <- 0
P4_Normal_EZH2_freq <- P4_Normal_EZH2_freq[order(P4_Normal_EZH2_freq$position),] 
P4_Tumor_EZH2_freq <- data.frame(position=148883700+P4_Tumor_EZH2$V2,CpG_methy_ratio=P4_Tumor_EZH2$V4/(P4_Tumor_EZH2$V4+P4_Tumor_EZH2$V5))
P4_Tumor_EZH2_freq[is.na(P4_Tumor_EZH2_freq)] <- 0
P4_Tumor_EZH2_freq <- P4_Tumor_EZH2_freq[order(P4_Tumor_EZH2_freq$position),] 
P4_data <- data.frame(position=P4_Normal_EZH2_freq$position,N_ratio=P4_Normal_EZH2_freq$CpG_methy_ratio,T_ratio=P4_Tumor_EZH2_freq$CpG_methy_ratio, P4_T_vs_N=P4_Tumor_EZH2_freq$CpG_methy_ratio-P4_Normal_EZH2_freq$CpG_methy_ratio)
P4_tmp <- subset(P4_data,P4_T_vs_N > 0 | P4_T_vs_N )
P4_tmp <- P4_tmp[order(P4_tmp$position),] 
list <- data.frame(chr7_pos=c(148883700:148885800),order=1:2101)
P4_tmp <- merge(list,P4_tmp,by.x="chr7_pos",by.y="position",all.x=TRUE)
P4_tmp[is.na(P4_tmp)] <- 0
P4_Normal_EZH2 <- data.frame(P4_tmp[,c("N_ratio")])
P4_Tumor_EZH2 <- data.frame(P4_tmp[,c("T_ratio")])
names(P4_Normal_EZH2) <- "Methy"
names(P4_Tumor_EZH2) <- "Methy"
P4_Normal_EZH2$Group <- "Normal"
P4_Tumor_EZH2$Group <- "Tumor"
rownames(P4_Normal_EZH2) <- P4_tmp$chr7_pos
rownames(P4_Tumor_EZH2) <- P4_tmp$chr7_pos
P4_Normal_EZH2$SampleID <- P4_tmp$chr7_pos
P4_Tumor_EZH2$SampleID <- P4_tmp$chr7_pos
P4_Normal_EZH2 <- P4_Normal_EZH2[order(P4_Normal_EZH2$SampleID),]
P4_Tumor_EZH2 <- P4_Tumor_EZH2[order(P4_Tumor_EZH2$SampleID),]
ff_P4_Normal <- ggplot(P4_Normal_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 
ff_P4_Tumor <- ggplot(P4_Tumor_EZH2, aes(x = SampleID, y = Methy, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +  # 使用原始数据绘制条形图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  # 设置颜色
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转 x 轴标签
        panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = 'SampleID', y = 'Methy', title = 'Methy Levels by SampleID') +  
  scale_y_continuous(limits = c(0, 0.4)) + 
  theme(legend.title = element_blank()) 


P4 <- ff_P4_Normal/ff_P4_Tumor
P5 <- ff_P5_Normal/ff_P5_Tumor
P8 <- ff_P8_Normal/ff_P8_Tumor
P13 <- ff_P13_Normal/ff_P13_Tumor
P14 <- ff_P14_Normal/ff_P14_Tumor
P15 <- ff_P15_Normal/ff_P15_Tumor
P18 <- ff_P18_Normal/ff_P18_Tumor
ggsave(P4, file="P4_EZH2_barplot.pdf",height=8,width=6)
ggsave(P5, file="P5_EZH2_barplot.pdf",height=8,width=6)
ggsave(P8, file="P8_EZH2_barplot.pdf",height=8,width=6)
ggsave(P13, file="P13_EZH2_barplot.pdf",height=8,width=6)
ggsave(P14, file="P14_EZH2_barplot.pdf",height=8,width=6)
ggsave(P15, file="P15_EZH2_barplot.pdf",height=8,width=6)
ggsave(P18, file="P18_EZH2_barplot.pdf",height=8,width=6)
```

