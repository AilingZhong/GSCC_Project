# WES analysis of GSCC tumors

```shell

vim GSCC_Tumor.csv
GSCC_1 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/3_1857655/3_1857655
GSCC_2 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/4_1830584/4_1830584
GSCC_3 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/7_1800393/7_1800393
GSCC_4 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/8_1808168/8_1808168
GSCC_5 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/9_2117121/9_2117121
GSCC_6 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/12_1728429/12_1728429
GSCC_7 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/14_2048944/14_2048944
GSCC_8 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/15_2151232/15_2151232
GSCC_9 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/16_1747087/16_1747087
GSCC_10 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/17_1818629/17_1818629
GSCC_11 /mnt/data3/sequencedata/WES/WES_12_ZMS_GSCC_20230727_11samples/CP2023051200040/H101SC23070355/RSHD00204/X101SC23070355-Z01/X101SC23070355-Z01-J001/00.CleanData/18_1731790/18_1731790
GSCC_12 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P12_Tumor_WES_H2O/P12_Tumor_WES_H2O_FIPM240041745-1A_E200032701_L1
GSCC_13 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P13_Tumor_WES_H2O/P13_Tumor_WES_H2O_FIPM240041747-1A_E200032701_L1
GSCC_14 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/Merge_P14_Tumor_WES/Merge_P14_Tumor_WES
GSCC_15 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P15_Tumor_WES_PET/P15_Tumor_WES_PET_FIPM240041749-1A_E200032701_L1
GSCC_16 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P16_Tumor_WES_PET/P16_Tumor_WES_PET_FIPM240041750-1A_E200032701_L1
GSCC_17 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P17_Tumor_WES_H2O/P17_Tumor_WES_H2O_FIPM240041748-1A_E200032701_L1
GSCC_18 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P18_Tumor_WES_PET/P18_Tumor_WES_PET_FIPM240041743-1A_E200032701_L1
GSCC_19 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P19_Tumor_WES_PET/P19_Tumor_WES_PET_FIPM240041746-1A_E200032701_L1
GSCC_20 /mnt/data3/sequencedata/WES/WES_22_ZhangMengsha_GSCC_human_9samples/CP2024093000033/H101SC24103959/RSHD00204/X101SC24103959-Z01/X101SC24103959-Z01-F021/01.RawData/P20_Tumor_WES_PET/P20_Tumor_WES_PET_FIPM240041751-1A_E200032701_L1


GENOME_index=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_bwa_index/hg19_bwa_index
GATK=/mnt/data2/userdata/shuhao/programme/GATK/gatk-4.4.0.0/gatk
cat GSCC_Tumor.csv | while read id ; do
arr=($id)
fq2=${arr[1]}'_1.fq.gz'
fq1=${arr[1]}'_2.fq.gz'
sample=${arr[0]}
echo ${fq2}
echo ${fq1}
echo $sample
bwa mem -t 30 -M -R "@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:Illumina" $GENOME_index $fq1 $fq2 > $sample.sam
$GATK  --java-options "-Xmx40G -Djava.io.tmpdir=./"  SortSam -SO coordinate  -I $sample.sam -O $sample.sort.bam > $sample.log.sort;
done



GENOME_index=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_bwa_index/hg19_bwa_index
GATK=/mnt/data2/userdata/shuhao/programme/GATK/gatk-4.4.0.0/gatk
cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK  --java-options "-Xmx50G -Djava.io.tmpdir=./"  MarkDuplicates \
-I $sample.sort.bam \
-O ${sample}_marked.bam \
--REMOVE_DUPLICATES true \
-M $sample.metrics 1>$sample.log.mark  2>&1;
done

GENOME_index=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_bwa_index/hg19_bwa_index
GATK=/mnt/data2/userdata/shuhao/programme/GATK/gatk-4.4.0.0/gatk
cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx60G -Djava.io.tmpdir=./" FixMateInformation \
-I ${sample}_marked.bam \
-O ${sample}_fix.bam \
-SO coordinate --CREATE_INDEX true;
done


GENOME=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_genome.fa
GENOME_index=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_bwa_index/hg19_bwa_index
GATK=/mnt/data2/userdata/shuhao/programme/GATK/gatk-4.4.0.0/gatk
hg19_vcf_1=/mnt/data2/userdata/abao/1_project/reference/with_chr_hg19_v0_Homo_sapiens_assembly19.dbsnp.vcf
hg19_vcf_2=/mnt/data2/userdata/abao/1_project/reference/with_chr_hg19_v0_Homo_sapiens_assembly19.dbsnp138.vcf

cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx40G -Djava.io.tmpdir=./" BaseRecalibrator \
-R $GENOME -I ${sample}_fix.bam \
--known-sites $hg19_vcf_1 \
--known-sites $hg19_vcf_2 \
-O $sample.wes.recal_data.table;
done



cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx50G -Djava.io.tmpdir=./" ApplyBQSR \
-R $GENOME -I ${sample}_fix.bam \
-bqsr $sample.wes.recal_data.table \
-O ${sample}_BQSR.bam ;
done



GENOME=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_genome.fa
GENOME_index=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_bwa_index/hg19_bwa_index
GATK=/mnt/data2/userdata/shuhao/programme/GATK/gatk-4.4.0.0/gatk
hg19_vcf_1=/mnt/data2/userdata/abao/1_project/reference/with_chr_hg19_v0_Homo_sapiens_assembly19.dbsnp.vcf
hg19_vcf_2=/mnt/data2/userdata/abao/1_project/reference/with_chr_hg19_v0_Homo_sapiens_assembly19.dbsnp138.vcf


********************************Using the normal PON_hg19 provided by the GATK database as the control*************************************************************************
********************************Using the normal PON_hg19 provided by the GATK database as the control*************************************************************************

PON=/mnt/data2/userdata/abao/8922/abao/project/ZMS_EZH2_WES/chr_hg19_PON.vcf.gz
Af_files=/mnt/data2/userdata/abao/1_project/reference/af-only-gnomad.raw.sites.hg19.vcf.gz
EXON_BED=/mnt/data2/userdata/abao/1_project/reference/hg19_Exome-Agilent_V6.bed
GENOME=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_genome.fa
GENOME_index=/mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_bwa_index/hg19_bwa_index
GATK=/mnt/data2/userdata/shuhao/programme/GATK/gatk-4.4.0.0/gatk
hg19_vcf_1=/mnt/data2/userdata/abao/1_project/reference/with_chr_hg19_v0_Homo_sapiens_assembly19.dbsnp.vcf
hg19_vcf_2=/mnt/data2/userdata/abao/1_project/reference/with_chr_hg19_v0_Homo_sapiens_assembly19.dbsnp138.vcf

cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx50G -Djava.io.tmpdir=./" Mutect2 \
-R $GENOME \
-L $EXON_BED \
-I ${sample}_BQSR.bam \
-germline-resource $Af_files \
--panel-of-normals $PON \
--f1r2-tar-gz ${sample}_f1r2.tar.gz \
-O ${sample}_somatic.vcf.gz ;
done


#Run GetPileupSummaries to summarize read support for a set number of known variant sites.
cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx60G -Djava.io.tmpdir=./" GetPileupSummaries \
-I ${sample}_BQSR.bam \
-V $Af_files \
-L $EXON_BED \
-O ${sample}_getpileupsummaries.table ;
done


cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx60G -Djava.io.tmpdir=./" CalculateContamination \
-I ${sample}_getpileupsummaries.table \
-tumor-segmentation ${sample}_segments.table \
-O ${sample}_calculatecontamination.table ;
done


cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx50G -Djava.io.tmpdir=./" LearnReadOrientationModel \
-I ${sample}_f1r2.tar.gz \
-O ${sample}_read_orientation_model.tar.gz ;
done


cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
$GATK --java-options "-Xmx40G -Djava.io.tmpdir=./" FilterMutectCalls \
-V /mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/batch2_unfilter_vcf/${sample}_filtered_with_PON.vcf \
-R $GENOME \
-L $EXON_BED \
--tumor-segmentation /mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/contamination_segments_getpileupsummaries/${sample}_segments.table \
--contamination-table /mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/contamination_segments_getpileupsummaries/${sample}_calculatecontamination.table \
--ob-priors /mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/LearnReadOrientationModel/${sample}_read_orientation_model.tar.gz \
-O /mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/${sample}_FilterMutectCalls.vcf;
done




**********************************VCF_to_MAF******************************************************
**********************************VCF_to_MAF******************************************************
**********************************VCF_to_MAF******************************************************
**********************************VCF_to_MAF******************************************************

cat GSCC_Tumor.csv | while read id ; do
arr=($id)
sample=${arr[0]}
echo $sample
perl /mnt/data2/userdata/xiangyu/programme/vcf2maf-master/vcf2maf.pl \
--input-vcf ${sample}_FilterMutectCalls.vcf \
--output-maf $sample.Final.vep.maf \
--tumor-id $sample \
--vcf-tumor-id $sample \
--ref-fasta /mnt/data2/userdata/abao/1_project/reference/bwa_reference/hg19_genome.fa \
--vep-path /mnt/data2/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/ \
--filter-vcf /mnt/data2/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
--vep-data /mnt/data2/userdata/xiangyu/programme/ensembl-vep-release-95.3/vep_all/.vep \
--vep-forks 40 ;
done

****************************CNV analysis******************************************************
****************************CNV analysis******************************************************
****************************CNV analysis******************************************************
cnvkit.py batch \
/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/All_bam_bai/*_BQSR.bam \
--normal /mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/Normal_as_control/*Normal_BQSR.bam \
--target-avg-size 500 \
--targets /mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/reference/Exome-Agilent_V6.bed \
--fasta /mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/reference/hg19.fa \
--annotate /mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/reference/refFlat.txt \
--access /mnt/data/user_data/xiangyu/programme/cnvkit-0.9.9/data/access-5k-mappable.hg19.bed \
--output-reference /mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/500bp_All_GSCC_add_PBMC_as_contrl_my_reference.cnn \
--output-dir /mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/All_GSCC_add_PBMC_as_contrl_CNVkit_500bp/ \
--diagram --scatter -p 30


```

```R

library(maftools)
library(ggplot2)

files <- c(
  "GSCC_20" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_20.Final.vep.maf",
  "GSCC_19" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_19.Final.vep.maf",
  "GSCC_18" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_18.Final.vep.maf",
  "GSCC_17" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_17.Final.vep.maf",
  "GSCC_16" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_16.Final.vep.maf",
  "GSCC_15" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_15.Final.vep.maf",
  "GSCC_14" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_14.Final.vep.maf",
  "GSCC_13" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_13.Final.vep.maf",
  "GSCC_12" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_12.Final.vep.maf",
  "GSCC_11" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_11.Final.vep.maf",
  "GSCC_10" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_10.Final.vep.maf",
  "GSCC_9" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_9.Final.vep.maf",
  "GSCC_8" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_8.Final.vep.maf",
  "GSCC_7" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_7.Final.vep.maf",
  "GSCC_6" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_6.Final.vep.maf",
  "GSCC_5" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_5.Final.vep.maf",
  "GSCC_4" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_4.Final.vep.maf",
  "GSCC_3" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_3.Final.vep.maf",
  "GSCC_2" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_2.Final.vep.maf",
  "GSCC_1" = "/mnt/data2/userdata/abao/ZhangMengsha/Combine_batch1_and_batch2_WES_GSCC/final_Filter_VCF/GSCC_1.Final.vep.maf"
)

all_samples <- list()

for (sample_name in names(files)) {
  maf_data <- read.maf(files[sample_name])
  maf_df <- data.frame(maf_data@data)
  maf_df$sample <- sample_name
  all_samples[[sample_name]] <- maf_df
}
All <- do.call(rbind, all_samples)
All$Frequency <- as.numeric(All$t_alt_count) / as.numeric(All$t_depth)
tmp <- subset(All, Frequency < 1 & t_depth >= 10 & t_alt_count >= 3)
PASS <- subset(tmp, FILTER == "PASS")
write.csv(PASS, file = "Final_filter_mutation_calling.csv", row.names = FALSE)


PASS_maf <- read.maf(PASS)
PASS_maf <- read.maf(PASS_maf)
pdf("EZH2_lollipopplot.pdf",width=8,height=3)
lollipopPlot( maf = PASS_maf,gene = 'EZH2',showMutationRate = TRUE)
dev.off()


&&&&&&&&&&&&&&&&&&&&&&&&&&&Integration_of_CNV_and_SNV&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&Integration_of_CNV_and_SNV&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&Integration_of_CNV_and_SNV&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&Integration_of_CNV_and_SNV&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

cns_dir <- "/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/All_GSCC_add_PBMC_as_contrl_CNVkit_500bp"
cns_files <- list.files(cns_dir, pattern = "*.call.cns$", full.names = TRUE)
combined_cns <- data.frame()

for (file in cns_files) {
  cns <- read.table(file, sep = "\t", header = TRUE)
  sample_name <- tools::file_path_sans_ext(basename(file))
  cns$sample <- sample_name
  combined_cns <- rbind(combined_cns, cns)
}

cnv_data <- subset(combined_cns, cn == "1" | cn == "0" | cn > 2 )
cnv_data_long <- cnv_data %>%
  separate_rows(gene, sep = ",") %>%
  select(gene, sample, cn)  # 添加 cn 列
cnv_data <- data.frame(cnv_data_long)

library(dplyr)
library(readr)

cnv_data <- cnv_data %>%
  mutate(Type = case_when(
    cn > 3  ~ "Amp",
    cn == 3 ~ "Gain",
    cn == 0 ~ "Deep_Del",
    cn == 1 ~ "Shallow_Del",
    cn == 2 ~ "Neutral" 
  ))

cnv_data <- cnv_data %>% mutate(sample = gsub("_BQSR.call", "", sample))
cnv_data$CN <- cnv_data$cn
CN_TABLE <- cnv_data[,c("gene","sample","Type")]
colnames(CN_TABLE) <- c("Gene", "Sample_name", "CN")

library(dplyr)
library(readr)


maf_data <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/Final_filter_mutation_calling.csv")
oncoKB_from_MSK_database <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/oncoKB_from_MSK_database.csv")
oncoKB_from_MSK_database <- subset(oncoKB_from_MSK_database, OncoKB.Annotated=="Yes")
GSCC_onco <- subset(maf_data, maf_data$Hugo_Symbol %in% oncoKB_from_MSK_database$Hugo.Symbol)
GSCC_onco <- GSCC_onco[order(-GSCC_onco$Frequency),]
GSCC_maf = read.maf(maf=GSCC_onco,cnTable = CN_TABLE)



maf_data <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/Final_filter_mutation_calling.csv")
oncoKB_from_MSK_database <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/1_Rebuttal_GSCC_for_NatureCommunications/WES_GSCC/oncoKB_from_MSK_database.csv")
oncoKB_from_MSK_database <- subset(oncoKB_from_MSK_database, OncoKB.Annotated=="Yes")
GSCC_onco <- subset(maf_data, maf_data$Hugo_Symbol %in% oncoKB_from_MSK_database$Hugo.Symbol)
GSCC_onco <- GSCC_onco[order(-GSCC_onco$Frequency),]
GSCC_maf = read.maf(maf=GSCC_onco,cnTable = CN_TABLE)

library(ggplot2)
pdf("GSCC_mut_plot_final.pdf",width=15,height=18)
oncoplot(maf = GSCC_maf, genes=c("KMT2A", "KMT2D", "KMT2C", "KMT2B", "TP53", "CREBBP", "EP300", "EZH2","PIK3CA", "MTOR", "PIK3R1", "TSC1", "TSC2", "PTEN", "AKT1", "AKT2"),clinicalFeatures ="Tumor_Sample_Barcode",keepGeneOrder = FALSE)
dev.off()
```

