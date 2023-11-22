#!/bin/bash
#STAR的基本使用流程分为两步：
#1.生成基因组索引文件。你需要提供基因组序列(FASTA)和注释文件(GTF)
#2.将读段回帖到基因组。这一步需要提供的是RNA-seq数据，格式目前都是FASTQ/FASTA, 最后会得到很多很多的文件，常规的SAM/BAM文件，可变剪切点文件，未回帖上的读段和常用于展示信号的WIG文件。

workdir=/data/heqiang/labmember/Zhengjun/RNA_seq
genome=/data/heqiang/labmember/Zhengjun/RNA_seq/genome/Zm-B73-REFERENCE-NAM-5.0.fa
gtf=/data/heqiang/labmember/Zhengjun/RNA_seq/genome/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.gtf
picard=/home/heqiang/biosoft/picard/picard.jar
thread=30

cd ${workdir}

# 建立基因组索引及STAR索引,首次运行时构建索引
# samtools faidx ${genome} 
# gatk CreateSequenceDictionary -R ${genome} -O ${workdir}/genome/Zm-B73-REFERENCE-NAM-5.0.dict
# STAR --runThreadN ${thread} \
#      --runMode genomeGenerate \
#      --genomeDir ${workdir}/genome/ \
#      --genomeFastaFiles ${genome} \
#      --sjdbGTFfile ${gtf}

# --runMode 运行模式，基因组索引构建
# --genomeDir 构建索引输出文件的目录, 一定要是存在的文件夹，需提前建好
# --genomeFastaFiles 你的基因组fasta文件所在的目录
# --sjdbGTFfile 指示的是基因组的注释文件位置，提高比对精确性，你也可以不加这一个，但是在后续比对的时候需要指定基因组注释文件，否则会报错，这里使用的是gtf文件，我尝试使用gff3文件会报错，原因是文件中缺乏strand链的+ -
# --sjdbOverhang 读段长度: 后续回帖读段的长度, 如果读长是PE 100， 则该值设为100-1=99，这个值为你测序read的长度减1，是在注释可变剪切序列的时候使用的最大长度值

# mkdir ${workdir}/json
# mkdir ${workdir}/html
# mkdir ${workdir}/first_star
# mkdir ${workdir}/second_star


cd ${workdir}
sample_name=$(ls *.fq.gz | awk -F '_' '{print $1 "_" $2}' | uniq)
for i in ${sample_name};
do
    echo ${i}
    fastp -w ${thread} -i ${workdir}/${i}_R1.fq.gz -o ${workdir}/${i}_R1.clean.fq.gz -I ${workdir}/${i}_R2.fq.gz -O ${workdir}/${i}_R2.clean.fq.gz -j ${workdir}/json/${i}.json -h ${workdir}/html/${i}.html

    # STAR第一次比对
    STAR --runThreadN ${thread} \
         --genomeDir ${workdir}/genome \
         --readFilesIn ${workdir}/${i}_R1.clean.fq.gz ${workdir}/${i}_R2.clean.fq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${workdir}/first_star/${i}_first
    # --genomeDir 上一步构建的索引的路径

    # 再次构建STAR索引，用第一次比对的SJ.out.tab文件
    STAR --runThreadN ${thread} \
         --runMode genomeGenerate \
         --genomeDir ${workdir}/second_star \
         --genomeFastaFiles ${genome} \
         --sjdbFileChrStartEnd ${workdir}/first_star/${i}_firstSJ.out.tab \
    # --genomeDir 索引输出的路径

    # 第二次比对
    STAR --runThreadN ${thread} \
         --genomeDir ${workdir}/second_star \
         --readFilesIn ${workdir}/${i}_R1.clean.fq.gz ${workdir}/${i}_R2.clean.fq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${workdir}/second_star/${i}_second
    # --genomeDir 索引的输入为第二次构建的索引

    # 添加标签、排序
    java -jar /home/heqiang/biosoft/picard/picard.jar AddOrReplaceReadGroups \
         I=${workdir}/second_star/${i}_secondAligned.out.sam \
         O=${workdir}/${i}_rg_added_sorted.bam \
         SO=coordinate \
         RGID=${i} \
         RGLB=rna \
         RGPL=ILLUMINA \
         RGPU=snpcall \
         RGSM=${i}

    # 去重复
    gatk MarkDuplicates \
         -I ${workdir}/${i}_rg_added_sorted.bam \
         --CREATE_INDEX true \
         -O ${workdir}/${i}_dedup.bam \
         --VALIDATION_STRINGENCY SILENT \
         -M ${workdir}/${i}_dedup.metrics

    # MAPQ同步和reads剪切，这一步是RNA-seq特异性的一步。因为mRNA转录本是主要由DNA的外显子exon可变剪切组合而成
    gatk SplitNCigarReads \
         -R ${genome} \
         -I ${workdir}/${i}_dedup.bam \
         -O ${workdir}/${i}_dedup_split.bam
    
    echo ${i}_finished

done
# gatk SplitNCigarReads 通常在处理RNA-seq数据时更为常见，尤其是在进行基因表达分析和变异检测时。在RNA-seq数据中，由于剪切变异（splice variants）的存在，reads的CIGAR字符串可能包含 "N"，表示某些碱基被剪切掉。
# 对于DNA数据，通常在分析结构变异时也可能会使用 gatk SplitNCigarReads，因为它同样能够处理包含 "N" 的CIGAR字符串的reads。结构变异可能导致DNA片段的插入、删除或移动，因此在处理DNA数据时，尤其是在进行全基因组分析时，这个工具也可能会有一些用处。

BES236_mut=/data/heqiang/labmember/Zhengjun/RNA_seq/BES236_mut_dedup_split.bam
BES236_WT=/data/heqiang/labmember/Zhengjun/RNA_seq/BES236_WT_dedup_split.bam
BES398_mut=/data/heqiang/labmember/Zhengjun/RNA_seq/BES398_mut_dedup_split.bam
BES398_WT=/data/heqiang/labmember/Zhengjun/RNA_seq/BES398_WT_dedup_split.bam

# 变异检测
gatk --java-options -Xmx20G HaplotypeCaller \
     -R ${genome} \
     -I ${BES236_mut} \
     -I ${BES236_WT} \
     -I ${BES398_mut} \
     -I ${BES398_WT} \
     --dont-use-soft-clipped-bases \
     -stand-call-conf 20.0 \
     -O ${workdir}/RNA_seq.gatk.vcf

# 过滤低质量
gatk VariantFiltration \
     -R ${genome} \
     -V ${workdir}/RNA_seq.gatk.vcf \
     -window 35 \
     -cluster 3 \
     --filter-expression "FS > 30.0" --filter-name "FS30.0" \
     --filter-expression "QD < 2.0" --filter-name "LowQD" \
     -O ${workdir}/RNA_seq.flt.vcf

grep -vP "\tFilter" ${workdir}/RNA_seq.flt.vcf > ${workdir}/RNA_seq.filter.vcf

# Select variants SNPs/INDEL
gatk SelectVariants \
     -R ${genome} \
     -select-type SNP \
     -V ${workdir}/RNA_seq.filter.vcf \
     -O ${workdir}/RNA_seq.filter.SNPs.vcf


gatk SelectVariants \
     -R ${genome} \
     -select-type INDEL \
     -V ${workdir}/RNA_seq.filter.vcf \
     -O ${workdir}/RNA_seq.filter.INDELs.vcf

