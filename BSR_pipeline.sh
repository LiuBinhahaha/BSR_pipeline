#!/bin/bash
#STAR的基本使用流程分为两步：
#1.生成基因组索引文件。你需要提供基因组序列(FASTA)和注释文件(GTF)
#2.将读段回帖到基因组。这一步需要提供的是RNA-seq数据，格式目前都是FASTQ/FASTA, 最后会得到很多很多的文件，常规的SAM/BAM文件，可变剪切点文件，未回帖上的读段和常用于展示信号的WIG文件。

work_dir=/data/heqiang/labmember/Zhengjun/SWL6/test_BSR
sample=/data/heqiang/labmember/Zhengjun/SWL6/test_BSR/sample.txt
# 2_mut 2_WT  注意"_" !!!
genome=/data/heqiang/labmember/Zhengjun/SWL6/test_BSR/ref/Zm-B73-REFERENCE-NAM-5.0.fa
gtf=/data/heqiang/labmember/Zhengjun/SWL6/test_BSR/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.gtf
picard=/home/heqiang/biosoft/picard/picard.jar
thread=70
filename=SWL6

cd ${work_dir}

# 建立基因组索引及STAR索引
# samtools faidx ${genome} 
# gatk CreateSequenceDictionary -R ${genome} -O ${genome}.dict
# STAR --runThreadN ${thread} \
#      --runMode genomeGenerate \
#      --genomeDir ${work_dir}/ref/ \
#      --genomeFastaFiles ${genome} \
#      --sjdbGTFfile ${gtf}

# --runMode 运行模式，基因组索引构建
# --genomeDir 构建索引输出文件的目录, 一定要是存在的文件夹，需提前建好
# --genomeFastaFiles 你的基因组fasta文件所在的目录
# --sjdbGTFfile 指示的是基因组的注释文件位置，提高比对精确性，你也可以不加这一个，但是在后续比对的时候需要指定基因组注释文件，否则会报错，这里使用的是gtf文件，我尝试使用gff3文件会报错，原因是文件中缺乏strand链的+ -
# --sjdbOverhang 读段长度: 后续回帖读段的长度, 如果读长是PE 100， 则该值设为100-1=99，这个值为你测序read的长度减1，是在注释可变剪切序列的时候使用的最大长度值

#IFS_OLD=$IFS
#IFS=$'\n'       IFS_OLD=$IFS 将当前的字段分隔符值保存在变量 IFS_OLD 中，以便在以后需要时可以还原它。IFS=$'\n' 将字段分隔符 IFS 设置为换行符（\n），这意味着 Bash 将使用换行符作为字段分隔符，而不再使用默认的空格、制表符和换行符。

mkdir ${work_dir}/json
mkdir ${work_dir}/html

cd ${work_dir}
for i in $(cat ${sample});
do
        mkdir ./$i
        fastp -w 5 \
              -i ${work_dir}/${i}_R1.fq.gz -o ${work_dir}/${i}/${i}_1.clean.fastq.gz \
              -I ${work_dir}/${i}_R2.fq.gz -O ${work_dir}/${i}/${i}_2.clean.fastq.gz \
              -j ${work_dir}/json/${i}.json \
              -h ${work_dir}/html/${i}.html
        
        # # STAR第一次比对
        STAR --runThreadN ${thread} \
             --genomeDir ${work_dir}/ref \
             --readFilesIn ${work_dir}/${i}/${i}_1.clean.fastq.gz ${work_dir}/${i}/${i}_2.clean.fastq.gz \
             --readFilesCommand zcat \
             --outFileNamePrefix ${work_dir}/${i}/${i}_first
        # --genomeDir 上一步构建的索引的路径

        cd ${work_dir}/${i}
        # 再次构建STAR索引，用第一次比对的SJ.out.tab文件
        STAR --runThreadN ${thread} \
             --runMode genomeGenerate \
             --genomeDir ${work_dir}/${i}/ \
             --genomeFastaFiles ${genome} \
             --sjdbFileChrStartEnd ${work_dir}/${i}/${i}_firstSJ.out.tab \
        # --genomeDir 索引输出的路径

        # 第二次比对
        STAR --runThreadN ${thread} \
             --genomeDir ${work_dir}/${i} \
             --readFilesIn ${work_dir}/${i}/${i}_1.clean.fastq.gz ${work_dir}/${i}/${i}_2.clean.fastq.gz \
             --readFilesCommand zcat \
             --outFileNamePrefix ${work_dir}/${i}/${i}
        # --genomeDir 索引的输入为第二次构建的索引


        cd ${work_dir}
        # 添加标签、排序
        java -jar /home/heqiang/biosoft/picard/picard.jar AddOrReplaceReadGroups \
             I=${work_dir}/${i}/${i}Aligned.out.sam \
             O=${work_dir}/${i}/${i}_rg_added_sorted.bam \
             SO=coordinate \
             RGID=${i} \
             RGLB=rna \
             RGPL=ILLUMINA \
             RGPU=snpcall \
             RGSM=${i}

        # 去重复
        gatk MarkDuplicates \
             -I ${work_dir}/${i}/${i}_rg_added_sorted.bam \
             --CREATE_INDEX true \
             -O ${work_dir}/${i}/${i}_dedup.bam \
             --VALIDATION_STRINGENCY SILENT \
             -M ${work_dir}/${i}/${i}_dedup.metrics

        # MAPQ同步和reads剪切
        gatk SplitNCigarReads \
             -R ${genome} \
             -I ${work_dir}/${i}/${i}_dedup.bam \
             -O ${work_dir}/${i}/${i}_dedup_split.bam

done


file_mut=/data/heqiang/labmember/Zhengjun/SWL6/test_BSR/2_mut/2_mut_dedup_split.bam
file_WT=/data/heqiang/labmember/Zhengjun/SWL6/test_BSR/2_WT/2_WT_dedup_split.bam
# 变异检测
gatk HaplotypeCaller \
     -R ${genome} \
     -I ${file_mut} \
     -I ${file_WT} \
     --dont-use-soft-clipped-bases \
     -stand-call-conf 20.0 \
     -O ${work_dir}/${filename}.gatk.vcf
        
# 过滤低质量
gatk VariantFiltration \
     -R ${genome} \
     -V ${work_dir}/${filename}.gatk.vcf \
     -window 35 \
     -cluster 3 \
     --filter-expression "FS > 30.0" --filter-name "FS30.0" \
     --filter-expression "QD < 2.0" --filter-name "LowQD" \
     -O ${work_dir}/${filename}.flt.vcf
      
grep -vP "\tFilter" ${work_dir}/${filename}.flt.vcf > ${work_dir}/${filename}.filter.vcf

# Select variants SNPs/INDEL
gatk SelectVariants \
     -R ${genome} \
     -select-type SNP \
     -V ${work_dir}/${filename}.filter.vcf \
     -O ${work_dir}/${filename}.filter.SNPs.vcf


gatk SelectVariants \
     -R ${genome} \
     -select-type INDEL \
     -V ${work_dir}/${filename}.filter.vcf \
     -O ${work_dir}/${filename}.filter.INDELs.vcf


# 过滤
# 去除缺失位点(./.)的行; 去除野生和突变同时包含(1/1)的行，1/1.*1/1表示去除同时包含两个1/1的行，他们之间可以有任何字符，包含零个字符; 将#CHROM替换为CHROM.
grep -v "##" ${work_dir}/${filename}.filter.SNPs.vcf | grep -v "\./\." | grep -v "1/1.*1/1" | sed 's/^#CHROM/CHROM/' > ${work_dir}/${filename}.filter.SNPs.txt

