#过滤
#去除缺失位点(./.)的行; 去除野生和突变同时包含(1/1)的行，1/1.*1/1表示去除同时包含两个1/1的行，他们之间可以有任何字符，包含零个字符; 将#CHROM替换为CHROM.
grep -v '##' SWL6_snpeff.filter.SNPs.vcf | grep -v '\./\.' | grep -v '0/0.*0/0' | grep -v '1/1.*1/1'| sed s'/^#CHROM/CHROM/' | sed s'/2_WT/WT/' | sed s'/2_mut/mut/' > SWL6_snpeff.filter.SNPs.txt

# -v表示反向过滤，即不包含(delete)目标所在的行
# \./\. : ./.
# '0/0.*0/0' : 因为这个项目的BSR-seq数据仅有两个池，没有亲本，所以可以通过判断两个池中的是否含有 要过滤 的基因型的SNP
