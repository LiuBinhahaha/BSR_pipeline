# build构建参考基因组注释文件：
java -jar snpEff.jar build -gff3 -v setaria

# 开始注释
java -jar snpEff.jar Setaria_italica sample.vcf > sample.ann.vcf
