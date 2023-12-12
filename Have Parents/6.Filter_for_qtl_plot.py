# 按QTL-seq软件的输入格式提取，没有做任何过滤
# QTL-seq计算SNP-index，在有父母本的情况下

file = open("/data/heqiang/labmember/Zhengjun/BSR_dna_rna/BES398_filter.vcf", 'r')
out = open("/data/heqiang/labmember/Zhengjun/BSR_dna_rna/BES398_filter_qtlplot.vcf", 'w')

out.write("##fileformat=VCFv4.2" + '\n')

# 一定注意亲本的表型要与后面第一个池的表型保持一致，这与后续我们想要使用qtlplot进行snp-index分析时有关。
# 仅需提取出BES236(不抗)、BES236_mut(不抗)、BES236_WT(抗)

# 仅需提取出BES398(不抗)、BES398_mut(不抗)、BES398_WT(抗)
# BES236: 10 13 12
# BES398: 10 13 12

for line in file:
    if line.startswith("#"):
        if line.startswith("#CHROM"):
            _line = line.rstrip().split('\t')
            out.write('\t'.join(_line[0:9]) + '\t' + _line[10] + '\t' + _line[13] + '\t' + _line[12] + '\n')
    else:
        _line = line.rstrip().split("\t")
        tmp = [_line[0], _line[1], _line[2], _line[3], _line[4], _line[5], _line[6], _line[7], _line[8], _line[10], _line[13], _line[12]]
        out.write('\t'.join(tmp) + '\n')

file.close()
out.close()

# QTL-seq: qtlplot
# install: conda/micromamba install -c bioconda qtlseq
# nohup qtlplot --vcf ./BES236_filter_qtlplot.vcf --out ./BES236 --N-bulk1 50 --N-bulk2 50 -F 2 -w 1000 -s 100 --min-depth 4 --max-depth 100 --N-rep 10000 --min-SNPindex 0.3 &
# nohup qtlplot --vcf ./BES398_filter_qtlplot.vcf --out ./BES398 --N-bulk1 50 --N-bulk2 50 -F 2 -w 1000 -s 100 --min-depth 4 --max-depth 100 --N-rep 10000 --min-SNPindex 0.3 &

# --N-bulk1/2，池的大小
# -w代表滑窗，单位Kb
# -s代表步长，单位Kb
# -N-rep代表迭代次数，用于计算阈值
# 一般来说，滑窗大小选择1Mb，步长为10%的滑窗
