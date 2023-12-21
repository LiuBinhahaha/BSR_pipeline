# 将上一步取交集后的vcf文件分成BES236和BES398两个vcf文件，用于后续的分析
import re
file = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/merged_dna_rna.vcf", 'r')
out = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/BES398.vcf", 'w')

out.write("##fileformat=VCFv4.2" + '\n')
# BES236: 9 10 12 13 14
# BES398: 9 11 12 15 16

for line in file:
    if line.startswith("#"):
        if line.startswith("#CHROM"):
            _line = line.rstrip().split('\t')
            out.write('\t'.join(_line[0:10]) + '\t' + _line[11] + '\t' + _line[12] + '\t' + _line[15] + '\t' + _line[16] + '\n')
    else:
        _line = line.rstrip().split('\t')
        tmp = [_line[0], _line[1], _line[2], _line[3], _line[4], _line[5], _line[6], _line[7], _line[8], _line[9], _line[11], _line[12], _line[15], _line[16]]
        out.write('\t'.join(tmp) + '\n')
out.close()
file.close()
