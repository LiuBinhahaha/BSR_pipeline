# 按照技术路线将RNA-seqvcf文件拆分
file = open("/data/heqiang/labmember/Zhengjun/BSR_RNA/RNA_seq.filter_ann.SNPs.vcf", 'r')
out = open("/data/heqiang/labmember/Zhengjun/BSR_RNA/BES236_RNA.txt", 'w')
# out = open("/data/heqiang/labmember/Zhengjun/BSR_RNA/BES398_RNA.txt", 'w')

for line in file:
    if line.startswith("#"):
        if line.startswith("#CHROM"):
            _line = line.rstrip().split("\t")
            # out.write('\t'.join(_line[0:9]) + '\t' + '\t'.join(_line[11:13]) + '\n')  # BES398
            out.write('\t'.join(_line[0:11]) + '\n')  # BES236
    else:
        _line = line.rstrip().split('\t')
        # out.write('\t'.join(_line[0:9]) + '\t' + '\t'.join(_line[11:13]) + '\n')
        out.write('\t'.join(_line[0:11]) + '\n')
file.close()
out.close()
