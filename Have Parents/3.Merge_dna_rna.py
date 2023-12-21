# 合并DNA-seq和RNA的vcf文件，取交集
import re
import os
# ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'B73', 'BES236', 'BES398', 'Zh58']

file_dna = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/DNA_seq_ann.filter.snp.vcf", 'r')
file_rna = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/RNA_seq.filter_ann.SNPs.vcf", 'r')
out = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/merged_dna_rna.vcf", 'w')

dna = {}

hd = []

for line_d in file_dna:
    if line_d.startswith("#"):
           if line_d.startswith("#CHROM"):
               _line_d = line_d.rstrip().split('\t')
               hd += _line_d
    else:
         _line_d = line_d.rstrip().split('\t')
         dna[(_line_d[0], _line_d[1], _line_d[2], _line_d[3], _line_d[4])] = _line_d[5:]

for line in file_rna:
     if line.startswith("#"):
          if line.startswith("#CHROM"):
               _line = line.rstrip().split('\t')
               hd += _line[9:]
               break
file_rna.seek(0,0)

out.write("##fileformat=VCFv4.2" + '\n')
out.write('\t'.join(hd[0:]) + '\n')

for line in file_rna:
     if line.startswith("#"):
          continue
     else:
          _line = line.rstrip().split('\t')
          tmp = (_line[0], _line[1], _line[2], _line[3], _line[4])  # rna
          if tmp in dna:
               dna_seq = dna[tmp]
               out.write('\t'.join(tmp) + '\t' + '\t'.join(dna_seq) + '\t' + '\t'.join(_line[9:]) + '\n')

file_dna.close()
file_rna.close()
out.close()

