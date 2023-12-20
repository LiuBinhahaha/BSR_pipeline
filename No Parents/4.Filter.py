# -*- coding: utf-8 -*-
# 分别对BES236_RNA.txt和BES398_RNA.txt两个文件进行过滤
# 这个脚本过滤后的数据还可用于最终位点的提取
import re
file = open("/data/heqiang/labmember/Zhengjun/BSR_RNA/BES398_RNA.txt", 'r')
out = open("/data/heqiang/labmember/Zhengjun/BSR_RNA/BES398_RNA_filter.txt", 'w')

# 过滤标准
# filter QUAL<300; 
# filter P1=='0/0' & P2=='0/0'; # 这里的P1 P2指两个混池
# filter P1=='1/1' & P2=='1/1';
# filter indel
# filter missing

for line in file:
    line = line.replace("|", "/")
    if line.startswith("#CHROM"):
        _line = line.rstrip().split('\t')
        out.write('\t'.join(_line[0:]) + '\n')
    else:
        _line = line.rstrip().split('\t')
        if re.search(',', _line[4]):
            continue
        elif float(_line[5]) < 300:
            continue
        elif any(re.search(r'\.', i) for i in _line[9:]):
            continue
        else:
            tmp = _line[0:9]
            P1, P2 = _line[9:11]
            gt_P1 = P1.rstrip().split(":")[0]
            gt_P2 = P2.rstrip().split(":")[0]
            if gt_P1 == '0/0' and gt_P2 == '0/0':
                continue
            elif gt_P1 == '1/1' and gt_P2 == '1/1':
                continue
            else:
                tmp.append(P1)
                tmp.append(P2)
            out.write('\t'.join(tmp) + '\n')

file.close()
out.close()

# grep -v '^sca' BES236_RNA_filter.txt | sed s'/^#CHROM/CHROM/' > ./BES236/BES236_RNA_filter2.txt
# grep -v '^sca' BES398_RNA_filter.txt | sed s'/^#CHROM/CHROM/' > ./BES398/BES398_RNA_filter2.txt

