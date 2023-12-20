import re
import os

# 上一步获得SNP-index的结果中不包含滑窗的win_start和win_end，需要自己根据窗口的大小去补充(-499999, +500000)
# 由于sliding_window含有overlap，所以使用bedtools去除overlap后再提取位点，以下是具体的去除过程：
# awk -F '\t' -v OFS='\t' '{print $1, $2, $3}' sliding_window.p95.tsv > tmp_sliding_window.p95.bed
# bedtools merge -i tmp_sliding_window.p95.bed > no_overlap_sliding_window.p95.bed

file_interval = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/BES236/no_overlap_sliding_window.p95.bed", 'r')
vcf_file = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/BES236_filter.vcf", 'r')  # 过滤后的vcf文件
out = open("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/BES236/sliding_window.p95_extract.csv", 'w')

hd = []
dic = {}  # 存储vcf文件的信息             

for lines in vcf_file:
    if lines.startswith("#"):
        if lines.startswith("#CHROM"):
            _lines = lines.rstrip().split("\t")
            hd += _lines

vcf_file.seek(0,0)
out.write(','.join(hd) +'\n')

# Chr
# 创建vcf_file字典
for i in range(1,11):
    inner_dic = {}
    pattern = f'chr{i}'

    vcf_file.seek(0,0)
    for linevcf in vcf_file:
        if linevcf.startswith("#"):
            if lines.startswith("#CHROM"):
                continue
        else:
            _linevcf = linevcf.rstrip().split('\t')
            if _linevcf[0] == pattern:
                if _linevcf[1] not in inner_dic:
                    inner_dic[_linevcf[1]] = []
                inner_dic[_linevcf[1]] += _linevcf
    if inner_dic:
        dic[pattern] = inner_dic
# print(dic['chr1']["38511"])  # 测试

for line in file_interval:
    _line = line.rstrip().split('\t')
    if line.startswith("CHROM"):  # bed文件无表头，CHROM忽略
        continue
    else:
        _line = line.rstrip().split('\t')
        chr = _line[0]
        start = int(_line[1])  # win_start
        end = int(_line[2])  # win_end
        if chr in dic:
            inner_dic = dic[chr]
            for pos, data in inner_dic.items():
                if start <= int(pos) <= end:
                    data = dic[chr][pos]
                    tmp = data[0:9]
                    for s in data[9:]:
                        gt = '-'
                        _s = s.rstrip().split(':')[0]
                        # print(_s)
                        if _s == '0/0':
                            gt = f'{data[3]}/{data[3]}'
                        elif _s == '0/1':
                            gt = f'{data[3]}/{data[4]}'
                        elif _s == '1/1':
                            gt = f'{data[4]}/{data[4]}'
                        tmp.append(gt)
                    out.write(','.join(tmp) + '\n')

file_interval.close()
vcf_file.close()
out.close()

