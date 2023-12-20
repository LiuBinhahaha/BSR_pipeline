file_interval = open("/data/heqiang/labmember/Zhengjun/SWL7/ED_99_filter.txt", 'r')
vcf_file = open("/data/heqiang/labmember/Zhengjun/SWL7/SWL7_snpeff.filter.SNPs.txt", 'r')  # 注释后的vcf文件
out = open("/data/heqiang/labmember/Zhengjun/SWL7/SWL7_ED_99_SNP.txt", 'w')

hd = []
dic = {}  # 存储vcf文件的信息             

for lines in vcf_file:
    if lines.startswith("CHROM"):
        _lines = lines.rstrip().split('\t')
        hd += _lines

hd += 'gt_WT', 'gt_mut'  # 添加基因型，AA、AT、TT...
vcf_file.seek(0,0)
out.write('\t'.join(hd) +'\n')

for i in range(1,11):
    inner_dic = {}
    pattern = f'chr{i}'

    vcf_file.seek(0,0)
    for linevcf in vcf_file:
        if linevcf.startswith("CHROM"):
            continue
        else:
            _linevcf = linevcf.rstrip().split('\t')
            if _linevcf[0] == pattern:
                if _linevcf[1] not in inner_dic:
                    inner_dic[_linevcf[1]] = []
                inner_dic[_linevcf[1]] += _linevcf
    if inner_dic:
        dic[pattern] = inner_dic
#print(dic['chr1']["9574645"])  # 测试

for line in file_interval:
    _line = line.rstrip().split('\t')
    if line.startswith("CHROM"):
        continue
    else:
        _line = line.rstrip().split('\t')
        chr = _line[0]
        start = int(_line[3])  # win_start
        end = int(_line[4])  # win_end
        if chr in dic:
            inner_dic = dic[chr]
            for pos, data in inner_dic.items():
                if start <= int(pos) <= end:
                    data = dic[chr][pos]
                    gt1 = data[9].rstrip().split(':')[0]
                    if gt1 == '0/0':
                        _gt1 = f'{data[3]}/{data[3]}'
                    elif gt1 == '0/1':
                        _gt1 = f'{data[3]}/{data[4]}'
                    elif gt1 == '1/1':
                        _gt1 = f'{data[4]}/{data[4]}'
                    gt2 = data[10].rstrip().split(':')[0]
                    _gt2 = '-'
                    if gt2 == '0/0':
                        _gt2 = f'{data[3]}/{data[3]}'
                    elif gt2 == '0/1':
                        _gt2 = f'{data[3]}/{data[4]}'
                    elif gt2 == '1/1':
                        _gt2 = f'{data[4]}/{data[4]}'
                    tmp = [str(e) for e in data]
                    tmp.append(_gt1)
                    tmp.append(_gt2)
                    out.write('\t'.join(tmp) + '\n')
