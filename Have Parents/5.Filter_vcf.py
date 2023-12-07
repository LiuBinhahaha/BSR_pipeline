# 分别对BES236.vcf和BES398.vcf两个文件进行过滤
# 这个脚本过滤后的数据还可用于最终位点的提取
import re
file = open("/data/heqiang/labmember/Zhengjun/BSR_dna_rna/BES398.vcf", 'r')
out = open("/data/heqiang/labmember/Zhengjun/BSR_dna_rna/BES398_filter.vcf", 'w')

out.write("##fileformat=VCFv4.2" + '\n')  # vcf文件的头格式
# 过滤标准
# filter QUAL<300; 
# filter P1=='0/0' & P2=='0/0'; 
# filter compose1=='0/0' & compose2=='0/0';
# filter compose1=='1/1' & compose2=='1/1'.
# filter indel
# filter missing

for line in file:
    line = line.replace("|", "/")
    if line.startswith("#"):
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
            tmp = _line[0:10]
            P1, P2, compose1, compose2 = _line[10:14]
            gt_P1 = P1.rstrip().split(":")[0]
            gt_P2 = P2.rstrip().split(":")[0]
            gt_compose1 = compose1.rstrip().split(":")[0]
            gt_compose2  =compose2.rstrip().split(":")[0]
            if gt_P1 == '0/0' and gt_P2 == '0/0':
                continue
            else:
                tmp.append(P1)
                tmp.append(P2)
            if gt_compose1 == '0/0' and gt_compose2 == '0/0':
                continue
            elif gt_compose1 == '1/1' and gt_compose2 == '1/1':
                continue
            else:
                tmp.append(compose1)
                tmp.append(compose2)
            out.write('\t'.join(tmp) + '\n')

file.close()
out.close()
