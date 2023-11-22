# BSR_pipeline_STAR
A pipeline of BSR-seq, with parents information
# 实验设计：
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/trail_design.png)
# 原理：
BSR-seq是将BSA与RNA-seq结合起来的分析方法，与重测序BSA不同的是，在分离群体中选择极端性状的个体构建两个池，提取两个池的总RNA，进行转录组测序，根据混池个数和物种的基因组大小确定测序的数据量。另外，BSR还能反映混池的基因表达量信息，可以根据候选区段内相关基因的表达差异筛选候选基因。基于欧几里得距离(Euclidena distance, ED)，也称为MMAPPR，或SNP-index找到与突变基因共分离的SNP位点和区间，最后预测候选基因。ED主要用于无亲本测序信息的BSA分析，通过计算不同混池间各突变型的频率距离，采用距离差异来反映标记与目标区域的连锁强度，画图时可以考虑使用ED^4以放大信号。SNP-index来表示子代群体与亲本之间的序列差异程度，指的是在特定位点上，携带有不同于参考亲本的SNP的reads数占比对到同一位点的所有reads总数的比值。将两个极端池的SNP-index相减获得ΔSNP-index，Marker与性状关联度越强，ΔSNP-index越接近于1。
ED公式：![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/ED1.png)
SNP-index公式：SNP-index = Alt_count/(Ref_count + Alt_count)

过滤标准：首先去掉低质量位点(QUAL<300)，然后去掉两池纯和且相同的位点和存在缺失值的位点。

# 结果：
## 1. Depth_density
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/depth_desity.png)

## 2. 查看SNP的分布
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SNP_distribution_histogram.png)

## 3. ED plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_ED.png)

## 4. ED<sup>4</sup> plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_ED4.png)

## 5. SNP-index_WT plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_WT_index.png)

## 6. SNP-index_mut plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_mut_index.png)

## 7. ΔSNP-index plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_delta_index.png)
