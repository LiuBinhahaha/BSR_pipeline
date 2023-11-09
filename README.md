# BSR_pipeline_STAR
A pipeline of BSR-seq, without parents information
# 原理：
BSR-seq属于BSA分析方法中的一种，只不过是使用转录组测序的数据进行分析。在这里只有两个混池的转录组测序数据而缺少亲本数据，因此使用欧几里得距离（欧氏距离，Euclidean Distance，ED）进行分析，即：![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/ED.png).同时，因为没有亲本信息，在SNP位点过滤也需要注意，首先去掉低质量位点，然后去掉两池纯和且相同的位点和存在缺失值的位点。在最后画图时可以考虑使用ED4以放大信号。
虽缺少亲本测序信息，但是也可以计算SNP-index来辅助判断目标区间。根据AD计算SNP-index，公式如下：Alt_count/(Ref_count + Alt_count)，值在0-1之间。设定一个最大和最小深度用于计算频率，太大的深度可能是同源基因或者是重复序列，太低的深度在计算的时候不太准确：min.depth = 10, max.depth = 100

# 结果：
## 1. Depth_density
![image](https://github.com/LiuBinhahaha/Figs/raw/BSR_pipeline/depth_desity.png)

## 2. 查看SNP的分布
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SNP_distribution_histogram_WSL6.png)

## 3. ED plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_ED_00.png)

## 4. ED<sup>4</sup> plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_ED4_00.png)

## 5. SNP-index_WT plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_WT_index_00.png)

## 6. SNP-index_mut plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_mut_index_00.png)

## 7. ΔSNP-index plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SWL6_delta_SNP-index_00.png)
