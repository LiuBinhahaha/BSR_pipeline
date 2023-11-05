# BSR_pipeline_STAR
A pipeline of BSR-seq, without parents information
# 原理：
BSR-seq属于BSA分析方法中的一种，只不过是使用转录组测序的数据进行分析。在这里只有两个混池的转录组测序数据而缺少亲本数据，因此使用欧几里得距离（欧氏距离，Euclidean Distance，ED）进行分析，即：![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/ED.png).同时，因为没有亲本信息，在SNP位点过滤也需要注意，首先去掉低质量位点，然后去掉两池纯和且相同的位点和存在缺失值的位点。在最后画图时可以考虑使用ED4以放大信号。
# 结果：
## 1. Depth_density
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/depth_desity.png)

## 2. 查看SNP的分布
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/SNP_distribution_histogram_WSL6.png)

## 3. ED plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/BSR_ED.png)

## 4. ED<sup>4</sup> plot
![image](https://github.com/LiuBinhahaha/Figs/blob/main/BSR_pipeline/BSR_ED4.png)
