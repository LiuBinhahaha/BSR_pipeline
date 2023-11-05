# BSR_pipeline_STAR
A pipeline of BSR-seq, without parents information
# 原理：
BSR-seq属于BSA分析方法中的一种，只不过是使用转录组测序的数据进行分析。在这里只有两个混池的转录组测序数据而缺少亲本数据，因此使用欧几里得距离（欧氏距离，Euclidean Distance，ED）进行分析，即：
