library(dplyr)  # deal data
library(tidyverse)
library(cowplot)
library(ggsci)
library(RColorBrewer)
#library(devtools)
#install_github('tavareshugo/windowscanr')
library(windowscanr)
# 注：脚本中的WT和mut是两个样本的名字

setwd("/data/heqiang/labmember/Zhengjun/SWL6/test")
input <- "./SWL6_snpeff.filter.SNPs.txt"
df <- read_tsv(input)
chromColor <- read_tsv("chromColor.txt")

dd <- df %>% filter(QUAL > 300) # 过滤掉QUAL小于300的位点

# 过滤去除Indel
len_ref <- nchar(dd$REF) # nchar的参数是一个字符向量，返回对应位置上字符串的长度
len_alt <- nchar(dd$ALT)
type <- len_ref == 1 & len_alt == 1
dd <- dd[type, ]

# 提取WT的基因型及其ref和alt的reads覆盖度
wt <- dd %>% select("WT") %>% separate(WT, c("WT.geno", "WT.depth"), sep=":") %>% separate(WT.depth, c("WT.ref.depth", "WT.alt.depth"), sep=",", convert = T)

# 提取mut的基因型及其ref和alt的reads覆盖度
mut <- dd %>% select("mut") %>% separate(mut, c("mut.geno", "mut.depth"), sep=":") %>% separate(mut.depth, c("mut.ref.depth", "mut.alt.depth"), sep=",", convert = T)

raw <- cbind(dd, wt, mut) %>% as.tibble()  # tidyverse中的as.tibble作用是转换数据框的格式

gt_group_type <- raw %>% group_by(WT.geno, mut.geno) %>% count()  # 查看ref和alt各种基因型组合及数目

raw <- raw %>% mutate(WT.sum = WT.ref.depth + WT.alt.depth, mut.sum = mut.ref.depth + mut.alt.depth)

# depth_density(reads) plot
pdf("depth_desity.pdf", height = 3, width = 6)
par(mfrow = c(1, 2))
plot(density(raw$WT.sum, width = 1), main = "WT", xlim = c(0, 100))
plot(density(raw$mut.sum, width = 1), main = "mut", xlim = c(0, 100))
dev.off()

# 过滤保留WT.sum和mut.sum小于10的位点
raw <- raw %>% filter(WT.sum >10, mut.sum > 10)

# 每个chrom的SNP数目
SNPnumber <- raw %>% group_by(CHROM) %>% count()


options(scipen = 200) # 在 R 中，默认情况下，当浮点数的小数部分为零时，它们将以科学计数法的形式显示。scipen 的默认值是 0，这意味着浮点数小数部分为零时将采用科学计数法。例如，0.0001 将显示为 1e-04，如果将 scipen 设置为较大的值，比如 200，那么 R 将更不容易采用科学计数法显示浮点数，即使小数部分为零，也会以普通的小数形式显示。


# SNP的染色体分布图
P1 <- chromColor %>% left_join(raw, by = "CHROM") %>% ggplot(aes(x=POS))+
  geom_histogram(aes(fill=LABEL), binwidth = 1000000)+  # binwidth 参数设置为 1000000。这意味着 x 轴上的连续范围被分成多个宽度为 1,000,000 的柱子。这将导致数据的 x 轴范围被分为多个区间，每个区间的宽度为 1,000,000，并且直方图将显示每个区间中数据点的频数。
  theme_half_open()+
  scale_x_continuous(expand = c(0, 0)) +  # expand = c(0, 0) 指定了 x 轴的显示范围不进行扩展,即设置图像是否紧贴原点
  theme(strip.text = element_text(color = NA, size = 0.1), # strip.text 和 strip.background 参数设置了分面（facet）标签的文本颜色和背景
        axis.text.y =element_text(size=5, colour = "black"), # x轴和y轴标注字体的大小
        strip.background = element_rect(color = NA, fill = NA))+
  labs(x = "POS", y = "SNP Count/1Mb")+
  facet_grid(LABEL ~ .)
P1
ggsave(P1, filename = "SNP_distribution_histogram.pdf", width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5)


# 计算ED(Euclidean distance)
raw <- raw %>% mutate(WT.ref.rate = WT.ref.depth/WT.sum,
                                              WT.alt.rate = WT.alt.depth/WT.sum,
                                              mut.ref.rate = mut.ref.depth/mut.sum,
                                              mut.alt.rate = mut.alt.depth/mut.sum,
                                              ED = sqrt((WT.ref.rate - mut.ref.rate)^2 + (WT.alt.rate - mut.alt.rate)^2),
                                              ED4 = ED^4)

# 计算SNP-index
raw <- raw %>% mutate(WT.index = WT.alt.depth/(WT.ref.depth + WT.alt.depth),
                                              mut.index = mut.alt.depth/(mut.ref.depth + mut.alt.depth),
                                              delta.index = WT.index - mut.index)

# 设定滑窗
slid_win <- winScan(x = raw,
                    groups = "CHROM",
                    position = "POS",
                    values = c("ED", "ED4", "WT.index", "mut.index","delta.index"),
                    win_size = 1000000,
                    win_step = 100000,
                    funs = "mean")

filter <- slid_win %>% select(CHROM, win_start, win_end, win_mid, ED = ED_mean, ED4 = ED4_mean, SNP_num = ED_n, WT_index = WT.index_mean, mut_index = mut.index_mean, delta_index = delta.index_mean)  # 从滑窗后的重新提取出数据
write_tsv(file = "BSR_ED_index.txt", x = slid_win)

# 添加颜色，10条染色体交替显示两种颜色
#cols =rep(c('#467aaa','#38539c'),10)
#chromColor['COLOR'] = cols[1:nrow(chromColor)]  # 将前面设定的颜色替换

raw <- chromColor %>% left_join(raw, by = "CHROM")


# 设定阈值筛选高可信窗口
# method_1 ED&ED4
# ED
EDh1 =  quantile(na.omit(filter$ED),probs=seq(0,1,0.01))['99%']
EDh2 =  quantile(na.omit(filter$ED),probs=seq(0,1,0.01))['95%']

EDh1_filter_data <- filter %>% filter(ED >= EDh1)
write_tsv(file="ED_99_filter.txt", x=EDh1_filter_data)

EDh2_filter_data <- filter %>% filter(ED >= EDh2)
write_tsv(file="ED_95_filter.txt", x=EDh2_filter_data)

# ED^4
ED4h1 =  quantile(na.omit(filter$ED4),probs=seq(0,1,0.01))['99%']
ED4h2 =  quantile(na.omit(filter$ED4),probs=seq(0,1,0.01))['95%']

ED4h1_filter_data <- filter %>% filter(ED4 >= ED4h1)
write_tsv(file="ED4_99_filter.txt", x=ED4h1_filter_data)

ED4h2_filter_data <- filter %>% filter(ED4 >= ED4h2)
write_tsv(file="ED4_95_filter.txt", x=ED4h2_filter_data)


# method_2 SNP-index
WT_index1 = quantile(na.omit(filter$WT_index), probs = seq(0,1,0.01))['99%']
WT_index2 = quantile(na.omit(filter$WT_index), probs = seq(0,1,0.01))['95%']

WT_index1_filter_data <- filter %>% filter(WT_index >= WT_index1)
write_tsv(file="WT_index_99_filter.txt", x=WT_index1_filter_data)

WT_index2_filter_data <- filter %>% filter(WT_index >= WT_index2)
write_tsv(file="WT_index_95_filter.txt", x=WT_index2_filter_data)

mut_index1 = quantile(na.omit(filter$mut_index), probs = seq(0,1,0.01))['99%']
mut_index2 = quantile(na.omit(filter$mut_index), probs = seq(0,1,0.01))['95%']

mut_index1_filter_data <- filter %>% filter(mut_index >= mut_index1)
write_tsv(file="mut_index_99_filter.txt", x=mut_index1_filter_data)

mut_index2_filter_data <- filter %>% filter(mut_index >= mut_index2)
write_tsv(file="mut_index_95_filter.txt", x=mut_index2_filter_data)

delta_index1 = quantile(na.omit(filter$delta_index), probs = seq(0,1,0.01))['99%']
delta_index2 = quantile(na.omit(filter$delta_index), probs = seq(0,1,0.01))['95%']

delta_index1_filter_data <- filter %>% filter(delta_index >= delta_index1)
write_tsv(file="delta_index_99_filter.txt", x=delta_index1_filter_data)

delta_index2_filter_data <- filter %>% filter(delta_index >= delta_index2)
write_tsv(file="delta_index_95_filter.txt", x=delta_index2_filter_data)


filter = chromColor %>% left_join(filter, by = "CHROM")

# 1.ED图(背景点用所有原始点绘制，而不是平均值)
P1 <- ggplot(na.omit(raw), aes(x = POS, y = ED)) +
  geom_point(aes(color = COLOR), size = 1.5,stroke=0,alpha=0.9) +scale_color_identity() +
  labs(x="", y="ED") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#添加window的平均值的折线图
P2 = P1+geom_line(data=filter(filter, SNP_num > 10),aes(x = win_mid, y = ED, group=LABEL),color = "black", size = 0.7)
#添加阈值线
P3 = P2+geom_hline(yintercept =EDh1, linetype = "dotdash", color ="red", size = 0.7)+geom_hline(yintercept =EDh2, linetype = "solid", color ='grey', size = 0.7)
P3
ggsave(P3, filename = "SWL6_ED.pdf", width = 15, height = 4)


#2.ED4点图同样利用所有信息绘制
P1 <- ggplot(na.omit(raw), aes(x = POS, y = ED4)) +
  geom_point(aes(color = COLOR), size = 1.5,stroke=0,alpha=0.9) +scale_color_identity() +
  labs(x="", y="ED^4") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#添加window的平均值的折线图
P2 = P1+geom_line(data=filter(filter, SNP_num > 10),aes(x = win_mid, y = ED4, group=LABEL),color = "black", size = 0.7)
#添加阈值线
P3 = P2+geom_hline(yintercept =ED4h1, linetype = "dotdash", color ="red", size = 0.7)+geom_hline(yintercept =ED4h2, linetype = "solid", color ='grey', size = 0.7)
P3
ggsave(P3, filename = "SWL6_ED4.pdf", width = 15, height = 5)


# 3. WT SNP-index plot
P1 <- ggplot(na.omit(raw), aes(x = POS, y = WT.index)) +
  geom_point(aes(color = COLOR), size = 1.5,stroke=0,alpha=0.9) +scale_color_identity() +
  labs(x="", y="WT_index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#添加window的平均值的折线图
P2 = P1+geom_line(data=filter(filter, SNP_num > 10),aes(x = win_mid, y = WT_index, group=LABEL),color = "black", size = 0.7)
#添加阈值线
P3 = P2+geom_hline(yintercept =WT_index1, linetype = "dotdash", color ="red", size = 0.7)+geom_hline(yintercept =WT_index2, linetype = "solid", color ='grey', size = 0.7)
P3
ggsave(P3, filename = "SWL6_WT_index.pdf", width = 15, height = 5)

# 4. mut SNP-index plot
P1 <- ggplot(na.omit(raw), aes(x = POS, y = mut.index)) +
  geom_point(aes(color = COLOR), size = 1.5,stroke=0,alpha=0.9) +scale_color_identity() +
  labs(x="", y="mut_index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#添加window的平均值的折线图
P2 = P1+geom_line(data=filter(filter, SNP_num > 10),aes(x = win_mid, y = mut_index, group=LABEL),color = "black", size = 0.7)
#添加阈值线
P3 = P2+geom_hline(yintercept =mut_index1, linetype = "dotdash", color ="red", size = 0.7)+geom_hline(yintercept =mut_index2, linetype = "solid", color ='grey', size = 0.7)
P3
ggsave(P3, filename = "SWL6_mut_index.pdf", width = 15, height = 4)


# 5 delta SNP-index plot
P1 <- ggplot(na.omit(raw), aes(x = POS, y = delta.index)) +
  geom_point(aes(color = COLOR), size = 1.5,stroke=0,alpha=0.9) +scale_color_identity() +
  labs(x="", y="delta_SNP-index") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
#添加window的平均值的折线图
P2 = P1+geom_line(data=filter(filter, SNP_num > 10),aes(x = win_mid, y = delta_index, group=LABEL),color = "black", size = 0.7)
#添加阈值线
P3 = P2+geom_hline(yintercept =delta_index1, linetype = "dotdash", color ="red", size = 0.7)+geom_hline(yintercept =delta_index2, linetype = "solid", color ='grey', size = 0.7)
P3
ggsave(P3, filename = "SWL6_delta_SNP-index.pdf", width = 15, height = 4)
