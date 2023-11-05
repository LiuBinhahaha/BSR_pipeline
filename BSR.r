library(tidyverse)
library(cowplot)
library(ggsci)
library(RColorBrewer)
#library(devtools)
#install_github('tavareshugo/windowscanr')
library(windowscanr)

setwd("/data/heqiang/labmember/Zhengjun/SWL7")

input <- "./SWL7.filter.SNPs.txt"
(df <- read_tsv(input))
(chromColor <- read_tsv("./chromColor.txt"))
dd <- df %>% filter(QUAL > 300)  # 过滤QUAL低于300的位点

l1<-nchar(dd$REF)	#求REF字符串长度
l2<-nchar(dd$ALT)	#求ALT字符串长度

#过滤indel行，保留SNP信息
lse <- l1 == 1 & l2 == 1
dd <- dd[lse, ]

dd1 <- dd %>% separate(WT, c("WT.geno", "WT"), sep = ":") %>%
  separate(mut, c("mut.geno", "mut"), sep = ":")
geno <- dd1 %>% select(WT.geno, mut.geno)  # geno为WT和mut每个位点的基因型

dd2 <- dd1 %>% select(WT, mut)  #筛选出来每个reads的覆盖度
dd3 <- dd2 %>% 
  separate(WT, c("WT.ref.depth", "WT.alt.depth"), sep = ",", convert = T) %>%
  separate(mut, c("mut.ref.depth", "mut.alt.depth"), sep = ",", convert = T) 

hh <- cbind(dd, geno, dd3) %>% as_tibble()
hh %>% group_by(WT.geno, mut.geno) %>% count()

hh <- hh %>% mutate(WT.sum = WT.ref.depth + WT.alt.depth,
                    mut.sum = mut.ref.depth + mut.alt.depth)  #加和

pdf("depth_desity.pdf", height = 3, width = 6)
par(mfrow = c(1, 2))
plot(density(hh$WT.sum, width = 1), main = "WT", xlim = c(0, 100))
plot(density(hh$mut.sum, width = 1), main = "mut", xlim = c(0, 100))
dev.off()
png("depth_desity.png", height = 3, width = 6, units = "in", res = 500)
par(mfrow = c(1, 2))
plot(density(hh$WT.sum, width = 1), main = "WT", xlim = c(0, 100))
plot(density(hh$mut.sum, width = 1), main = "mut", xlim = c(0, 100))
dev.off()
print(hh$WT.sum)
hh <- hh %>% filter(WT.sum > 10, mut.sum > 10)  # 过滤保留总depth大于10的位点

SNPnumber <- hh %>% group_by(CHROM) %>% count()
write_tsv(x = SNPnumber, path = "SNP_number_per_chr.txt")
write_csv(x = SNPnumber, path = "SNP_number_per_chr.csv")


options(scipen = 200)  # 在 R 中，默认情况下，当浮点数的小数部分为零时，它们将以科学计数法的形式显示。scipen 的默认值是 0，这意味着浮点数小数部分为零时将采用科学计数法。例如，0.0001 将显示为 1e-04，如果将 scipen 设置为较大的值，比如 200，那么 R 将更不容易采用科学计数法显示浮点数，即使小数部分为零，也会以普通的小数形式显示。
colourCount = dim(chromColor)[[1]]
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # brewer.pal(8, "Set1") 是使用 R 包中的 RColorBrewer 中的函数 brewer.pal，它用于生成颜色调色板。在这里，它要求生成 8 种颜色，从 "Set1" 调色板中选择颜色。
Phist <- chromColor %>% left_join(hh, by = "CHROM") %>% ggplot(aes(x = POS)) +  # 直方图。chromColor 和 hh 是两个数据框，通过 left_join 函数使用 "CHROM" 列将它们合并
  geom_histogram(aes(fill = LABEL), color = NA, binwidth = 1000000) +
  labs(x = NULL, y = "SNP Count / 1Mb") +
  scale_x_continuous(expand = c(0, 0)) +  # expand = c(0, 0) 指定了 x 轴的显示范围不进行扩展
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_half_open() +
  theme(strip.text = element_text(color = NA, size = 0.1), # strip.text 和 strip.background 参数设置了分面（facet）标签的文本颜色和背景
        axis.text.y =element_text(size=5, colour = "black"), # x轴和y轴标注字体的大小
        strip.background = element_rect(color = NA, fill = NA)) +
  facet_grid(LABEL ~ .)
Phist
ggsave(Phist, filename = "SNP_distribution_histogram.pdf", width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5)
ggsave(Phist, filename = "SNP_distribution_histogram.png", width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5, dpi = 500)

# 欧氏距离（Euclidean Distance）/ED 欧氏距离是用于度量多维空间中点之间的距离的指标，通常用于比较不同位点的等位基因比例之间的差异。这里，它计算了两个维度的差异，即参考等位基因和变异等位基因的比例。
hh <- hh %>% mutate(WT.ref.rate = WT.ref.depth / WT.sum,
                    WT.alt.rate = WT.alt.depth / WT.sum,
                    mut.ref.rate = mut.ref.depth / mut.sum,
                    mut.alt.rate = mut.alt.depth / mut.sum,
                    ED = sqrt((WT.ref.rate - mut.ref.rate)^2 + (WT.alt.rate - mut.alt.rate)^2),  # sqrt 是数学函数，代表 "square root"（平方根）。它用于计算一个数字的平方根
                    ED4 = ED^4)

w <- winScan(x = hh,
             groups = "CHROM",  # 这指定了按照 "CHROM" 列的值对数据进行分组。窗口扫描将在不同的染色体（"CHROM" 列的不同值）上执行。
             position = "POS",
             values = c("ED", "ED4"), 
             win_size = 500000,
             win_step = 250000,
             funs = c("mean"))  # funs = c("mean"): 这是你要应用于每个窗口的统计函数，这里指定了 "mean"，表示计算每个窗口中的平均值。

w <- w %>% select(CHROM, win_start, win_end, win_mid, ED = ED_mean, ED4 = ED4_mean, SNP_number = ED_n)
write_tsv(file = "BSR_ED.txt", x = w)
write_csv(file = "BSR_ED.csv", x = w)

d <- chromColor %>% left_join(w, by = "CHROM")  # left_join 操作将在 "CHROM" 列上匹配两个数据框，并将匹配的行合并到一个新的数据框 d 中。连接两个数据框
d
P1 <- ggplot(filter(d, SNP_number > 10), aes(x = win_mid, y = ED)) +
  geom_point(aes(color = COLOR), size = 0.7) +
  #ylim(-1, 1) +
  labs(x="", y="ED") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
P1
ggsave(P1, filename = "BSR_ED.pdf", width = 15, height = 4)
ggsave(P1, filename = "BSR_ED.png", width = 15, height = 4, dpi = 500)

P2 <- ggplot(filter(d, SNP_number > 10), aes(x = win_mid, y = ED4)) +
  geom_point(aes(color = COLOR), size = 0.7) +
  #ylim(-1, 1) +
  labs(x="", y="ED4") +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  scale_color_aaas() +
  theme_cowplot() +
  theme(legend.position = "none") +
  #mytheme +
  facet_grid(. ~ LABEL, scales = "free_x", space = "free_x")
P2
ggsave(P2, filename = "BSR_ED4.pdf", width = 15, height = 4)
ggsave(P2, filename = "BSR_ED4.png", width = 15, height = 4, dpi = 500)

###
# chrA10
P_A10_1 <- ggplot(filter(d, SNP_number > 10, CHROM == "scaffoldA10"), aes(x = win_mid, y = ED)) +
  geom_point(color = "orange", size = 1) +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0M", "5M", "10M", "15M", "20M", "25M")) +
  labs(x = "Position", y = "ED") +
  theme_half_open()
P_A10_1
ggsave(P_A10_1, filename = "A10_ED.pdf", height = 3, width = 4)
ggsave(P_A10_1, filename = "A10_ED.png", height = 3, width = 4, dpi = 500)

P_A10_2 <- ggplot(filter(d, SNP_number > 10, CHROM == "scaffoldA10"), aes(x = win_mid, y = ED4)) +
  geom_point(color = "orange", size = 1) +
  scale_x_continuous(breaks = c(0, 5000000, 10000000, 15000000, 20000000, 25000000),
                     labels = c("0M", "5M", "10M", "15M", "20M", "25M")) +
  labs(x = "Position", y = "ED4") +
  theme_half_open()
P_A10_2
ggsave(P_A10_2, filename = "A10_ED4.pdf", height = 3, width = 4)
ggsave(P_A10_2, filename = "A10_ED4.png", height = 3, width = 4, dpi = 500)



x <- hh %>% filter((WT.geno == "0/0" & mut.geno == "1/1") | (WT.geno == "1/1" & mut.geno == "0/0"))
x <- chromColor %>% left_join(x, by = "CHROM")
ggplot(x, aes(x = POS)) + 
  geom_histogram() + 
  facet_grid(LABEL ~ ., scales = "free_x", space = "free_x")
