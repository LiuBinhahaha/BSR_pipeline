# grep -v '##' BES398_filter.vcf | sed s'/#CHROM/CHROM/' > BES398_filter.txt

library(tidyverse)  # separate, filter, select...
library(cowplot)
library(ggsci)
library(RColorBrewer)
#library(devtools)
#install_github('tavareshugo/windowscanr')
library(windowscanr)

setwd("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA")

# BES236
# BES398
input <- "BES398_filter.txt"
df <- read_tsv(input)
chromColor <- read_tsv("chromColor.txt")

# 提取
df_tmp1 <- df %>% select(CHROM, POS, REF, ALT, B73, BES398, Zh58, BES398_WT, BES398_mut)
df_tmp2 <- df_tmp1 %>% separate(BES398, c("BES398.gt"), sep = ":") %>%
  separate(Zh58, c("Zh58.gt"), sep = ":") %>%
  separate(BES398_WT, c("BES398_WT.gt", "BES398_WT.depth"), sep = ":") %>%
  separate(BES398_mut, c("BES398_mut.gt", "BES398_mut.depth"), sep = ":")

# 分割
df_tmp3 <- df_tmp2 %>% separate(BES398_WT.depth, c("AD_REF.BES398_WT", "AD_ALT.BES398_WT"), sep = ",") %>%
                       separate(BES398_mut.depth, c("AD_REF.BES398_mut", "AD_ALT.BES398_mut"), sep = ",")


# 过滤掉父母本中有缺失、WT和mut池中有缺失的位点、WT和mut基因型同时是0/0的位点、WT和mut基因型同时是1/1的位点
# python脚本已过滤！！！
'''
df_tmp4 <- df_tmp3 %>% filter(!(BES398_WT.gt == '0/0' & BES398_mut.gt == '0/0')) %>%
  filter(BES398.gt != './.' & Zh58.gt != './.') %>%
  filter(BES398_WT.gt != './.' & BES398_mut.gt != './.') %>%
  filter(!(BES398_WT.gt == '1/1' & BES398_mut.gt == '1/1'))
'''

# 添加reads数sum
df_tmp5 <- df_tmp3 %>% mutate(AD_REF.BES398_WT = as.numeric(AD_REF.BES398_WT),
                              AD_ALT.BES398_WT = as.numeric(AD_ALT.BES398_WT),
                              AD_REF.BES398_mut = as.numeric(AD_REF.BES398_mut),
                              AD_ALT.BES398_mut = as.numeric(AD_ALT.BES398_mut),
                              BES398_WT.sum = AD_REF.BES398_WT + AD_ALT.BES398_WT,
                              BES398_mut.sum = AD_REF.BES398_mut + AD_ALT.BES398_mut)

# reads密度图
pdf("/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/BES398/BES398_reads_desity.pdf", height = 3, width = 6)
par(mfrow = c(1, 2))
plot(density(df_tmp5$BES398_WT.sum, width = 1), main = "BES398_WT", xlim = c(0, 100))
plot(density(df_tmp5$BES398_mut.sum, width = 1), main = "BES398_mut", xlim = c(0, 100))
dev.off()

# 统计两池的gene type
gene_type_BES398 <- df_tmp5 %>% group_by(BES398_WT.gt, BES398_mut.gt) %>% count()
write.table(gene_type_BES398, file = "/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/BES398/BES398_gene_type.tsv", row.names = F, quote = F)

# SNP_distribution_fig
options(scipen = 200)  # 在 R 中，默认情况下，当浮点数的小数部分为零时，它们将以科学计数法的形式显示。scipen 的默认值是 0，这意味着浮点数小数部分为零时将采用科学计数法。例如，0.0001 将显示为 1e-04，如果将 scipen 设置为较大的值，比如 200，那么 R 将更不容易采用科学计数法显示浮点数，即使小数部分为零，也会以普通的小数形式显示。
colourCount = dim(chromColor)[[1]]
getPalette = colorRampPalette(brewer.pal(9, "Set1"))  # brewer.pal(8, "Set1") 是使用 R 包中的 RColorBrewer 中的函数 brewer.pal，它用于生成颜色调色板。在这里，它要求生成 8 种颜色，从 "Set1" 调色板中选择颜色。

Phist <- chromColor %>% left_join(df_tmp5, by = "CHROM") %>% ggplot(aes(x = POS)) +  # 直方图。chromColor 和 hh 是两个数据框，通过 left_join 函数使用 "CHROM" 列将它们合并
  geom_histogram(aes(fill = LABEL), color = NA, binwidth = 1000000) +  # binwidth 参数设置为 1000000。这意味着 x 轴上的连续范围被分成多个宽度为 1,000,000 的柱子。这将导致数据的 x 轴范围被分为多个区间，每个区间的宽度为 1,000,000，并且直方图将显示每个区间中数据点的频数。
  labs(x = NULL, y = "SNP Count / 1Mb") +
  scale_x_continuous(expand = c(0, 0)) +  # expand = c(0, 0) 指定了 x 轴的显示范围不进行扩展
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_half_open() +
  theme(strip.text = element_text(color = NA, size = 0.1), # strip.text 和 strip.background 参数设置了分面（facet）标签的文本颜色和背景
        axis.text.y =element_text(size=5, colour = "black"), # x轴和y轴标注字体的大小
        strip.background = element_rect(color = NA, fill = NA)) +
  facet_grid(LABEL ~ .)
Phist
ggsave(Phist, filename = "/data/heqiang/labmember/Zhengjun/BSR_DNA_RNA/BES398/BES398_SNP_distribution_histogram.pdf", width = 9, height = dim(chromColor)[[1]] * 0.6 + 0.5)


