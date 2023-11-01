# TO compare the correlation 
COSM script of YASARA generate different results within two server with same input and yasara.ini

Marco have done CPU and GPU calculation, and all results differ (as expected for simulated annealing), but show in general a sufficient correlation of bad and good variants.

Try to figure out the problem.

I think the result is acceptable, at least not bad, only a few outliers are challenging.

```R
# AUTHOR: zhu
# DATE:20203.11.1


# load package and data
setwd('/home/zhu/save/cosmcor/')
dat = read.table('cov.tab', sep='\t', header = T)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(ggpubr)


# Scatterplot
theme_set(theme_classic())  # pre-set the bw theme.


g1 = ggplot(dat, aes(PSE_earth, PSE_RTX1)) + 
  geom_count(alpha=0.7) + 
  geom_smooth(method="lm", se=F) +  theme_classic() + theme(legend.position = 'none') + theme(panel.grid = element_blank()) +
  stat_cor(method = 'pearson', label.x=0, label.y=1500)

p1 = ggMarginal(g1, type = "histogram", fill="transparent")
# ggMarginal(g, type = "boxplot", fill="transparent")

g2 = ggplot(dat, aes(LSE_earth, LSE_RTX1)) + 
  geom_count(alpha=0.7) + 
  geom_smooth(method="lm", se=F) +  theme_classic() + theme(legend.position = 'none') + theme(panel.grid = element_blank())+ stat_cor(label.x=-60,label.y=80,method = 'spearman')
p2 = ggMarginal(g2, type = "histogram", fill="transparent") 
require(gridExtra)
grid.arrange(p1,p2)
ggsave('cor.png',width = 8,height = 16)

```