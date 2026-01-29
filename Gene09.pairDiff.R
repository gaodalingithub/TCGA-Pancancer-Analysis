#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

#引用包
library(limma)
library(ggplot2)
expFile="geneExp.txt"      #表达数据文件

#读取表达输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
geneName=colnames(rt)[1]

#提取配对样品的数据
normalData=rt[rt$Type=="Normal",1,drop=F]
normalData=as.matrix(normalData)
rownames(normalData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(normalData))
normalData=avereps(normalData)
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
tumorData=avereps(tumorData)
sameSample=intersect(row.names(normalData), row.names(tumorData))
data=cbind(normalData[sameSample,,drop=F], tumorData[sameSample,,drop=F])
colnames(data)=c("Normal", "Tumor")
data=as.data.frame(data)

#绘制配对差异分析的图形
pdf(file="pairDiff.pdf", width=5, height=4.5)

# 转换数据格式用于ggplot2
data_long <- data.frame(
  sample = rep(rownames(data), 2),
  condition = rep(c("Normal", "Tumor"), each = nrow(data)),
  expression = c(data$Normal, data$Tumor)
)

p <- ggplot(data_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_line(aes(group = sample), color = "gray70", alpha = 0.7) +
  geom_point(size = 2, shape = 21, color = "black") +
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = "", y = paste0(geneName, " expression"), fill = "Type") +
  theme_bw() +
  theme(legend.position = "top") +
  # 添加配对统计比较
  stat_compare_means(paired = TRUE,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                       symbols = c("***", "**", "*", "ns")),
                     label = "p.signif", label.x = 1.5)

print(p)
dev.off()
