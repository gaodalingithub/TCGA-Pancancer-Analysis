#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

# 参数设置
pvalueFilter=0.05      #p值过滤条件
qvalueFilter=0.05      #矫正后的p值过滤条件

# 检查参数是否定义
if(!exists("pvalueFilter")) pvalueFilter=0.05
if(!exists("qvalueFilter")) qvalueFilter=0.05

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

# 检查是否有富集结果
if(nrow(KEGG)==0){
  stop("错误: 没有找到显著富集的KEGG通路，请检查pvalueFilter和qvalueFilter设置")
}

cat("找到", nrow(KEGG), "个显著富集的KEGG通路\n")
cat("将显示前", showNum, "个通路\n")

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
# 使用suppressWarnings抑制可能的参数警告
bar <- suppressWarnings(
  barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, color=colorSel)
)
print(bar)
dev.off()

#气泡图
pdf(file="bubble.pdf", width = 9, height = 7)
# 使用suppressWarnings抑制可能的参数警告
bub <- suppressWarnings(
  dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, color=colorSel)
)
print(bub)
dev.off()
