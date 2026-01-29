#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)

gene="VCAN"              #目标基因的名称
expFile="symbol.txt"     #表达输入文件
logFCfilter=1            #logFC过滤条件
fdrFilter=0.05           #fdr过滤条件

#读取文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
rt=t(avereps(data))

#按照基因表达的中位值对样品分组
low=rt[gene,]<=median(rt[gene,])
high=rt[gene,]>median(rt[gene,])
lowRT=rt[,low]
highRT=rt[,high]
conNum=ncol(lowRT)        #低表达组样品数目
treatNum=ncol(highRT)     #高表达组样品数目
data=cbind(lowRT,highRT)
Type=c(rep(1,conNum), rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

#输出差异的结果
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="diff.txt", sep="\t", row.names=F, quote=F)

#绘制差异基因的热图
geneNum=50    #设置显示基因的数目
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()

# 筛选热图显示的基因
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName[1:min(diffLength, 2*geneNum)]
}

# 检查hmGene是否存在于data中
hmGene=hmGene[hmGene %in% rownames(data)]
cat("选择的热图基因数量:", length(hmGene), "\n")

if(length(hmGene)>0){
  # 提取热图数据
  hmData=data[hmGene,,drop=F]

  # 检查数据是否有问题
  if(any(is.na(hmData)) || any(is.infinite(hmData))){
    cat("警告: 数据包含NA或Inf值，进行清理\n")
    hmData[is.na(hmData) | is.infinite(hmData)] = 0.01
  }

  # 确保所有值都为正数再取log
  hmData[hmData <= 0] = 0.01
  hmExp=log2(hmData + 0.01)

  # 再次检查hmExp是否有问题
  if(any(is.na(hmExp)) || any(is.infinite(hmExp))){
    cat("警告: log转换后仍有问题，使用原始数据\n")
    hmExp=hmData
  }

  Type=c(rep("Low",conNum),rep("High",treatNum))
  Type=factor(Type, levels=c("Low","High"))
  names(Type)=colnames(hmData)
  Type=as.data.frame(Type)

  # 定义绘制热图的函数
  drawHeatmap <- function() {
    tryCatch({
      pheatmap(hmExp,
               annotation=Type,
               color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
               cluster_cols =F,
               show_colnames = F,
               scale="row",
               fontsize = 8,
               fontsize_row=5,
               fontsize_col=8)
    }, error=function(e){
      cat("pheatmap出错，使用简化参数重试:", e$message, "\n")
      # 使用简化参数重试
      pheatmap(hmExp,
               annotation=Type,
               color = colorRampPalette(c("blue", "white", "red"))(50),
               cluster_cols =F,
               show_colnames = F,
               scale="none",  # 改为不标准化
               fontsize = 8,
               fontsize_row=5,
               fontsize_col=8)
    })
  }

  # 保存为PDF格式
  pdf(file="heatmap.pdf", width=10, height=7)
  drawHeatmap()
  dev.off()
  cat("热图已保存到 heatmap.pdf\n")

  # 保存为PNG格式 (高分辨率)
  png(file="heatmap.png", width=2000, height=1400, res=200)
  drawHeatmap()
  dev.off()
  cat("热图已保存到 heatmap.png\n")

  # 保存为JPEG格式
  jpeg(file="heatmap.jpg", width=2000, height=1400, quality=90)
  drawHeatmap()
  dev.off()
  cat("热图已保存到 heatmap.jpg\n")

  # 保存为TIFF格式 (高质量)
  tiff(file="heatmap.tiff", width=2000, height=1400, res=200)
  drawHeatmap()
  dev.off()
  cat("热图已保存到 heatmap.tiff\n")
}else{
  cat("错误: 没有找到合适的差异基因来绘制热图\n")
}
