#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("corrplot")
#install.packages("circlize")


#引用包
library(limma)
library(circlize)
library(corrplot)

expFile="symbol.txt"         #表达数据文件
corFile="corResult.txt"      #共表达的结果文件

#读取输入文件，对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
rt=data[rowMeans(data)>0,]

#提取共表达基因的表达量
corRT=read.table(corFile, header=T, sep="\t", check.names=F)
corRT=corRT[order(corRT$cor),]
corRT=corRT[c(1:5,(nrow(corRT)-5):nrow(corRT)),]
corGene=unique(c(corRT[,1],corRT[,2]))
data=t(rt[corGene,])
cor1=cor(data)

#设置图形颜色
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(data))

#绘制圈图
pdf(file="circos.pdf", width=6.5, height=6.5)
par(mar=c(2,2,2,4))
circos.par(gap.degree =c(3,rep(2, nrow(cor1)-1)),start.degree = 180)
chordDiagram(cor1,grid.col=rainbow(ncol(data)),transparency = 0.5,col=col1,symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))       #绘制图例
dev.off()
circos.clear()
