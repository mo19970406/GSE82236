# GSE82236
#BiocManager::install('DESeq2')
#DESeq2需要输入的数值是整数
library(DESeq2)
library(GEOquery)
mirna=read.table('GSE82236_miRNA_counts.txt.gz',header=T,row.names=1,sep="\t")
mirna[1:4,1:4]
# 对所有数四舍五入取整
mirna=round(mirna,0)
head(mirna)
mirna_matrix=as.matrix(mirna)
# 进行分组
org=factor(c(rep('CC.CR',3),rep('CC',3)))
coldata=data.frame(row.names = colnames(mirna),org)

#source("http//biocondutor.org/biocLite.R")
#biocLite('DESeq2')
# 差异表达分析
dds=DESeqDataSetFromMatrix(mirna_matrix,DataFrame(coldata),design = ~org)
# 标准化
dds=DESeq(dds)
res=results(dds,alpha = 0.001)
res=res[order(res$padj),]
summary(res)
# 输出结果
write.csv(res,file = 'E:\\heatmap\\new\\results.csv')

# 绘制MA-PLOT
plotMA(res,ylim=c(-3,3))
#绘制火山图
plot(res$log2FoldChange,res$padj)
# 用ggplot2画图
library(ggplot2)
data_res <- as.data.frame(res)
data_res$color <- ifelse(data_res$padj<0.05&abs(data_res$log2FoldChange)>1,'LFC>1,p<0.05',ifelse(abs(data_res$log2FoldChange)>1,'LFC>1,p>0.05','LFC<1,p>0.05'))
ggplot(data_res,aes(log2FoldChange,padj,col=factor(color))) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  geom_vline(xintercept=c(-2,2) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("red", "green", "black"), limits = c("LFC>1,p<0.05", "LFC>1,p>0.05", "LFC<1,p>0.05")) 
# 提取差异基因并写出
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file = 'DEG.csv')

# 绘制热图

# 从mydata中提取差异最显著的前12个miRNA的表达数据
library('pheatmap')
diff_miRNA <- row.names(data_res)[1:12] # 提取前12的差异miRNA
# 选择差异最显著的前12个miRNA绘制热图
myvars <- row.names(mirna) %in% diff_miRNA
heatmap_data <- mirna[myvars,]
mymatrix <- as.matrix(heatmap_data)
annotation_col = data.frame(
  CellType = factor(c(rep("CC.CR",12),rep("CC",12)))
)
rownames(annotation_col) <- c(paste('CC.CR-',1:12,sep = ''),paste('CC-',1:12,sep = ''))
pheatmap(heatmap_data,scale = 'row', annotation_col = annotation_col, show_rownames = FALSE)


