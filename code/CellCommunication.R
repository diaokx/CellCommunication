rm(list = ls())  # 清除所有对象

library(Matrix)

# 将工作空间设置到此处,后面用到的文件都在此空间中（根据自己需求设置）
setwd("D:/personal_file/Bioinfomatics Study/data")
getwd()
# [1] "D:/personal_file/Bioinfomatics Study/data"  


# 读入数据并制作表达矩阵
matrix_data <- Matrix::readMM("matrix.txt")
print(object.size(matrix_data), unit="MB")
new_matrix_data <- Matrix::as.matrix(matrix_data)
col <- read.table("barcodes.txt")
col <- as.character(col$V1)
row <- read.table("features.txt")
row <- as.character(row$V1)
row <- toupper(row)
dimnames(new_matrix_data)=list(row,col)

exp <- as.data.frame(new_matrix_data)  # 转化为数据框便于操作
print(object.size(exp), unit="GB")
rm(new_matrix_data,matrix_data)
gc()


# 筛选出配受体基因并整理exp
LRgene <- read.table("geneOfpairs.txt", col.names = "Gene")
exp <- exp[rownames(exp) %in% LRgene$Gene,]


# 计算NK、Tumor和MDSC类型细胞的基因表达水平
ann <- read.table("mdsc_tumor_nk_celltype.txt", head=T, row.names = 1)
MDSC_c <- ann$cell[ann$celltype=="MDSC"]  #MDSC类型细胞
NK_c <- ann$cell[ann$celltype=="NK"]
Tumor_c <- ann$cell[ann$celltype=="Tumor"]
MDSC_Gene_Mean <- rowMeans(exp[MDSC_c])    #计算MDSC类型细胞的基因表达水平
NK_Gene_Mean <- rowMeans(exp[NK_c])
Tumor_Gene_Mean <- rowMeans(exp[Tumor_c])

geneMean <- cbind(row.names(exp), MDSC_Gene_Mean, NK_Gene_Mean, Tumor_Gene_Mean)
rownames(geneMean) <- seq(1, length(row.names(geneMean)))
colnames(geneMean) <- c("gene", "MDSC", "NK", "Tumor")
geneMean <- as.data.frame(geneMean)

# cal0 <- function(data){  # 计算表达相应基因的比例
#   x <- data[data!=0]
#   p <- length(x) / length(data)
#   return(p)
# }   #计算向量中非零的比率
# calGenePro <- function(data){
#   rowname <- rownames(data)  #数据的行名
#   ret <- apply(data, 1, cal0)
#   names(ret) <- rowname
#   return(ret)
# }
# MDSC_Gene_Pro <- calGenePro(exp[MDSC_c])   #MDSC类型细胞中表达相应基因的比例
# NK_Gene_Pro <- calGenePro(exp[NK_c])
# Tumor_Gene_Pro <- calGenePro(exp[Tumor_c])
# geneMean$MDSC <- round(MDSC_Gene_Mean*(ifelse(MDSC_Gene_Pro >= 0.1, 1, 0)), 3)
# geneMean$NK <- round(NK_Gene_Mean*(ifelse(NK_Gene_Pro >= 0.1, 1, 0)), 3)
# geneMean$Tumor <- round(Tumor_Gene_Mean*(ifelse(Tumor_Gene_Pro >= 0.1, 1, 0)), 3)

write.table(geneMean, "./output/gene_means.txt", sep = "\t", quote = F, row.names = F)


# 计算各个配受体对的平均表达量
LRpairs <- read.table("L-Rpairs.txt", header = T)
noLR <- LRgene$Gene[!(LRgene$Gene %in% geneMean$gene)]
LRpairs <- LRpairs[!(LRpairs$Ligand %in% noLR),]  # 去除用不到的LRpairs
LRpairs <- LRpairs[!(LRpairs$Receptor %in% noLR),]

typePairs <- c("Tumor-NK","Tumor-MDSC", "NK-Tumor", "NK-MDSC", "MDSC-NK", "MDSC-Tumor")
temp <- c(0, 0)                                  # 临时存储配受体表达量
names(temp) <- paste("celltype", 1:2, sep = "")
means <- matrix(0, nrow = length(rownames(LRpairs)), ncol = 7, 
                dimnames = list(1:length(rownames(LRpairs)), c("LR", typePairs)))
for (p in 1:length(rownames(LRpairs))) {  # 计算平均表达量
  pair <- as.character(LRpairs[p,])
  n <- 2
  for (typePair in typePairs) {
    typePair <- strsplit(typePair, split = "-")[[1]]
    for (i in 1:length(pair)) {
      index <- grep(pair[i], geneMean$gene, fixed = TRUE)
      temp[i] <- geneMean[index, typePair[i]]
    }
    temp <- as.numeric(temp)
    expLevel <- ifelse(!(temp[1] & temp[2]), 0, (temp[1]+temp[2])/2)
    means[p, n] <- round(expLevel, 3)
    n <- n + 1
  }
  means[p, 1] <- paste(pair[1], pair[2], sep = "-")
}
means <- as.data.frame(means)
write.table(means, "./output/means.txt", sep = "\t", quote = F, row.names = F)


# 计算各个配受体的p-values
p_values <- matrix(0, nrow = 129, ncol = 6, 
                   dimnames = list(1:length(rownames(LRpairs)), typePairs))
for (i in 1:1000) {
  ann_Random <- cbind(sample(ann$celltype, length(ann$celltype), replace = F), ann$cell)  # shuffle
  ann_Random <- as.data.frame(ann_Random)
  colnames(ann_Random) <- c("celltype", "cell")
  
  MDSC_c_Random <- ann_Random$cell[ann_Random$celltype=="MDSC"]  #MDSC类型细胞
  NK_c_Random <- ann_Random$cell[ann_Random$celltype=="NK"]
  Tumor_c_Random <- ann_Random$cell[ann_Random$celltype=="Tumor"]
  MDSC_Gene_Mean_Random <- round(rowMeans(exp[MDSC_c_Random]), 3)    #计算MDSC类型细胞的基因表达水平
  NK_Gene_Mean_Random <- round(rowMeans(exp[NK_c_Random]), 3)
  Tumor_Gene_Mean_Random <- round(rowMeans(exp[Tumor_c_Random]), 3)
  
  geneMean_Random <- cbind(row.names(exp), MDSC_Gene_Mean_Random, NK_Gene_Mean_Random, Tumor_Gene_Mean_Random)
  rownames(geneMean_Random) <- seq(1, length(row.names(geneMean_Random)))
  colnames(geneMean_Random) <- c("gene", "MDSC", "NK", "Tumor")
  geneMean_Random <- as.data.frame(geneMean_Random)
  
  means_Random <- matrix(0, nrow = length(rownames(LRpairs)), ncol = 7, 
                         dimnames = list(1:length(rownames(LRpairs)), c("LR", typePairs)))
  for (p in 1:length(rownames(LRpairs))) {  # 计算平均表达量
    pair <- as.character(LRpairs[p,])
    n <- 2
    for (typePair in typePairs) {
      typePair <- strsplit(typePair, split = "-")[[1]]
      for (i in 1:length(pair)) {
        index <- grep(pair[i], geneMean_Random$gene, fixed = TRUE)
        temp[i] <- geneMean_Random[index, typePair[i]]
      }
      temp <- as.numeric(temp)
      expLevel <- ifelse(!(temp[1] & temp[2]), 0, (temp[1]+temp[2])/2)
      means_Random[p, n] <- round(expLevel, 3)
      n <- n + 1
    }
    means_Random[p, 1] <- paste(pair[1], pair[2], sep = "-")
  }
  means_Random <- as.data.frame(means_Random)
  
  TEMP <- means_Random[2:7] > means[2:7]
  p_values <- p_values + TEMP
}
p_values <- p_values/1000  # 计算P值
p_values <- as.data.frame(p_values)
p_values$`Tumor-NK`[means$`Tumor-NK`==0] <- 1.0
p_values$`Tumor-MDSC`[means$`Tumor-MDSC`==0] <- 1.0
p_values$`NK-Tumor`[means$`NK-Tumor`==0] <- 1.0
p_values$`NK-MDSC`[means$`NK-MDSC`==0] <- 1.0
p_values$`MDSC-NK`[means$`MDSC-NK`==0] <- 1.0
p_values$`MDSC-Tumor`[means$`MDSC-Tumor`==0] <- 1.0
p_values$LR <- means$LR
p_values <- p_values[c(7, 1, 2, 3, 4, 5, 6)]
write.table(p_values, "./output/p-values.txt", sep = "\t", quote = F, row.names = F)


# 制作significant_means.txt  阈值为0.05
significant_means <- means
significant_means$`Tumor-NK` <- as.numeric(p_values$`Tumor-NK` <=0.05) * as.numeric(significant_means$`Tumor-NK`)
significant_means$`Tumor-MDSC` <- as.numeric(p_values$`Tumor-MDSC` <=0.05) * as.numeric(significant_means$`Tumor-MDSC`)
significant_means$`NK-Tumor` <- as.numeric(p_values$`NK-Tumor` <=0.05) * as.numeric(significant_means$`NK-Tumor`)
significant_means$`NK-MDSC` <- as.numeric(p_values$`NK-MDSC` <=0.05) * as.numeric(significant_means$`NK-MDSC`)
significant_means$`MDSC-NK` <- as.numeric(p_values$`MDSC-NK` <=0.05) * as.numeric(significant_means$`MDSC-NK`)
significant_means$`MDSC-Tumor` <- as.numeric(p_values$`MDSC-Tumor` <=0.05) * as.numeric(significant_means$`MDSC-Tumor`)
write.table(significant_means, "./output/significant_means.txt", sep = "\t", quote = F, row.names = F)


# 绘制circos弦图
library(circlize)
data <- matrix(0, nrow = nrow(significant_means)*(ncol(significant_means)-1), ncol = 3, 
               dimnames = list(rep(significant_means$LR, each = ncol(significant_means)-1), 
                               c("cell1", "cell2", "relation")))

t <- 1
for (i in 1:nrow(significant_means)) {  # 为data填写数据
  n <- 2
  for (typePair in typePairs) {
    typePair <- strsplit(typePair, split = "-")[[1]]
    data[t, 1] <- typePair[1]
    data[t, 2] <- typePair[2]
    data[t, 3] <- significant_means[i, n]
    n <- n + 1
    t <- t + 1
  }
}
data <- as.data.frame(data)
data$relation <- as.numeric(data$relation)

grid.col = NULL  # 开始绘制
grid.col[c("MDSC", "NK", "Tumor")] <- c("red", "grey", "grey")
visible <- rep(FALSE, nrow(data))
visible[grep("CD274.PDCD1", rownames(data), fixed = TRUE)] <- TRUE

pdf(file="circlize.pdf", width=8, height=8, pointsize=8)
chordDiagram(data, grid.col = grid.col, order = c("NK", "MDSC", "Tumor"), 
             directional = 1, link.visible = visible, 
             direction.type = c("diffHeight", "arrows"), 
             link.arr.type = "big.arrow", diffHeight = uh(0.05, "mm"), 
             annotationTrack = c("name", "grid"), 
             link.sort = TRUE, link.decreasing = FALSE)
title(main = "CD274-PDCD1", cex.main = 1.8)
dev.off()

