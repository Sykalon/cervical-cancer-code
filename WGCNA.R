# Source-code archive: read KNOWN_LIMITATIONS.md before running.
# Scientific settings retained from the supplied script; not validated against study data.
library("WGCNA")
options(stringsAsFactors = FALSE)
rm(list=ls())
# Run from a working directory containing the required input files; see README.md.
femData <- read.csv("data.csv",header = T,row.names = 1,check.names = F)
dim(femData)	# dim查看矩阵形状
#行为基因名，列为样本名
names(femData)	# names查看列标，即每一列的标题

tData <- t(femData)
##### 软阈值β Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(tData, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
pdf("beita.pdf",height = 5,width = 7.5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels=powers, cex=cex1,col="red")
dev.off()
##计算power值
power_select=sft[['powerEstimate']]
power_select  #这里最佳软阈值是5，可以通过代码得到
cor <- WGCNA::cor
net = blockwiseModules(tData,power= power_select,	 # 表达矩阵，软阈值
                       TOMType ="unsigned", minModuleSize = 50,	# 数据为无符号类型，最小模块大小为50
                       reassignThreshold = 0, mergeCutHeight = 0.15,	#mergeCutHeight合并模块的阈值，越大模块越少
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase ="TOM",
                       verbose = 3)
table(net$colors)
# 可视化模块
# Interactive graphics window omitted for command-line execution.
# 将标签转换为颜色
moduleColors = labels2colors(net$colors,zeroIsGrey=T)
# 绘制树状图和模块颜色图
pdf("module.pdf",width = 10,height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off() 
table(net$color)#查看有几个模块

moduleLabels = net$colors # 获取基因模块的代表编号
moduleColors = labels2colors(moduleLabels,zeroIsGrey=T) # 获取基因模块的代表颜色

####分析####
ADJ=adjacency(tData,power = power_select) ## 计算每个基因之间的相关性
vis=exportNetworkToCytoscape(ADJ,edgeFile = 'edge1.txt',nodeFile = 'node.txt',threshold = 0.4)
library("stringr")
MEs=net$MEs
MEs_col=MEs
colnames(MEs_col)=paste0('ME',labels2colors(
  as.numeric(str_replace_all(colnames(MEs),'ME',''))))
MEs_col=orderMEs(MEs_col)
### 各模块之间的相关性热图

pdf("模块间相关性热图.pdf",width = 10,height = 10)
plotEigengeneNetworks(MEs_col,'Eigengene adjacency heatmap',
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2),plotDendrograms = T,
                      xLabelsAngle=90)
dev.off() 

#与表型结合

traits <- read.csv("trait.csv",header = T,row.names = 1,check.names = F)
# 计算eigengenes
me_list <- moduleEigengenes(tData, colors = moduleColors)
mes <- me_list$eigengenes

n_genes <- ncol(tData)
n_samples <- nrow(tData)

module_trait_cor <- cor(mes, traits, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor,n_samples)
text_matrix <- paste(signif(module_trait_cor, 2), "\n(", signif(module_trait_pvalue, 1), ")", sep = "")
dim(text_matrix) = dim(module_trait_cor)

pdf("heat2.pdf",width = 8,height =10)

par(mar = c(3, 8.5, 3, 3))
labeledHeatmap(Matrix = module_trait_cor, xLabels = names(traits), 
               yLabels = names(mes), ySymbols = names(mes), 
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = text_matrix, 
               xLabelsAngle = 0,
               setStdMargins = FALSE, cex.text = 0.8,
               zlim = c(-1,1))

dev.off()

# 计算KME值  模块隶属度
kme <- signedKME(tData, mes, outputColumnName="kME")
# 提取感兴趣的module，这里以magenta为例
filtered <- abs(as.numeric(kme$kMEgreen)) > 0.8
               hubgene_magenta <- rownames(kme)[filtered]
hubgene_magenta
# Original console output: [1] "BAG3"    "BICD2"   "CSNK1A1" "DSG3"   
# Original console output: [5] "DSP"     "FABP5"   "JUP"     "KRT16"  
# Original console output: [9] "KRT5"    "KRT6A"   "MICALL1" "SFN"    
# Original console output: [13] "TMEM40"  "ZNF185" 
table(net$color)
####结束####
#数据导出
wgcna_result <- data.frame(gene_id = names(net$colors),module = net$colors,color=moduleColors)
write.csv(wgcna_result, file="WGCNA_result.csv",row.names = F,quote = F)

#输出到Cytoscape
TOM = TOMsimilarityFromExpr(tData, power = 4);
annot = read.csv(file = "WGCNA_result.csv");
#选择模块
modules = "greenyellow"
probes <- colnames(tData)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# 选择相关的 TOM矩阵
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


MEs <- net$MEs # 修改名称
colnames(MEs) <- str_remove(colnames(MEs),"ME")
save(MEs,design,file="df.Rdata") # 结果保存
