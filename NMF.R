# Source-code archive: read KNOWN_LIMITATIONS.md before running.
# Scientific settings retained from the supplied script; not validated against study data.
library(survival)
library(survminer)
library(NMF)
library("doMPI")
rm(list = ls())
# Run from a working directory containing the required input files; see README.md.
dat <- read.csv("sd.csv",header = T,row.names = 1,check.names = F)
x.rank <- nmf(dat,2:5,nrun=200,seed=2583,method = 'brunet') 
summary(x.rank)
pdf('rank_survey.pdf',width = 7,height = 6,onefile = F)
plot(x.rank)
dev.off() 
#plot(x.rank) #选择最佳rank值，最大值或波动最大处的最高点

coph <- x.rank$measures$cophenetic


coph_diff <- NULL
for (i in 2:length(coph))
{coph_diff <- c(coph_diff, abs(coph[i-1]-coph[i]))}

k.best <- which.max(coph_diff)+1


res.improve <- nmf(dat, k.best, nrun = 200, seed = 2583)
# coefmap from multiple run fit: includes a consensus track
#coefmap(res.improve)

library(doMPI)

res.improve <- nmf(dat,k.best , nrun = 200, seed = 2583)
pdf('coefmap.pdf',width = 8,height = 8,onefile = F)
coefmap(res.improve,
        annRow = NA,
        annCol = NA,
        main = "Metagene contributions in each sample",
        info = FALSE)
dev.off() 
pdf('basismap.pdf',width = 8,height = 8,onefile = F)
basismap(res.improve,
         annRow = NA,
         annCol = NA,
         main = "Metagenes",
         info = FALSE)
dev.off() 
pdf('consensusmap.pdf',width = 8,height = 8,onefile = F)
consensusmap(res.improve,
             annRow = NA,
             annCol = NA,
             main = "Consensus matrix",
             info = FALSE)
dev.off() 
res.multi.method <- nmf(dat, 3, list('brunet', 'lee', 'ns'), seed=2583)
compare(res.multi.method)
##查看分组情况
group <- predict(res.improve)
group <- as.data.frame(group)
group$group <- paste0('Cluster',group$group)
group$sample <- rownames(group)
group<- group[order(group$group),]
table(group$group)
head(group)
#保存分组情况
#save(group,file = 'expr_group.rdata')  #建立稳定模型后，可以直接打开rdata对数据进行分类分析
write.csv(group,'expr_group2.csv') #也可以以csv格式导出到本地，查看具体分类情况

A <- read.csv("expr_group2.csv",check.names = F)
B <- read.csv("总生存率.csv",check.names = F)
A <- A[,-1]
colnames(A) <- c("group","ID")
data<- merge(A,B,by.x='ID',by.y = 'ID')
write.csv(data,"redult.csv")
data <- read.csv("redult.csv")
data$Status <- as.numeric(data$Status)
data$days <- as.numeric(data$days)
data$group <- as.factor(data$group)
survdata <- Surv(time = data$days,            #生存时间数据
                 event = data$Status)
#data <- as.data.frame(data)

#data <- as.data.frame(data[,c(1,2,4,5)])
fit_Surgery <- survfit(Surv(days,Status) ~ group,data = data) 
ggsurvplot(fit_Surgery, 
           data = data, 
          pval.method=T,#添加计算p值的统计方法
           pval = TRUE,#添加P值
           #risk.table = TRUE,#添加风险表
           xlab = "Follow up time(Months)",
           legend = c(0.8,0.15), #指定风险表图例位置
           break.x.by = 12) #设置X轴刻度间距


tiff(filename = "w.tif",width=12,height=12,
     units="cm",res=1000,pointsize=12)
dev.off()
res.improve <- nmf(dat, 4, nrun = 10, seed = 123)

