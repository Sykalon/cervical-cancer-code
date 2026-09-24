# Source-code archive: read KNOWN_LIMITATIONS.md before running.
# Scientific settings retained from the supplied script; not validated against study data.
#清空
rm(list=ls())
gc()

library(IOBR)

tme_deconvolution_methods
# Run from a working directory containing the required input files; see README.md.
expr<- read.csv("ydata.csv",header = T,row.names = 1,check.names = F)
#列名是名字，行名是基因
# MCPcounter
mcpcounter_immo <- deconvo_tme(eset = expr,
                               #project="TCGA-COAD",#项目名称
                               method = "mcpcounter"#使用方法
)
write.csv(mcpcounter_immo,"MCPcounter.csv")
# EPIC
epic_immo <- deconvo_tme(eset = expr,
                         method = "epic",
                         tumor = T)
write.csv(epic_immo,"epic.csv")
# xCell
xCell_immo <- deconvo_tme(eset = expr,
                          method = "xcell",
                          arrays = F)
write.csv(xCell_immo,"xCell.csv")
# CIBERSORT
cibersort_immo <- deconvo_tme(eset = expr,
                              method = "cibersort",
                              arrays = F,
                              perm = 1000)#设置排列数量
write.csv(cibersort_immo,"cibersort.csv")



# IPS
ips_immo <- deconvo_tme(eset = expr,
                        method = "ips",
                        plot = F)
write.csv(ips_immo,"ips.csv")
# quanTIseq
quantiseq_immo <- deconvo_tme(eset = expr,
                              method = "quantiseq",
                              scale_mrna = T)#如果为FALSE，则不能矫正不同细胞类型的mRNA含量
write.csv(quantiseq_immo,"quantiseq.csv")
# ESTIMATE
estimate_immo <- deconvo_tme(eset = expr,
                             method = "estimate")
write.csv(estimate_immo ,"estimate.csv")
# TIMER
timer_immo <- deconvo_tme(eset = expr,
                          method = "timer",
                          group_list = rep("coad",dim(expr)[2]))
write.csv(timer_immo,"timer.csv")
