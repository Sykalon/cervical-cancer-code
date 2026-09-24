# Source-code archive: read KNOWN_LIMITATIONS.md before running.
# Scientific settings retained from the supplied script; not validated against study data.
library(randomForest)
library(caret)
library(dplyr)


# Run from a working directory containing the required input files; see README.md.
rm(list = ls())
data <- read.csv("scale.csv")
X <- data %>% select(-Group)  # 特征（蛋白质数据）
y <- data$Group              # 标签（三组分类）
new_data <- read.csv("data.csv",header = T)
set.seed(656)
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_idx, ]
y_train <- y[train_idx]
X_test <- X[-train_idx, ]
y_test <- y[-train_idx]
data_train <- data[train_idx, ]
# 训练随机森林
rf_model <- randomForest(
  x = X_train,
  y = as.factor(y_train),  # 必须将标签转为因子
  ntree = 800,             # 树的数量
  mtry = floor(sqrt(ncol(X_train))),  # 每棵树随机选择的特征数（默认值）
  importance = TRUE        # 计算特征重要性
)


# 预测测试集
y_pred <- predict(rf_model, X_test)

# 计算准确率
accuracy <- mean(y_pred == y_test)
cat("测试集准确率:", accuracy, "\n")



importance_scores <- importance(rf_model, type = 2)
mean_importance <- mean(importance_scores[, "MeanDecreaseGini"])
selected_features <- which(importance_scores[, "MeanDecreaseGini"] > mean_importance)
cat("筛选出的变量数量:", length(selected_features), "\n")
cat("变量名称:", names(selected_features))

# 按重要性降序排序
sorted_features <- importance_scores[order(importance_scores[, "MeanDecreaseGini"], decreasing = TRUE), ]

# 选择前5个特征（原始代码 top_n = 5）
top_n <- 5
selected_features2 <- names(sorted_features)[1:top_n]
# 加载新数据（假设为 new_data.csv）


selected_features<- selected_features2[selected_features2 %in% colnames(new_data)]
#cat("筛选出的变量数量:", top_n, "\n")
#cat("变量名称:",selected_features,"\n")

# 提取筛选后的特征数据
X_train_selected <- X_train[, selected_features]

# 重新训练模型
set.seed(184)
rf_model_selected <- randomForest(x = X_train_selected, 
                                  y = as.factor(y_train), 
                                  ntree = 800)

# 评估性能（对比原始模型）
y_pred_selected <- predict(rf_model_selected, X_test[, selected_features])
accuracy_selected <- mean(y_pred_selected == y_test)
accuracy_selected
cat("筛选后模型准确率:", accuracy_selected, "\n")

#新数据验证
predicted_groups <- predict(rf_model_selected, new_data[, selected_features])
names(predicted_groups) <- new_data$name
#cat("新数据的预测组别:", predicted_groups, "\n")
write.csv(predicted_groups,"predicted_groups.csv")
A1 <- read.csv("predicted_groups.csv")
A2 <- read.csv("PFS.csv")
colnames(A1)[1] <- "name"
colnames(A1)[2] <- "group"
re<- merge(A1,A2,by.x='name',by.y = 'name')
library(survival)
library(survminer)

surv_obj <- Surv(time = re$Days, event = re$Status)
fit <- survfit(surv_obj ~ group, data = re)
logrank_test <- survdiff(surv_obj ~ group, data = re)

# 绘制生存曲线
ggsurvplot(
  fit,
  data = re,
  pval = TRUE,          # 显示 Logrank 检验的 p值
  pval.method = TRUE,   # 显示检验方法
  conf.int = FALSE,     # 是否显示置信区间
  risk.table = F,    # 显示风险表
  #legend.title = "Group",
  #legend.labs = c("S-I", "S-II", "S-III"),  # 确保标签顺序正确
  palette = c("#3179B5", "#F07A00","#D02C1B")    # 自定义颜色
)
accuracy_selected
ggplot2::ggsave("km.png", plot = km_plot$plot, width = 5, height = 4)  ##ggplot 中直接保存
#两两比较（示例：Group1 vs Group2）
#pairwise_test1 <- survdiff(Surv(Days, Status) ~ group, data = re, subset = group %in% c("S_I", "S_II"))
#pairwise_test1
#pairwise_test2 <- survdiff(Surv(Days, Status) ~ group, data = re, subset = group %in% c("S_I", "S_III"))
#pairwise_test2
#pairwise_test3 <- survdiff(Surv(Days, Status) ~ group, data = re, subset = group %in% c("S_II", "S_III"))
#pairwise_test3
#dev.off()
# 输出概率（可选）
#probabilities <- predict(rf_model_selected, new_data, type = "prob")
#print("属于每组的概率:")
#print(probabilities)


# 定义特征选择后的数据
data_selected <- data[, c(selected_features, "Group")]

#使用交叉验证训练模型
#ctrl <- trainControl(method = "cv", number = 10)
#model_cv <- train(Group ~ ., data = data_selected, 
#method = "rf", trControl = ctrl)
#best <- model_cv[["bestTune"]][["mtry"]]
#print(model_cv)
#rf_model <- randomForest(
# x = X_train_selected, 
#y = as.factor(y_train),
# mtry = best,  # 手动指定
# ntree = 100
#)


print(accuracy)
print(accuracy_selected)
print(logrank_test[["pvalue"]])



# 打印混淆矩阵和分类报告
confusionMatrix(y_pred_selected, as.factor(y_test))
# 数据重塑
library(reshape2)
library(dplyr)
cm <- confusionMatrix(y_pred_selected, as.factor(y_test))
confusion_matrix_df<-as.data.frame.matrix(cm$table)
colnames(confusion_matrix_df)<-c("S-I","S-II","S-III")
rownames(confusion_matrix_df)<-c("S-I","S-II","S-III")
draw_data<-round(confusion_matrix_df/rowSums(confusion_matrix_df),2)
draw_data$real<-rownames(draw_data)
draw_data<-melt(draw_data)
library(ggprism)
P2 <- ggplot(draw_data,aes(real,variable,fill=value))+
  geom_tile()+
  geom_text(aes(label=scales::percent(value)))+
  scale_fill_gradient(low="#F0F0F0",high="#911D22")+
  labs(x="True",y="Guess",title="Confusion matrix")+
  theme_prism(border=T)+
  theme(panel.border=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

pdf("confusion_matrix.pdf",width=3.5,height=3.5,pointsize=10,onefile=FALSE)
P2
dev.off()



p1 <- ggsurvplot (fit,
                  #ggtheme=theme_bw(),#调整背景风格 https://ggplot2.tidyverse.org/reference/ggtheme.html
                  #title="",#设定标题
                  font.title=c(14,"bold","black"),#调整标题的字体大小、颜色等
                  xlab='Follow-up time (month)',#设定x轴标签
                  ylab="Survival probability (%)",#设定y轴标签
                  font.x=c(13,"plain","black"), #调整x轴的字体大小、颜色等
                  font.y=c(13,"plain","black"),#调整y轴的字体大小、颜色等
                  xlim=c(0,50),#设定x轴范围
                  ylim=c(0,1),#设定y轴范围
                  font.tickslab=10,#调整刻度的字体大小、颜色等
                  #break.x.by=365.25,#设定x轴间距,365.25代表1年
                  break.y.by=0.25,#设定y轴间距
                  # xscale='m_y',#表示将时间日转化为年，还可以"m_y"等，具体见帮助文档
                  surv.scale="percent",#以百分数表示
                  censor=T, #显示删失点
                  censor.shape=124, #删失点的形状
                  censor.size=2,#设定删失点大小
                  conf.int=F, #添加图形加置信区间
                  palette=c(#"#2c3385",
                    '#00558e', 
                    '#ca1d04',
                    '#eda163'),#自定义曲线颜色,也可选内置配色
                  pval = TRUE,#添加p值,也可以自定义,如pval="P < 0.001"
                  pval.size=5,#值字体大小
                  pval.method = TRUE, #添加检验方法,默认是log-rank
                  pval.method.size=4, #log-rank字体大小
                  pval.coord=c(0.5,0.15), #P值显示位置
                  # surv.median.line = "hv",#添加中位生存时间
                  risk.table =F, #添加风险表
                  risk.table.col="black",#风险表设定颜色
                  risk.table.height=0.2,#风险表高度（0-1之间）
                  risk.table.fontsize=4,#字体大小
                  risk.table.y.text=FALSE, #如为TRUE,风险表纵轴显示为male和female
                  risk.table.y.text.col=TRUE, #如为FALSE,风险表纵轴不显示为图例颜色
                  risk.table.pos="out",#设定风险表的位置,在生存曲线内("in")还是外("out")
                  risk.table.title="Number at risk ",#设定风险表的标题
                  tables.theme = theme_void(),##调整风险表的风格
                  legend.title="",#设定图例标题
                  legend=c(0.8,0.3),#设定图例位置
                  legend.labs=c("S-I","S-II","S-III"), #设定图例标签
                  font.legend=c(11,"plain","black"),#设定图例字体大小、颜色等
                  size = 0.8,#调整生存曲线的粗细
                  linetype="solid"##调整生存曲线
)
p1
