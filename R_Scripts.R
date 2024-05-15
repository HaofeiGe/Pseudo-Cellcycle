rm(list = ls())
setwd("D:/360MoveData/Users/Haofei Ge/Desktop/该读书了/细胞周期伪时间课题（毕业论文版）/活细胞与伪时间的验证数据/伪时间处理/PseudoTime")

##读取数据
Nuc <- read.csv("MyExpt_Cellpose_DNA.csv")

##提取Cells和Nuc数据并修改列名
Select_parameter <- function(data){
  data <- data[c(grep("ImageNumber|ObjectNumber|AreaShape_Area|AreaShape_Perimeter|
                 |AreaShape_FormFactor|AreaShape_Compactness|AreaShape_Eccentricity|
                 |AreaShape_MajorAxisLength|AreaShape_MinorAxisLength|AreaShape_BoundingBox|
                 |AreaShape_MaxFeretDiameter|AreaShape_MinFeretDiameter|Location_Center_X|
                 |Location_Center_Y|Intensity_IntegratedIntensity|Intensity_MaxIntensity|
                 |Intensity_MeanIntensity|Intensity_MedianIntensity|Intensity_MADIntensity|
                 |Intensity_StdIntensity|Intensity_UpperQuartileIntensity|Parent_Cellpose_IF",
                      colnames(data)))]
  data <- data[-c(grep("Edge",colnames(data)))]
  colnames(data)[colnames(data)==colnames(data[c(grep("Intensity_IntegratedIntensity",colnames(data)))])] <- "Intensity_IntegratedIntensity"
  colnames(data)[colnames(data)==colnames(data[c(grep("Intensity_MaxIntensity",colnames(data)))])] <- "Intensity_MaxIntensity"
  colnames(data)[colnames(data)==colnames(data[c(grep("Intensity_MeanIntensity",colnames(data)))])] <- "Intensity_MeanIntensity"
  colnames(data)[colnames(data)==colnames(data[c(grep("Intensity_MedianIntensity",colnames(data)))])] <- "Intensity_MedianIntensity"
  colnames(data)[colnames(data)==colnames(data[c(grep("Intensity_MADIntensity",colnames(data)))])] <- "Intensity_MADIntensity"
  colnames(data)[colnames(data)==colnames(data[c(grep("Intensity_StdIntensity",colnames(data)))])] <- "Intensity_StdIntensity"
  colnames(data)[colnames(data)==colnames(data[c(grep("Intensity_UpperQuartileIntensity",colnames(data)))])] <- "Intensity_UpperQuartileIntensity"
  colnames(data)[colnames(data)==colnames(data[c(grep("Parent_Cellpose_IF",colnames(data)))])] <- "Parent_ObjectNumber"
  return(data)
}#该函数用于提取后续细胞周期伪时间的必要参数
Nuc <- Select_parameter(data = Nuc)#使用方法：Select_parameter(data = data)
colnames(Nuc) <- paste("Nuc", colnames(Nuc), sep = "_")#加上前缀

#初步过滤非正常对象
Nuc <- Nuc[Nuc[,colnames(Nuc)[grep("AreaShape_Area",colnames(Nuc))]]>=5,]

##批次校正
data <- Nuc[Nuc$Nuc_ImageNumber==1 | Nuc$Nuc_ImageNumber==2,]
#看一下批次效应
library(ggplot2)
library(ggpubr)
ggplot() + 
  geom_density(data, mapping = aes(x=Nuc_Intensity_IntegratedIntensity, 
                                   color=as.factor(Nuc_ImageNumber))) + 
  theme_classic() + ggtitle("Nuc") + labs(color="ImageNumber") + 
  theme(plot.title = element_text(hjust = 0.5))
#校正批次效应方法1，适用于Hela等拥有显著G1期峰的数据类型
#在该方法中，计算最高峰所在位置，将第一峰对齐到参考点
{
  Prop_calculate <- function(data){
    prop_result <- as.data.frame(matrix(data=NA,
                                        nrow=length(unique(data$Nuc_ImageNumber)),
                                        ncol=2))
    colnames(prop_result) <- c("Nuc_ImageNumber","Proportion")
    prop_result$Nuc_ImageNumber <- paste("Image",1:nrow(prop_result),sep = "")
    for (i in 1:length(unique(data$Nuc_ImageNumber))) {
      d <- density(data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_IntegratedIntensity"])
      t <- data.frame(x=d$x,y=d$y)
      Peak <- t[which.max(t$y),][1,1]
      prop = 2/Peak
      prop_result[i,"Proportion"] <- prop
    }
    return(prop_result)
  }#该函数计算校正因子
  Prop_Normalization <- function(data, Prop_Result){
    for (i in 1:length(unique(data$Nuc_ImageNumber))) {
      data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_IntegratedIntensity"] <- data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_IntegratedIntensity"]*Prop_Result[i,"Proportion"]
      data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MaxIntensity"] <- data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MaxIntensity"]*Prop_Result[i,"Proportion"]
      data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MeanIntensity"] <- data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MeanIntensity"]*Prop_Result[i,"Proportion"]
      data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MedianIntensity"] <- data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MedianIntensity"]*Prop_Result[i,"Proportion"]
      data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MADIntensity"] <- data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_MADIntensity"]*Prop_Result[i,"Proportion"]
      data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_StdIntensity"] <- data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_StdIntensity"]*Prop_Result[i,"Proportion"]
      data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_UpperQuartileIntensity"] <- data[which(data$Nuc_ImageNumber==i),"Nuc_Intensity_UpperQuartileIntensity"]*Prop_Result[i,"Proportion"]
    }
    return(data)
  }#该函数将校正因子应用到数据
  #校正荧光强度分布以及Mean/Median/MAX/Std/MAD/UpperQuartileIntensity
  data_correct <- data
  Prop_result <- Prop_calculate(data = data_correct)
  data_correct <- Prop_Normalization(data = data_correct, Prop_Result = Prop_result)
  #检查校正效果
  ggplot() + 
    geom_density(data_correct, mapping = aes(x=Nuc_Intensity_IntegratedIntensity, 
                                             color=as.factor(Nuc_ImageNumber))) + 
    theme_classic() + ggtitle("Nuc") + labs(color="ImageNumber") + 
    theme(plot.title = element_text(hjust = 0.5))
  }

##############################细胞周期伪时间分析###################################
{
  Exp <- data_correct
  hist(Exp$Nuc_Intensity_IntegratedIntensity, breaks = 300)
  Exp_DDR <- Exp[,c("Nuc_ImageNumber","Nuc_ObjectNumber",
                    "Nuc_Intensity_IntegratedIntensity","Nuc_AreaShape_Area",
                    "Nuc_Intensity_MADIntensity","Nuc_Intensity_MaxIntensity",
                    "Nuc_Intensity_MeanIntensity","Nuc_Intensity_MedianIntensity",
                    "Nuc_Intensity_StdIntensity","Nuc_Intensity_UpperQuartileIntensity")]
  #无量纲化
  for (i in 5:ncol(Exp_DDR)) {
    Exp_DDR[,i] <- (Exp_DDR[,i]-min(Exp_DDR[,i]))/(max(Exp_DDR[,i])-min(Exp_DDR[,i]))
  }
  #像Monocle那样走一遍
  FM <- as.matrix(Exp_DDR[,5:ncol(Exp_DDR)])
  #FM <- subset(FM,select = -c(Nuc_Intensity_DNA_IntegratedIntensity))
  FM <- t(FM)
  library(igraph)
  library(monocle)
  library(DDRTree)
  ## cal_ncenter目的是筛选出 ncenter ，即低维空间（Z）中心或者根位置的一小群细胞的数量
  cal_ncenter <- function(ncells, ncells_limit = 100){
    round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
  }
  ncenter <- cal_ncenter(ncol(FM))
  ddrtree_res <- DDRTree(FM,dimensions = 1 , ncenter = 2)
  Z <- ddrtree_res[["Z"]]
  Zt <- as.data.frame(t(Z))
  Zt <- cbind(Zt,Exp)
  ggplot(Zt, aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity,
                       color=Nuc_Intensity_UpperQuartileIntensity)) +
    geom_point() + scale_color_gradient(low = "yellow",high = "blue") + 
    xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
    labs(colour = "UpperQuartile") +
    theme(plot.title = element_text(size=18,hjust=0.5)) +
    ggtitle("Hela") + theme_classic() + ylim(1,4.8)
  
  #筛选细胞
  library(iplots)
  with(Zt, iplot(x=V1, y=Nuc_Intensity_IntegratedIntensity))
  selected_cells <- read.table("../selected.txt", header = T)
  colnames(selected_cells) <- c("V1","Nuc_Intensity_IntegratedIntensity")
  selected_cells <- merge(selected_cells, Zt, by=colnames(selected_cells))
  ggplot(selected_cells, aes(x=V1,y=Nuc_Intensity_IntegratedIntensity,
                 color=Nuc_Intensity_UpperQuartileIntensity)) +
    geom_point() + scale_color_gradient(low = "yellow",high = "blue") + 
    xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
    labs(colour = "UpperQuartile") +
    theme(plot.title = element_text(size=18,hjust=0.5)) +
    ggtitle("Hela ") + theme_classic()
  ggplot(selected_cells, aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity,
                             color=factor(Nuc_ImageNumber))) +
    geom_point() + xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
    labs(colour = "UpperQuartile") +
    theme(plot.title = element_text(size=18,hjust=0.5)) +
    ggtitle("Hela ") + theme_classic()
 # write.csv(selected_cells, file = "selected_cells.csv", row.names = F)
  }


#可视化
selected_cells <- read.csv("selected_cells.csv")
selected_cells_2 <- read.table("selected_cells.txt", header = T)
selected_cells <- merge(selected_cells, selected_cells_2, by=colnames(selected_cells_2))
live_cells <- selected_cells
selected_cells <- selected_cells[,1:26]
live_cells <- na.omit(live_cells)




chooserows <- rownames(live_cells[live_cells$Third_Origin==0,])
live_cells[chooserows,"Third_Origin"] <- live_cells[chooserows,"Second_Origin"]
live_cells[chooserows,"Second_Origin"] <- live_cells[chooserows,"First_Origin"]
live_cells[chooserows,"First_Origin"] <- 0
live_cells[chooserows,"Third_Divide"] <- live_cells[chooserows,"Second_Divide"]
live_cells[chooserows,"Second_Divide"] <- live_cells[chooserows,"First_Divide"]
live_cells[chooserows,"First_Divide"] <- "/"

live_cells[,28:35] <- as.data.frame(lapply(live_cells[,28:35], as.numeric))

library(ggplot2)
ggplot() + 
  geom_point(selected_cells, mapping = aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity),
             size = 2, alpha = 1, color = "grey") + 
  geom_point(live_cells, mapping = aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity,
                                       color=Third_Divide),
             size = 2) +
  scale_color_gradient(low = "blue",high = "yellow") +
  xlab("DDRTree axis") + labs("Absolute Time") + theme_classic()
live_cells$AbsoluteTime <- (241-live_cells$Third_Divide)*10
ggplot() + 
  geom_point(selected_cells, mapping = aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity),
             size = 2, alpha = 0.5, color = "grey") + 
  geom_point(live_cells, mapping = aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity,
                                       color=AbsoluteTime),
             size = 2) +
  scale_color_gradient(low = "yellow",high = "blue") +
  xlab("DDRTree axis") + labs("Absolute Time") + theme_classic()

live_cells$Second_cyclelength <- (live_cells$Third_Divide - live_cells$Second_Divide)*10
live_cells$RelativeTime <- live_cells$AbsoluteTime/live_cells$Second_cyclelength
live_cells$First_cyclelength <- (live_cells$Second_Divide - live_cells$First_Divide)*10
library(dplyr)
live_cells_cyclelength <- as.data.frame(matrix(NA, nrow=205*2,ncol=2))
live_cells_cyclelength$V1 <- c(live_cells$First_cyclelength,live_cells$Second_cyclelength)
live_cells_cyclelength$V2 <- c(rep("First_cyclelength",205),rep("Second_cyclelength",205))
colnames(live_cells_cyclelength) <- c("CycleLenth","group")
library(cols4all)
library(ggsignif)
mycol <- c4a('10',6)
ggplot() + 
  geom_boxplot(live_cells_cyclelength, mapping=aes(x=group, y=CycleLenth,color=group),fill=NA) + 
  geom_jitter(live_cells_cyclelength, mapping=aes(x=group, y=CycleLenth,color=group),alpha=0.5,size=3) +
  stat_boxplot(live_cells_cyclelength, mapping=aes(x=group, y=CycleLenth,color=group),
               geom="errorbar",width=0.1,size=0.8,linetype=2,) + 
  scale_color_manual(values = mycol) +
  stat_summary(fun=median,geom="point",size=3,shape=21,color="black",fill="white") +
  ylim(700,1200) + theme_classic() + xlab(NULL) + ylab("Cell Cycle Length") +
  geom_signif(comparisons =  list(c("First_cyclelength", "Second_cyclelength")))
tmp <- na.omit(live_cells_cyclelength)
t.test(tmp[tmp$group=="First_cyclelength","CycleLenth"],tmp[tmp$group=="Second_cyclelength","CycleLenth"])

ggplot() + 
  geom_point(selected_cells, mapping = aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity),
             size = 2, alpha = 1, color = "grey") + 
  geom_point(live_cells, mapping = aes(x=abs(V1),y=Nuc_Intensity_IntegratedIntensity,
                                       color=RelativeTime),
             size = 2) +
  scale_color_gradient(low = "yellow",high = "blue") +
  xlab("DDRTree axis") + labs("Relative Time") + theme_classic()


Exp_DDR <- selected_cells[,c("Nuc_ImageNumber","Nuc_ObjectNumber",
                  "Nuc_Intensity_IntegratedIntensity","Nuc_AreaShape_Area",
                  "Nuc_Intensity_MADIntensity","Nuc_Intensity_MaxIntensity",
                  "Nuc_Intensity_MeanIntensity","Nuc_Intensity_MedianIntensity",
                  "Nuc_Intensity_StdIntensity","Nuc_Intensity_UpperQuartileIntensity")]
#无量纲化
Exp_DDR[,4:ncol(Exp_DDR)] <- scale(Exp_DDR[,4:ncol(Exp_DDR)])
#像Monocle那样走一遍
FM <- as.matrix(Exp_DDR[,4:ncol(Exp_DDR)])
#FM <- subset(FM,select = -c(Nuc_Intensity_DNA_IntegratedIntensity))
FM <- t(FM)
library(igraph)
library(monocle)
library(DDRTree)
## cal_ncenter目的是筛选出 ncenter ，即低维空间（Z）中心或者根位置的一小群细胞的数量
cal_ncenter <- function(ncells, ncells_limit = 30){
  round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
}
ncenter <- cal_ncenter(ncol(FM))
ddrtree_res <- DDRTree(FM,dimensions = 2 , ncenter = ncenter)
Z <- ddrtree_res[["Z"]]
Zt <- as.data.frame(t(Z))
colnames(Zt) <- c("V1.1","V2.1")
Zt <- cbind(Zt,selected_cells)
Zt_livecells <- merge(Zt, live_cells ,by = "V1")
ggplot(Zt, aes(x=V1.1,y=V2.1,
               color=Nuc_Intensity_UpperQuartileIntensity)) +
  geom_point() + scale_color_gradient(low = "yellow",high = "blue") + 
  xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
  labs(colour = "UpperQuartile") +
  theme(plot.title = element_text(size=18,hjust=0.5)) +
  ggtitle("Hela") + theme_classic()

ggplot() + 
  geom_point(Zt, mapping = aes(x=V1.1,y=V2.1),
             size = 2, alpha = 1, color = "grey") + 
  geom_point(Zt_livecells, mapping = aes(x=V1.1,y=V2.1,
                                       color=Third_Divide),
             size = 2) +
  scale_color_gradient(low = "blue",high = "yellow") +
  labs("Absolute Time") + theme_classic()

#第一个周期长度和第二个周期长度是否相等
cycleLength <- live_cells
cycleLength <- na.omit(cycleLength)
cycleLength$cycleLength1 <- cycleLength$Second_Divide - cycleLength$First_Divide
cycleLength$cycleLength2 <- cycleLength$Third_Divide - cycleLength$Second_Divide
cycleLength$dealta_cycleLength <- cycleLength$cycleLength2 - cycleLength$cycleLength1
ggplot(cycleLength, aes(y=dealta_cycleLength, x = 1)) + 
  geom_jitter() + geom_boxplot(fill = NA) + theme_classic()
ggplot(cycleLength, aes(y=cycleLength1, x = 1)) + 
  geom_jitter() + geom_boxplot(fill = NA) + theme_classic()
ggplot(cycleLength, aes(y=cycleLength2, x = 1)) + 
  geom_jitter() + geom_boxplot(fill = NA) + theme_classic()
hist(cycleLength$dealta_cycleLength, breaks = 20)
mean(cycleLength$dealta_cycleLength)
sd(cycleLength$dealta_cycleLength)

ggplot(cycleLength, aes(x=dealta_cycleLength)) + 
  geom_histogram(aes(y=..count..), bins = 20, fill = NA, color = "black") + 
  geom_density(aes(y=..scaled..*25), bw = 2, kernel="gaussian") + 
  ylab("Frequency") + xlab("Differences in cycle length over two consecutive divisions") +
  scale_y_continuous(expand = c(0,0),limits = c(0,30)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic()

#子细胞的周期长度是否一致
cycleLength$cellcode1 <- paste0(cycleLength$Nuc_ImageNumber,
                                cycleLength$First_Origin,
                                cycleLength$Second_Origin)
cycleLength$cellcode2 <- paste0(cycleLength$Nuc_ImageNumber,
                                cycleLength$First_Origin,
                                cycleLength$Second_Origin,
                                cycleLength$Third_Origin)
cycleLength$childcyclelength <- NA
for (i in 1:length(unique(cycleLength$cellcode1))) {
  cycleLength[cycleLength$cellcode1==unique(cycleLength$cellcode1)[i],"childcyclelength"] <- abs(unique(cycleLength[cycleLength$cellcode1==unique(cycleLength$cellcode1)[i],"dealta_cycleLength"])[1]-unique(cycleLength[cycleLength$cellcode1==unique(cycleLength$cellcode1)[i],"dealta_cycleLength"])[2])
    }
ggplot(cycleLength, aes(x=childcyclelength)) + 
  geom_histogram(bins = 20, fill = NA, color = "black") + 
  ylab("Difference in daughter cells' cycle length") +
  scale_y_continuous(expand = c(0,0),limits = c(0,30)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic()

#图像信息与细胞周期之间的相关性
live_cells_cor2meta <- live_cells
live_cells_cor2meta <- live_cells_cor2meta[!(live_cells_cor2meta$Nuc_ImageNumber==1&live_cells_cor2meta$Nuc_ObjectNumber==1302),]
live_cells_cor2meta <- live_cells_cor2meta[!(live_cells_cor2meta$Nuc_ImageNumber==2&live_cells_cor2meta$Nuc_ObjectNumber==367),]
live_cells_cor2meta <- live_cells_cor2meta[!(live_cells_cor2meta$Nuc_ImageNumber==2&live_cells_cor2meta$Nuc_ObjectNumber==358),]
for (i in 11:24) {
  live_cells_cor2meta[,i] <- (live_cells_cor2meta[,i]-min(live_cells_cor2meta[,i]))/(max(live_cells_cor2meta[,i])-min(live_cells_cor2meta[,i]))
}
live_cells_cor2meta[,2] <- (live_cells_cor2meta[,2]-min(live_cells_cor2meta[,2]))/(max(live_cells_cor2meta[,2])-min(live_cells_cor2meta[,2]))
live_cells_cor2meta[,5] <- (live_cells_cor2meta[,5]-min(live_cells_cor2meta[,5]))/(max(live_cells_cor2meta[,5])-min(live_cells_cor2meta[,5]))
colnames(live_cells_cor2meta)
setwd("D:/360MoveData/Users/Haofei Ge/Desktop/该读书了/细胞周期伪时间课题（毕业论文版）/活细胞与伪时间的验证数据/各参数与细胞周期之间的相关性")
write.csv(live_cells_cor2meta, file = "live_cells_cor2meta.csv")
ggplot(live_cells_cor2meta, aes(x=abs(V1), y=Nuc_Intensity_IntegratedIntensity)) +
  geom_point(color="#5E6C82") + geom_smooth(color="#5E6C82") + theme_classic()

ggplot(live_cells_cor2meta, aes(x=AbsoluteTime, y=Nuc_Intensity_IntegratedIntensity)) +
  geom_point(color="#5E6C82") + geom_smooth(color="#5E6C82") + theme_classic()
ggsave("Nuc_Intensity_IntegratedIntensity.pdf",device = cairo_pdf,width =10, height =4)
ggplot(live_cells_cor2meta, aes(x=AbsoluteTime, y=Nuc_AreaShape_Area)) +
  geom_point(color="#899FB0") + geom_smooth(color="#899FB0") + theme_classic()
ggsave("Nuc_AreaShape_Area.pdf",device = cairo_pdf,width =10, height =4)
ggplot(live_cells_cor2meta, aes(x=AbsoluteTime, y=Nuc_AreaShape_Compactness)) +
  geom_point(color="#81B3A9") + geom_smooth(color="#81B3A9") + theme_classic()
ggsave("Nuc_AreaShape_Compactness.pdf",device = cairo_pdf,width =10, height =4)
ggplot(live_cells_cor2meta, aes(x=AbsoluteTime, y=Nuc_Intensity_UpperQuartileIntensity)) +
  geom_point(color="#8BB5D1") + geom_smooth(color="#8BB5D1") + theme_classic()
ggsave("Nuc_Intensity_UpperQuartileIntensity.pdf",device = cairo_pdf,width =10, height =4)
ggplot(live_cells_cor2meta, aes(x=AbsoluteTime, y=Nuc_Intensity_MeanIntensity)) +
  geom_point(color="#D6CDBE") + geom_smooth(color="#D6CDBE") + theme_classic()
ggsave("Nuc_Intensity_MeanIntensity.pdf",device = cairo_pdf,width =10, height =4)
ggplot(live_cells_cor2meta, aes(x=AbsoluteTime, y=Nuc_Intensity_StdIntensity)) +
  geom_point(color="#F8D793") + geom_smooth(color="#F8D793") + theme_classic()
ggsave("Nuc_Intensity_StdIntensity.pdf",device = cairo_pdf,width =10, height =4)

ggplot(live_cells_cor2meta, aes(x=AbsoluteTime, y=colnames(live_cells_cor2meta)[5])) +
  geom_point(color="#5E6C82") + geom_smooth(color="#5E6C82") + theme_classic()

tmp <- live_cells_cor2meta[,c("V1",                             "Nuc_Intensity_IntegratedIntensity",   
                              "Nuc_AreaShape_Area",             "Nuc_AreaShape_BoundingBoxArea",       
                              "Nuc_AreaShape_Compactness",      "Nuc_AreaShape_Eccentricity",          
                              "Nuc_AreaShape_FormFactor",       "Nuc_AreaShape_MajorAxisLength",       
                              "Nuc_AreaShape_MaxFeretDiameter", "Nuc_AreaShape_MinFeretDiameter",      
                              "Nuc_AreaShape_MinorAxisLength",  "Nuc_AreaShape_Perimeter",             
                              "Nuc_Intensity_MADIntensity",     "Nuc_Intensity_MaxIntensity",          
                              "Nuc_Intensity_MeanIntensity",    "Nuc_Intensity_MedianIntensity",       
                              "Nuc_Intensity_StdIntensity",     "Nuc_Intensity_UpperQuartileIntensity",
                              "AbsoluteTime")]
for (i in 1:ncol(tmp)) {
  p <- ggplot(tmp, aes(x = AbsoluteTime, y = tmp[[i]])) +
    geom_point(color = "#5E6C82") +
    geom_smooth(color = "#5E6C82") +
    theme(text = element_text(size = 20)) +
    theme_classic() +
    ylab(colnames(tmp)[i])
  
  # 保存图片
  ggsave(paste0("Absolute_", colnames(tmp)[i], ".pdf"), p, width = 4, height = 3,
         path = "D:/360MoveData/Users/Haofei Ge/Desktop/该读书了/细胞周期伪时间课题（毕业论文版）/作图/补充图部分")
}


#各参数之间的相关性
live_cells_cor2meta2 <- live_cells_cor2meta 
library(GGally)
library(rstatix)
library(ggsci)
cols = pal_jco()(7) # 设置绘图颜色
my_density <- function(data, mapping, values = cols ) {
  ggplot(data = data, mapping = mapping) +
    geom_density() +
    scale_color_manual(values = values)
} # 自定义函数
live_cells_cor2meta2$Nuc_ImageNumber <- as.character(live_cells_cor2meta2$Nuc_ImageNumber)
live_cells_cor2meta2$Sample <- "Asyn"
library(stringr)
colnames(live_cells_cor2meta2) <- paste0(str_split(colnames(live_cells_cor2meta2), "Nuc_",simplify = T)[,1],str_split(colnames(live_cells_cor2meta2), "Nuc_",simplify = T)[,2])
colnames(live_cells_cor2meta2) <- paste0(str_split(colnames(live_cells_cor2meta2), "AreaShape_",simplify = T)[,1],str_split(colnames(live_cells_cor2meta2), "AreaShape_",simplify = T)[,2])
colnames(live_cells_cor2meta2) <- paste0(str_split(colnames(live_cells_cor2meta2), "Intensity_",simplify = T)[,1],str_split(colnames(live_cells_cor2meta2), "Intensity_",simplify = T)[,2])
live_cells_cor2meta3 <- live_cells_cor2meta2[,c("Area","Compactness","Eccentricity","FormFactor","MajorAxisLength","MaxFeretDiameter",
                                                "Perimeter","IntegratedIntensity","MADIntensity","MaxIntensity",
                                                "MeanIntensity","MedianIntensity","StdIntensity","UpperQuartileIntensity")]
live_cells_cor2meta3 <- live_cells_cor2meta2[,c("Area","Compactness",
                                                "IntegratedIntensity","UpperQuartileIntensity",
                                                "MeanIntensity","StdIntensity")]
live_cells_cor2meta3$Sample <- "Asyn"


ggpairs(live_cells_cor2meta3,columns = 1:6, 
        ggplot2::aes(color = Sample), #按分组信息分别着色、计算相关性系数和拟合
        upper = list(continuous = "cor"), 
        lower = list(continuous = wrap("smooth",method = "loess",se = FALSE)), #se = FALSE,不展示置信区间 
        diag = list(continuous = my_density)) +
  ggplot2::scale_color_manual(values = "#7AA6DCFF") +
  theme_bw()
ggsave("六参数之间相关性.pdf", device = cairo_pdf, width = 20, height = 20)

tmp1 <- tmp
library(stringr)
colnames(tmp1)[2:18] <- str_split(colnames(tmp1)[2:18],"_",simplify = T)[,3]
colnames(tmp1)[1] <- "V1"
colnames(tmp1)[2] <- "Integrated"
colnames(tmp1)[8] <- "MajorAxis"
colnames(tmp1)[9] <- "MaxFeret"
colnames(tmp1)[10] <- "MinFeret"
colnames(tmp1)[11] <- "MinorAxis"
colnames(tmp1)[13] <- "MAD"
colnames(tmp1)[14] <- "MAX"
colnames(tmp1)[15] <- "Mean"
colnames(tmp1)[16] <- "Median"
colnames(tmp1)[17] <- "Std"
colnames(tmp1)[18] <- "UpperQuartile"

tmp1 <- tmp1[,c("V1", "Area", "BoundingBoxArea",       
                "Compactness", "Eccentricity",  "FormFactor",  "MajorAxis",   
                "MinorAxis", "MaxFeret", "MinFeret", "Perimeter", 
                "Integrated", "MAD", "MAX", "Mean", "Median",       
                "Std", "UpperQuartile", "AbsoluteTime")]
ggpairs(tmp1,columns = 1:18, 
        upper = list(continuous = "cor"), 
        lower = list(continuous = wrap("smooth",method = "loess",se = FALSE)), #se = FALSE,不展示置信区间 
        diag = list(continuous = my_density)) +
  theme(text = element_text(size = 8)) 
ggsave("各参数之间相关性.pdf", device = cairo_pdf, width = 11, height = 8,
       path = "D:/360MoveData/Users/Haofei Ge/Desktop/该读书了/细胞周期伪时间课题（毕业论文版）/作图/补充图部分")
#"#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" "#003C67FF" "#8F7700FF"

#细胞位置密度图
ggplot() + 
  geom_point(Nuc[Nuc$Nuc_ImageNumber==1,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="grey") + 
  geom_point(live_cells_cor2meta[live_cells_cor2meta$Nuc_ImageNumber==1,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="#5E6C82") + 
  geom_density2d(Nuc[Nuc$Nuc_ImageNumber==1,], 
                 mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),
                 color="black", alpha = 0.5, h=110, bins=8, n =150) +
  theme_classic() 
ggplot() + 
  geom_point(Nuc[Nuc$Nuc_ImageNumber==1,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="grey") + 
  geom_point(live_cells_cor2meta[live_cells_cor2meta$Nuc_ImageNumber==1,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="#5E6C82") + 
  geom_density2d(Nuc[Nuc$Nuc_ImageNumber==1,], 
                 mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),
                 color="black", alpha = 0.5, h=110, bins=8, n =150) +
  theme_classic() + 
  xlim(150,950) + ylim(200,1100)

ggplot() + 
  geom_point(Nuc[Nuc$Nuc_ImageNumber==2,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="grey") + 
  geom_point(live_cells_cor2meta[live_cells_cor2meta$Nuc_ImageNumber==2,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="#5E6C82") + 
  geom_density2d(Nuc[Nuc$Nuc_ImageNumber==2,], 
                 mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),
                 color="black", alpha = 0.5, h=120, bins=8, n =150) +
  theme_classic() 
ggplot() + 
  geom_point(Nuc[Nuc$Nuc_ImageNumber==2,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="grey") + 
  geom_point(live_cells_cor2meta[live_cells_cor2meta$Nuc_ImageNumber==2,], 
             mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),color="#5E6C82") + 
  geom_density2d(Nuc[Nuc$Nuc_ImageNumber==2,], 
                 mapping=aes(x=Nuc_Location_Center_X, y=Nuc_Location_Center_Y),
                 color="black", alpha = 0.5, h=120, bins=8, n =150) +
  theme_classic() + 
  xlim(150,1300) + ylim(80,1100)













































