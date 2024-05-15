rm(list = ls())
setwd("D:/360MoveData/Users/Haofei Ge/Desktop/该读书了/细胞周期伪时间课题（毕业论文版）/smFISH in hESC/Pseudo/伪时间和smFISH匹配/")

#smFISH数据准备
Transcription_stats <- read.csv("Transcripts_Stats.csv")
metadata <- readxl::read_xlsx("伪时间和smFISH细胞匹配.xlsx")

Transcription_stats <- merge(Transcription_stats, metadata, by="ObjectNumber")
Transcription_stats <- na.omit(Transcription_stats)

#伪时间数据准备
#读取数据
Nuc <- read.csv("../Result_DNA_Cellpose.csv")
#################################整理数据格式#####################################
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
colnames(Nuc) <- paste("Nuc", colnames(Nuc), sep = "_")#加上前缀，以便于区分Dapi与IF

#初步过滤非正常对象
Nuc <- Nuc[Nuc[,colnames(Nuc)[grep("AreaShape_Area",colnames(Nuc))]]>=5,]

data <- Nuc
#可视化批次效应
library(ggplot2)
library(ggpubr)
ggplot() + 
  geom_density(data, mapping = aes(x=Nuc_Intensity_IntegratedIntensity, 
                                   color=as.factor(Nuc_ImageNumber))) + 
  theme_classic() + ggtitle("Nuc") + labs(color="ImageNumber") + 
  theme(plot.title = element_text(hjust = 0.5))

#校正批次效应方法2，适用于hESC等拥有2N和4N双峰且双峰哪一个更高是不确定的数据类型
#在该方法中，识别双峰的第一峰，将第一峰对齐到参考点
{
  #计算校正因子
  Prop_calculate <- function (data, m = 3){
    #输出表格格式
    prop_result <- as.data.frame(matrix(data=NA,
                                        nrow=length(unique(data$Nuc_ImageNumber)),
                                        ncol=2))
    colnames(prop_result) <- c("Nuc_ImageNumber","Proportion")
    prop_result$Nuc_ImageNumber <- paste("Image",1:nrow(prop_result),sep = "")
    for (i in 1:length(unique(data$Nuc_ImageNumber))) {
      tmp <- density(data[data$Nuc_ImageNumber==i,"Nuc_Intensity_IntegratedIntensity"])
      t <- data.frame(x=tmp$x,y=tmp$y)
      shape <- diff(sign(diff(t$y, na.pad = FALSE)))
      pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(t$y), w, length(t$y))
        if(all(t$y[c(z : i, (i + 2) : w)] <= t$y[i + 1])) return(i + 1) else return(numeric(0))
      })
      pks <- unlist(pks)
      pks_select <- t[pks,]
      dpks_select <- pks_select[order(pks_select$y, decreasing = T),][1:2,]
      dpks_select <- dpks_select[order(dpks_select$x, decreasing = F),"x"]
      #看看
      #plot(t,type="l")
      #abline(v=dpks_select)
      #输出校正因子
      dpks_firstpk <- dpks_select[1]
      prop = 2/dpks_firstpk
      prop_result[i,"Proportion"] <- prop
    }
    return(prop_result)
  }#该函数识别双峰，并计算第一峰的校正因子
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
  Prop_result <- Prop_calculate(data = data_correct, m=3)
  data_correct <- Prop_Normalization(data = data_correct, Prop_Result = Prop_result)
  #检查校正效果
  library(ggpubr)
  p1 <- ggplot() + 
    geom_density(data_correct, mapping = aes(x=Nuc_Intensity_IntegratedIntensity, 
                                             color=as.factor(Nuc_ImageNumber))) + 
    theme_classic() + ggtitle("Nuc") + labs(color="ImageNumber") + 
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggplot() + 
    geom_point(data_correct, mapping = aes(x=Nuc_Intensity_MeanIntensity,
                                           y=Nuc_Intensity_IntegratedIntensity,
                                           color=as.factor(Nuc_ImageNumber))) + 
    theme_classic() + ggtitle("Nuc") + labs(color="ImageNumber") + 
    theme(plot.title = element_text(hjust = 0.5))
  ggarrange(p1,p2)
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
    #cal_ncenter <- function(ncells, ncells_limit = 100){
    #  round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
    #}
    #ncenter <- cal_ncenter(ncol(FM))
    ddrtree_res <- DDRTree(FM,dimensions = 1 , ncenter = 2)
    Z <- ddrtree_res[["Z"]]
    Zt <- as.data.frame(t(Z))
    Zt <- cbind(Zt,Exp)
    ggplot(Zt, aes(x=V1,y=Nuc_Intensity_IntegratedIntensity,
                         color=Nuc_Intensity_UpperQuartileIntensity)) +
      geom_point() + scale_color_gradient(low = "yellow",high = "blue") + 
      xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
      theme(plot.title = element_text(size=18,hjust=0.5)) + ggtitle("Hela ") + 
      theme_classic()
    
    #过滤细胞(手动选择)
    library(igraph)
    library(iplots)
    cells_select <- Zt
    with(cells_select,iplot(V1,Nuc_Intensity_IntegratedIntensity))
    cells_select <- read.table("../cells_select.txt", header = T)
    cells_select <- merge(cells_select,Zt, by=colnames(cells_select))
    ggplot(cells_select, aes(x=V1,y=Nuc_Intensity_IntegratedIntensity,
                                   color=Nuc_Intensity_UpperQuartileIntensity)) +
      geom_point(alpha=0.5) + scale_color_gradient(low = "yellow",high = "blue") + 
      xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
      theme(plot.title = element_text(size=18,hjust=0.5)) + ggtitle("Hela ") + theme_classic()
}

#smFISH展示在伪时间上
cells_select <- merge(cells_select, Transcription_stats, 
                      by = c("Nuc_ImageNumber","Nuc_ObjectNumber"), 
                      all = T)
cells_select_smFISH <- na.omit(cells_select)
ggplot() +
  geom_point(cells_select, mapping=aes(x=V1,y=Nuc_Intensity_IntegratedIntensity),
             alpha=0.7,color="grey",size=1.5) + 
  geom_point(cells_select_smFISH, mapping=aes(x=V1,y=Nuc_Intensity_IntegratedIntensity,
                                              color=OCT4),alpha=1,size=1.5) + 
  scale_color_gradient(low = "yellow",high = "blue") + 
  xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
  theme_bw()
ggplot() +
  geom_point(cells_select, mapping=aes(x=V1,y=Nuc_Intensity_IntegratedIntensity),
             alpha=0.7,color="grey",size=1.5) + 
  geom_point(cells_select_smFISH, mapping=aes(x=V1,y=Nuc_Intensity_IntegratedIntensity,
                                              color=SOX2),alpha=1,size=1.5) + 
  scale_color_gradient(low = "yellow",high = "blue") + 
  xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
  theme_bw()
ggplot() +
  geom_point(cells_select, mapping=aes(x=V1,y=Nuc_Intensity_IntegratedIntensity),
             alpha=0.7,color="grey",size=1.5) + 
  geom_point(cells_select_smFISH, mapping=aes(x=V1,y=Nuc_Intensity_IntegratedIntensity,
                                              color=OCT4/SOX2),alpha=1) + 
  scale_color_gradient(low = "yellow",high = "blue") + 
  xlab("DDRTree_axis") + ylab("DNA IntegratedIntensity") + 
  theme_bw()

#转录本水平随细胞周期的变化定量分析
#先给细胞排序
used2sort <- cells_select_smFISH
write.csv(used2sort, file = "used2sort.csv")

used2sort <- read.csv("used2sort.csv")
colnames(used2sort)[29] <- "OCT4"

ggplot() +
  geom_point(used2sort, mapping=aes(x=sort,y=Nuc_Intensity_IntegratedIntensity),
             alpha=0.7,size=1.5, color="grey") + 
  geom_smooth(used2sort, mapping=aes(x=sort,y=Nuc_Intensity_IntegratedIntensity),
              color="grey") +
  xlab("Pseudo-Time") +
  theme_classic()
ggplot() +
  geom_point(used2sort, mapping=aes(x=sort,y=OCT4),
             alpha=0.7,size=1.5, color="#4A7298") + 
  geom_smooth(used2sort, mapping=aes(x=sort,y=OCT4),
              color="#4A7298") + 
  geom_point(used2sort, mapping=aes(x=sort,y=SOX2),
             alpha=0.7,size=1.5, color="#F3C846") + 
  geom_smooth(used2sort, mapping=aes(x=sort,y=SOX2),
              color="#F3C846") + 
  xlab("Pseudo-Time") + 
  theme_classic()


#新生OCT4和SOX2转录本呈正相关
plot(cells_select_smFISH$OCT4,cells_select_smFISH$SOX2)
cor(cells_select_smFISH$OCT4,cells_select_smFISH$SOX2)





