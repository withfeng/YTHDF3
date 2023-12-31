df <- read.table("input.txt",head=T,sep="\t",check.names = F)
head(df)
library("pROC")

#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("#a17db4", "#8EA5C8", "#B3D6AD","#ada579","#D2AEAC","#FF9900","#223D6C","#D20A13")

pdf("ROC2.pdf",height=6,width=6)
auc.out <- c()


#先画第一条线，此处是miRNA1
x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              ci=TRUE, 
              main="",
              #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
              col=mycol[2],#线的颜色
              lwd=2, #线的粗细
              legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1

ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限

auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
auc.out <- rbind(auc.out,auc.ci)


#再用循环画第二条和后面更多条曲线
for (i in 3:ncol(df)){
  x <- plot.roc(df[,1],df[,i],
                add=T, #向前面画的图里添加
                smooth=F,
                ci=TRUE,
                col=mycol[i],
                lwd=2,
                legacy.axes=T)
  
  ci.lower <- round(as.numeric(x$ci[1]),3)
  ci.upper <- round(as.numeric(x$ci[3]),3)
  
  auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
}


#对比多条曲线
#在参数`method=`后面，有三种方法可选“delong”, “bootstrap”或“venkatraman”，计算p值
p.out <- c()
for (i in 2:(ncol(df)-1)){
  for (j in (i+1):ncol(df)){
    p <- roc.test(df[,1],df[,i],df[,j], method="bootstrap")
    p.tmp <- c(colnames(df)[i],colnames(df)[j],p$p.value)
    p.out <- rbind(p.out,p.tmp)
  }
}

#输出p value到文件
p.out <- as.data.frame(p.out)
colnames(p.out) <- c("ROC1","ROC2","p.value")
write.table(p.out,"pvalue_output.xls",sep="\t",quote=F,row.names = F,col.names = T)

#还可以把p value写在图上
#这里有4条线，6组对比。太多，就不写了吧。
#如果只对比两条线，就运行下面这行
#text(0.4, 0.3, labels=paste("miRNA1 vs. miRNA2\np-value =", p.out[1,3]), adj=c(0, .5))


# 输出AUC、AUC CI到文件
auc.out <- as.data.frame(auc.out)
colnames(auc.out) <- c("Name","AUC","AUC CI")
write.table(auc.out,"auc_output.xls",sep="\t",quote = F,row.names = F,col.names = T)


#绘制图例
legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:length(df)],
       lwd = 2,
       bty="n")
dev.off()