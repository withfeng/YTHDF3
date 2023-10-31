
#install.packages("rms")
#install.packages("rmda")


#???ð?
library(rms)
library(rmda)

inputFile="normalize.txt"             #?????????ļ?
geneFile="importanceGene.XGB.txt"     #?????б??ļ?
#setwd("E:\\???ŷ???\\2023??5?·????ŷ???\\10.DNA?׻??????ݷ???\\3.ģ?͹???")      #???ù???Ŀ¼

#??ȡ?????ļ?
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

#??ȡ?????б??ļ?,??ȡ?????????????ı???��
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#??ȡ??Ʒ??????Ϣ
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

#??
# 查看哪些基因名在data的行名中不存在
not_in_data <- geneRT[,1][!(geneRT[,1] %in% row.names(data))]
print(not_in_data)
?ݴ???
ddist=datadist(rt)
options(datadist="ddist")

#????ģ?ͣ?????????ͼ
lrmModel=lrm(Type~ ZC3HYTHDF3+METTL3+METTL14+ALKBH5+YTHDC2ata=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.1,0.3,0.99),
	lp=F, funlabel="Risk of Disease")
#????????ͼ
pdf("Nomo.pdf", width=8, h17ight=6)
plot(nomo)
dev.off()

#????У׼????
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()

#???ƾ???????
rt$Type=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ ZC3HYTHDF3+METTL3+METTL14+ALKBH5+YTHDC2ta=rt, 
	family = binomial(link ='logit'),
	thresholds= seq(0,1,by = 0.01),
	confidence.intervals = 0.95)
#????DCAͼ??
pdf(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
	curve.names="Model",
	xlab="Threshold probability",
	cost.benefit.axis=T,
	col="red",
	confidence.intervals=FALSE,
	standardize=FALSE)
dev.off()


