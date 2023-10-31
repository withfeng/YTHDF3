library(ggstatsplot)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
```



rm(list = ls())


tcga_gsva <- read.csv("CIBERSORT.csv",row.names = 1)
tcga_gsva[1:3,1:3]

#rownames(tcga_gsva) <- gsub("\\.","-",rownames(tcga_gsva))
#tcga_gsva[1:3,1:3]

#基因表达矩阵
tcga_expr <- read.table("YTH.txt", row.names = 1, header = T, check.names = F)
tcga_expr[,1:3]


tcga_gsva <- tcga_gsva[colnames(tcga_expr),]

index <- rownames(tcga_expr) #基因名
y <- as.numeric(tcga_expr)
head(y)
```


colnames <- colnames(tcga_gsva)
data <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(tcga_gsva[,i]),y, method="spearman")
  data[i,2] <- test$estimate                                            
  data[i,3] <- test$p.value
}
names(data) <- c("symbol","correlation","pvalue")
head(data)


write.table(data, "output_cor.txt", sep = "\t", quote = F, row.names = F)
```





data <- read.table("output_cor.txt", sep = "\t", header = T)
head(data)

data<- na.omit(data)
data %>% 
  #filter(pvalue <0.05) %>% # 如果不想把p值大于0.05的放在图上，去掉最前面的#号
  ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
  geom_segment(aes(xend=0,yend=symbol)) +
  geom_point(aes(col=pvalue,size=abs(correlation))) +
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
  #scale_color_viridis_c(begin = 0.5, end = 1) +
  scale_size_continuous(range =c(2,8))  +
  theme_minimal() +
  ylab(NULL)

ggsave("gene_Xcell.pdf")
```