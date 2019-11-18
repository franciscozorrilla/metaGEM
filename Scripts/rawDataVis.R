library(ggplot2)

data = read.delim("dataset.stats",stringsAsFactors = FALSE,header = FALSE,sep = " ")
colnames(data) = c("sample","reads","length")

readsplot= ggplot() + 
  geom_density(data=data,aes(reads)) + 
  geom_vline(aes(xintercept=mean(data$reads)), linetype ="dashed",color = "red") + 
  geom_vline(aes(xintercept=median(data$reads)), linetype ="dashed",color = "blue") + 
  ggtitle("Total reads across all samples in dataset")

ggsave("dataset_reads.pdf",plot= readsplot,device = "pdf",dpi = 300)
