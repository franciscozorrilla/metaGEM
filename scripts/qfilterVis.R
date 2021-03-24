library(gridExtra)
library(dplyr)
library(ggplot2)

qfilter = read.delim("qfilter.stats",stringsAsFactors = FALSE, header = FALSE, sep = " ")
colnames(qfilter) = c("ID","readsBF","readsAF","basesBF","basesAF","percentReads","q20BF","q20AF","q30BF","q30AF")

reads = ggplot(data = qfilter) + 
  geom_density(aes(readsBF,fill="Pre-filtering"),alpha=0.8) + 
  geom_density(aes(readsAF,fill="Post-filtering"),alpha=0.8) +
  ggtitle("Number of reads across samples") + 
  xlab("Total reads") + 
  theme(legend.title = element_blank())

bases = ggplot(data = qfilter) + 
  geom_density(aes(basesBF,fill="Pre-filtering"),alpha=0.8) + 
  geom_density(aes(basesAF,fill="Post-filtering"),alpha=0.8) +
  ggtitle("Number of bases across samples") + 
  xlab("Total bases") + 
  theme(legend.title = element_blank())

q20 = ggplot(data = qfilter) + 
  geom_density(aes(q20BF*100,fill="Pre-filtering"),alpha=0.8) + 
  geom_density(aes(q20AF*100,fill="Post-filtering"),alpha=0.8) +
  ggtitle("Percent Q20 bases across samples") + 
  xlab("Percent of Q20 bases") + 
  theme(legend.title = element_blank())

q30 = ggplot(data = qfilter) + 
  geom_density(aes(q30BF*100,fill="Pre-filtering"),alpha=0.8) + 
  geom_density(aes(q30AF*100,fill="Post-filtering"),alpha=0.8) +
  ggtitle("Percent Q30 bases across samples") + 
  xlab("Percent of Q30 bases") + 
  theme(legend.title = element_blank())

bar = ggplot(data = qfilter) +
  geom_bar(aes(x=reorder(ID,-basesBF),y=basesBF,fill="Pre-filtering"),color="black",stat = "identity") + 
  geom_bar(aes(x=reorder(ID,-basesBF),y=basesAF,fill="Post-filtering"),color="black",stat = "identity") + 
  geom_bar(aes(x=reorder(ID,-basesBF),y=q20AF*basesAF,fill="Q20 bases"),color="black",stat = "identity") + 
  geom_bar(aes(x=reorder(ID,-basesBF),y=q30AF*basesAF,fill="Q30 bases"),color="black",stat = "identity") + 
  coord_flip() +
  ggtitle("Raw read QC summary stacked bar plot") + 
  xlab("Sample ID") + 
  ylab("Base pairs") +
  theme(legend.title = element_blank())

qfilt=grid.arrange(bar,arrangeGrob(reads,bases,q20,q30,nrow=4,ncol=1),ncol =2,nrow=1)
ggsave("qfilterVis.pdf",plot= qfilt,device = "pdf",height = 6, width=8)
