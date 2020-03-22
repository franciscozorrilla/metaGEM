library(gridExtra)
library(dplyr)
library(ggplot2)

gems = read.delim("GEMs.stats",stringsAsFactors = FALSE,header=FALSE, sep = " ")
gems$V5 = gsub("_.*$","",gems$V1)
colnames(gems) = c("bin","mets","rxns","genes","sample")

samplesplot = gems %>% 
  count(sample) %>% 
  ggplot(aes(x=reorder(sample,-n),y=n)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  ggtitle("Number of GEMs across samples") + 
  ylab("Number of GEMs carved") +
  xlab("Sample ID")

metplot = ggplot() + 
  geom_density(data=gems,aes(mets),fill="#7fc97f") +
  ggtitle("Unique metabolites across GEMs") + 
  theme(legend.position = "none") +
  theme(axis.text.y=element_blank())

rxnplot = ggplot() + 
  geom_density(data=gems,aes(rxns),fill="#beaed4") +
  ggtitle("Reactions across GEMs") + 
  theme(legend.position = "none") +
  theme(axis.text.y=element_blank())

geneplot = ggplot() + 
  geom_density(data=gems,aes(genes),fill="#fdc086") +
  ggtitle("Genes across GEMs") + 
  theme(legend.position = "none") +
  theme(axis.text.y=element_blank())

plot=grid.arrange(samplesplot,arrangeGrob(metplot,rxnplot,geneplot,nrow=3,ncol=1),nrow=1,ncol=2,heights=c(60),widths=c(30,30))

ggsave("modelVis.pdf",plot=plot, height = 8, width = 12)