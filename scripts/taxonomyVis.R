library(gridExtra)
library(dplyr)
library(ggplot2)

classification = read.delim("classification.stats",stringsAsFactors = FALSE,header = FALSE)
colnames(classification)=c("fasta","NCBI","taxonomy","motu","detect","map","percent","cog")
classification$percent=as.numeric(classification$percent)
classification$taxonomy=substr(classification$taxonomy,1,40)
classification$percent[is.na(classification$percent)]  <- 0
classification$fasta=gsub(" $","",classification$fasta)

taxplot = classification %>% 
  count(taxonomy) %>% filter(n>10) %>% 
  ggplot(aes(x=reorder(taxonomy,-n),y=n)) +
  geom_bar(stat = "identity") + 
  ggtitle("Taxonomy of reconstructed MAGs") +
  xlab("Taxonomy") +
  ylab("Count") +
  coord_flip()

mapplot=ggplot(classification)+
  geom_density(aes(map),fill="#7fc97f") + 
  ggtitle("Density of marker genes mapped ") +
  xlab("Number of marker genes mapped to MAG") +
  ylab("Density")

detplot=ggplot(classification)+
  geom_density(aes(detect),fill="#beaed4") + 
  ggtitle("Density of marker genes detected") +
  xlab("Number of marker genes detected in MAG") +
  ylab("Density")

perplot=ggplot(classification)+
  geom_density(aes(percent),fill="#fdc086") + 
  ggtitle("Density of agreeing percentage of marker genes") +
  xlab("Percentage of mapped marker genes agreeing with assigned taxomy") +
  ylab("Density")

plotax=grid.arrange(taxplot,arrangeGrob(detplot,mapplot,perplot,nrow=3,ncol = 1),nrow=1,ncol=2)

ggsave("taxonomyVis.pdf",plot= plotax,device = "pdf",dpi = 300, width = 40, height = 20, units = "cm")
