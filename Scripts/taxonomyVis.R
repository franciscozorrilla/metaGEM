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
  count(taxonomy) %>% filter(n>15) %>% 
  ggplot(aes(x=reorder(taxonomy,-n),y=n)) +
  geom_bar(stat = "identity") + 
  coord_flip()

mapplot=ggplot(classification)+
  geom_density(aes(map),fill="#7fc97f")

detplot=ggplot(classification)+
  geom_density(aes(detect),fill="#beaed4")

perplot=ggplot(classification)+
  geom_density(aes(percent),fill="#fdc086")

plotax=grid.arrange(taxplot,arrangeGrob(detplot,mapplot,perplot,nrow=3,ncol = 1),nrow=1,ncol=2)

ggsave("taxonomyVis.pdf",plot= plotax,device = "pdf",dpi = 300, width = 40, height = 20, units = "cm")
