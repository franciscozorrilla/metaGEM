library(tidyverse)
library(tidytext)

taxonomy=read.delim("GTDBTk.stats",header=TRUE) %>% 
  select(user_genome,classification) %>% 
  separate(.,classification,into = c("kingdom","phylum","class","order","family","genus","species"),sep = ";")

abundance=read.delim("abundance.stats",header=FALSE)
colnames(abundance)=c("user_genome","absolute_ab","rel_ab")

taxab = left_join(taxonomy,abundance,by="user_genome")
taxab$sample = gsub("\\..*$","",taxab$user_genome)
taxab$species = gsub("s__$","Undefined sp.",taxab$species)
taxab$species = gsub("s__","",taxab$species)

ggplot(taxab%>% filter(species!="Undefined sp.")) +
  geom_bar(aes(x=reorder_within(species,-rel_ab,sample),y=rel_ab*100),stat="identity") + 
  scale_x_reordered() +
  facet_wrap(~sample,scales = "free") + 
  ylab("Relative abundance (%)") + 
  xlab("Species") +
  coord_flip()

ggsave("compositionVis.pdf",width = 12,height=8)