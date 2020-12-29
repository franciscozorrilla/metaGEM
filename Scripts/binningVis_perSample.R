library(gridExtra)
library(dplyr)
library(ggplot2)

concoctCheckm = read.delim("concoct.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(concoctCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
concoctBins= read.delim("concoct_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(concoctBins) = c("bin","contigs","length")
concoct = left_join(concoctCheckm,concoctBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50) %>% distinct() %>% select(-set)
concoct$sample = gsub("\\..*$","",concoct$bin)
concoct$binner = "CONCOCT"

metabatCheckm = read.delim("metabat.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(metabatCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
metabatBins= read.delim("metabat_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(metabatBins) = c("bin","contigs","length")
metabat = left_join(metabatCheckm,metabatBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct() %>% select(-set)
metabat$sample = gsub("\\..*$","",metabat$bin)
metabat$binner = "MetaBAT2"

maxbinCheckm = read.delim("maxbin.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(maxbinCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
maxbinBins= read.delim("maxbin_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(maxbinBins) = c("bin","contigs","length")
maxbin = left_join(maxbinCheckm,maxbinBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct() %>% select(-set)
maxbin$contigs = as.numeric(maxbin$contigs)
maxbin$sample = gsub("\\..*$","",maxbin$bin)
maxbin$binner = "MaxBin2"

refinedCheckm = read.delim("refined.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(refinedCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
refinedBins= read.delim("refined_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(refinedBins) = c("bin","contigs","length")
refined = left_join(refinedCheckm,refinedBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct() %>% select(-set)
refined$sample = gsub("\\..*$","",refined$bin)
refined$binner = "metaWRAP_refined"

reassembledCheckm = read.delim("reassembled.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(reassembledCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size")
reassembledBins= read.delim("reassembled_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(reassembledBins) = c("bin","contigs","length")
reassembled = left_join(reassembledCheckm,reassembledBins%>%select(-length),by="bin") %>% filter(contamination<=10,completeness>=50)%>% distinct()
reassembled$sample = gsub("\\..*$","",reassembled$bin)
reassembled$binner = "metaWRAP_reassembled"

#bins <- as.data.frame(matrix(0,nrow = 5,ncol=2))
#colnames(bins) = c("variable","value")
#bins$variable = c("maxbin2","refined","CONCOCT","metabat2","reassembled")
#bins$value = c(as.numeric(dim(maxbin)[1]),as.numeric(dim(refined)[1]),as.numeric(dim(concoct)[1]),as.numeric(dim(metabat)[1]),as.numeric(dim(reassembled)[1]))

rbind(concoct,metabat,maxbin,refined,reassembled) %>% group_by(binner,sample) %>% summarize(count=n()) -> bins

binplot = ggplot(data = bins,aes(x=reorder(binner,-count),y=count,fill= binner)) +
  geom_bar(stat = "identity",color="black") +
  ylab("Generated bins") + 
  xlab("Binning tool") +
  theme(legend.title = element_blank()) +
  ggtitle("Number of MQ bins") + 
  coord_flip() + 
  theme(legend.position = "none")+ 
  facet_wrap(~sample,ncol=1)

compplot = ggplot() + 
  geom_density(data=concoct,aes(completeness,color="CONCOCT")) +
  geom_density(data=maxbin,aes(completeness,color="maxbin2")) +
  geom_density(data=metabat,aes(completeness,color="metabat2")) + 
  geom_density(data=refined,aes(completeness,color="refined")) + 
  geom_density(data=reassembled,aes(completeness,color="reassembled")) +
  ggtitle("Completeness") +
  theme(axis.text.y=element_blank()) + 
  theme(legend.position = "none")

contplot = ggplot() + 
  geom_density(data=concoct,aes(contamination,color="CONCOCT")) +
  geom_density(data=maxbin,aes(contamination,color="maxbin2")) +
  geom_density(data=metabat,aes(contamination,color="metabat2")) + 
  geom_density(data=refined,aes(contamination,color="refined")) + 
  geom_density(data=reassembled,aes(contamination,color="reassembled")) +
  ggtitle("Contamination") +
  theme(axis.text.y=element_blank()) + 
  theme(legend.position = "none")

lengthplot = ggplot() + 
  geom_density(data=concoct,aes(size,color="CONCOCT")) +
  geom_density(data=maxbin,aes(size,color="maxbin2")) +
  geom_density(data=metabat,aes(size,color="metabat2")) + 
  geom_density(data=refined,aes(size,color="refined")) + 
  geom_density(data=reassembled,aes(size,color="reassembled")) +
  ggtitle("BP Length") + 
  theme(legend.position = "none") +
  theme(axis.text.y=element_blank())

contigplot = ggplot() + 
  geom_density(data=concoct,aes(contigs,color="CONCOCT")) +
  geom_density(data=maxbin,aes(contigs,color="maxbin2")) +
  geom_density(data=metabat,aes(contigs,color="metabat2")) + 
  geom_density(data=refined,aes(contigs,color="refined")) + 
  geom_density(data=reassembled,aes(contigs,color="reassembled")) +
  ggtitle("Number of contigs") + 
  theme(legend.position = "none") +
  theme(axis.text.y=element_blank())

densities1= grid.arrange(compplot,lengthplot,nrow=2,ncol=1)
densities2=grid.arrange(contplot,contigplot,nrow=2,ncol=1)

plot=grid.arrange(binplot,densities1,densities2,nrow=1,ncol=3)

ggsave("binningVis.pdf",plot=plot, height = 6, width = 12)