library(gridExtra)
library(dplyr)
library(ggplot2)

concoctCheckm = read.delim("concoct.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(concoctCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
concoctBins= read.delim("concoct_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(concoctBins) = c("bin","contigs","length")
concoct = left_join(concoctCheckm,concoctBins,by="bin") %>% filter(contamination<=10,completeness>=50)

metabatCheckm = read.delim("metabat.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(metabatCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
metabatBins= read.delim("metabat_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(metabatBins) = c("bin","contigs","length")
metabat = left_join(metabatCheckm,metabatBins,by="bin") %>% filter(contamination<=10,completeness>=50)

maxbinCheckm = read.delim("maxbin.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(maxbinCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
maxbinBins= read.delim("maxbin_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(maxbinBins) = c("bin","contigs","length")
maxbin = left_join(maxbinCheckm,maxbinBins,by="bin") %>% filter(contamination<=10,completeness>=50)
maxbin$contigs = as.numeric(maxbin$contigs)
maxbin$length = as.numeric(maxbin$length)

refinedCheckm = read.delim("refined.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(refinedCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size","set")
refinedBins= read.delim("refined_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(refinedBins) = c("bin","contigs","length")
refined = left_join(refinedCheckm,refinedBins,by="bin") %>% filter(contamination<=10,completeness>=50)

reassembledCheckm = read.delim("reassembled.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(reassembledCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size")
reassembledBins= read.delim("reassembled_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(reassembledBins) = c("bin","contigs","length")
reassembled = left_join(reassembledCheckm,reassembledBins,by="bin") %>% filter(contamination<=10,completeness>=50)

bins <- as.data.frame(matrix(0,nrow = 5,ncol=2))
colnames(bins) = c("variable","value")
bins$variable = c("maxbin2","refined","CONCOCT","metabat2","reassembled")
bins$value = c(as.numeric(dim(maxbin)[1]),as.numeric(dim(refined)[1]),as.numeric(dim(concoct)[1]),as.numeric(dim(metabat)[1]),as.numeric(dim(reassembled)[1]))

binplot = ggplot(data = bins,aes(x=reorder(variable,-value),y=value,fill= variable)) +
  geom_bar(stat = "identity",color="black") +
  ylab("Number of bins generated") + 
  xlab("Binning tool") +
  theme(legend.title = element_blank()) +
  ggtitle("Number of MQ bins generated") + 
  coord_flip() + 
  theme(legend.position = "none")

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

covplot = ggplot() + 
  geom_density(data=concoct,aes(legnth,color="CONCOCT")) +
  geom_density(data=maxbin,aes(legnth,color="maxbin2")) +
  geom_density(data=metabat,aes(legnth,color="metabat2")) + 
  geom_density(data=refined,aes(legnth,color="refined")) + 
  geom_density(data=reassembled,aes(legnth,color="reassembled")) +
  ggtitle("legnth") + 
  theme(legend.position = "none") +
  theme(axis.text.y=element_blank())

gcplot = ggplot() + 
  geom_density(data=concoct,aes(GC,color="CONCOCT")) +
  geom_density(data=maxbin,aes(GC,color="maxbin2")) +
  geom_density(data=metabat,aes(GC,color="metabat2")) + 
  geom_density(data=refined,aes(GC,color="refined")) + 
  geom_density(data=reassembled,aes(GC,color="reassembled")) +
  ggtitle("GC content") + 
  theme(legend.position = "none") +
  theme(axis.text.y=element_blank())

densities1= grid.arrange(compplot,lengthplot,covplot,nrow=3,ncol=1)
densities2=grid.arrange(contplot,contigplot,gcplot,nrow=3,ncol=1)

grid.arrange(binplot,densities1,densities2,nrow=1,ncol=3)

plot=grid.arrange(binplot,arrangeGrob(compplot,lengthplot,covplot,nrow=3,ncol=1),arrangeGrob(contplot,contigplot,gcplot,nrow=3,ncol=1),nrow=1,ncol=3,heights=c(60),widths=c(30,30,30))

ggsave("binningVis.pdf",plot=plot, height = 8, width = 12)