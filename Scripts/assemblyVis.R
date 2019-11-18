library(ggplot2)
library(gridExtra)

data = read.delim("dataset.stats",stringsAsFactors = FALSE,header = FALSE,sep = " ")
colnames(data) = c("sample","reads","length")

readsplot= ggplot() + 
  geom_density(data=data,aes(reads)) + 
  geom_vline(aes(xintercept=mean(data$reads)), linetype ="dashed",color = "red") + 
  geom_vline(aes(xintercept=median(data$reads)), linetype ="dashed",color = "blue") + 
  ggtitle("Reads across samples")

assembly = read.delim("assembly.stats",stringsAsFactors = FALSE,header = FALSE,sep = " ")
colnames(assembly) = c("sample","contigs","length","coverage")

assemblyplotC= ggplot() + 
  geom_density(data=assembly,aes(contigs)) + 
  geom_vline(aes(xintercept=mean(assembly$contigs)), linetype ="dashed",color = "red") + 
  geom_vline(aes(xintercept=median(assembly$contigs)), linetype ="dashed",color = "blue") + 
  ggtitle("Contigs across assemblies")

assemblyplotL= ggplot() + 
  geom_density(data=assembly,aes(length)) + 
  geom_vline(aes(xintercept=mean(assembly$length)), linetype ="dashed",color = "red") + 
  geom_vline(aes(xintercept=median(assembly$length)), linetype ="dashed",color = "blue") + 
  ggtitle("Length across assemblies")

assemblyplotCov= ggplot() + 
  geom_density(data=assembly,aes(coverage)) + 
  geom_vline(aes(xintercept=mean(assembly$coverage)), linetype ="dashed",color = "red") + 
  geom_vline(aes(xintercept=median(assembly$coverage)), linetype ="dashed",color = "blue") + 
  ggtitle("Ave. cov. across assemblies")

grid.arrange(readsplot,assemblyplotC,assemblyplotL,assemblyplotCov,ncol=2,nrow=2)