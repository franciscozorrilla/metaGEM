library(ggplot2)
library(gridExtra)

assembly = read.delim("assembly.stats",stringsAsFactors = FALSE,header = FALSE,sep = " ")
colnames(assembly) = c("sample","contigs","length_total","ave_contig","gt1000","gt1000_total")

contplot = ggplot(data=assembly) +
  geom_density(aes(x=contigs,fill="Total contigs")) +
  geom_density(aes(x=gt1000,fill="Contigs >= 1000")) +
  ggtitle("Contigs across samples") +
  theme(legend.title = element_blank()) + 
  scale_x_log10()

contplot2 = ggplot(data=assembly) +
  geom_point(aes(x=contigs,y=gt1000)) +
  xlab("Total contigs") + 
  ylab("Contigs >= 1000bp") +
  ggtitle("Total contigs vs >=1000bp across samples") +
  geom_abline(slope = 1,intercept=0) +
  expand_limits(x = 0, y = 0)

lenplot = ggplot(data=assembly) +
  geom_density(aes(x=length_total,fill="Total length")) +
  geom_density(aes(x=gt1000_total,fill="Length >= 1000")) +
  xlab("Length")
  ggtitle("Length across samples") +
  theme(legend.title = element_blank()) + 
  scale_x_log10()

lenplot2= ggplot(data=assembly) +
  geom_point(aes(x=length_total,y=gt1000_total)) +
  ggtitle("Total length vs >= 1000bp across samples")+
  xlab("Total length") +
  ylab("Length of contigs >= 1000 bp") +
  geom_abline(slope = 1,intercept=0) +
  expand_limits(x = 0, y = 0)

fracplot = ggplot(data=assembly) +
  geom_density(aes(100*gt1000/contigs,fill="Contigs")) +
  geom_density(aes(100*gt1000_total/length_total,fill="Length")) +
  ggtitle("Percent of information represented by contigs >= 1000 bp") + 
  xlab("% Information") +
  theme(legend.title = element_blank())

assemblyplot=grid.arrange(fracplot,arrangeGrob(contplot2,contplot,lenplot2,lenplot,nrow=2,ncol=2),ncol=2,nrow=1)

ggsave("assemblyVis.pdf",plot= assemblyplot,device = "pdf",dpi = 300)