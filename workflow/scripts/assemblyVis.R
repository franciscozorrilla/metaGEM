library(ggplot2)
library(gridExtra)

assembly = read.delim("assembly.stats",stringsAsFactors = FALSE,header = FALSE,sep = " ")
colnames(assembly) = c("sample","contigs","length_total")

assembly$ave = assembly$length_total/assembly$contigs

aveplot = ggplot(data=assembly) +
  geom_density(aes(x=ave,fill="All contigs")) +
  xlab("Average contig length") +
  ggtitle("Average contig length across samples") +
  theme(legend.title = element_blank())+ 
  theme(legend.position = "none")

contplot = ggplot(data=assembly) +
  geom_density(aes(x=contigs,fill="Total contigs")) +
  ggtitle("Contigs across samples") +
  theme(legend.title = element_blank()) + 
  scale_x_log10() + 
  theme(legend.position = "none")

scatplot = ggplot(data=assembly) +
  geom_point(aes(x=length_total,y=contigs)) +
  ggtitle("Total length vs number of contigs")+
  xlab("Total length") +
  ylab("Number of contigs") +
  expand_limits(x = 0, y = 0)

barplot = ggplot(data=assembly) +
  geom_bar(aes(x=reorder(sample,-length_total),y=length_total,fill="Contigs"),stat = "identity",color="black",size=0.2) +
  ggtitle("Total length across assemblies") + 
  ylab("Length (bp)") +
  xlab("Sample ID") +
  #scale_y_log10() +
  coord_flip() + 
  theme(legend.position = "none")

assemblyplot=grid.arrange(barplot,arrangeGrob(scatplot,aveplot,contplot,nrow=3),ncol=2,nrow=1)

ggsave("assemblyVis.pdf",plot= assemblyplot,device = "pdf",dpi = 300, width = 30, height = 40, units = "cm")
