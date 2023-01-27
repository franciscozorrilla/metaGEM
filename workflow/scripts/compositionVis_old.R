library(gridExtra)
library(dplyr)
library(ggplot2)

classification = read.delim("classification.stats",stringsAsFactors = FALSE,header = FALSE)
colnames(classification)=c("fasta","NCBI","taxonomy","motu","detect","map","percent","cog") # add descriptive column names based on classify-genomes output
classification$percent=as.numeric(classification$percent) # force percentage to be numeric, will default to chr if any bin in dataset cannot be assigned a taxonomy
classification$taxonomy=substr(classification$taxonomy,1,50) # subset taxonomy name, some can get very long, making plot labels obscene
classification$percent[is.na(classification$percent)]  <- 0 # replace NAs with zero's (occurs when bin taxonomy cannot be assigned due to no marker genes)
classification$fasta=gsub(" $","",classification$fasta) # make sure there are not trailing white spaces
classification$sample=gsub("\\..*$","",classification$fasta) # add sample info from bin name

abundance = read.delim("abundance.stats",stringsAsFactors = FALSE,header=FALSE)
colnames(abundance) = c("fasta","ab","abNorm")

classification = left_join(classification,abundance,by="fasta")

plot = ggplot(classification, aes(x = sample, y = abNorm, fill = taxonomy)) + 
  geom_bar(stat = "identity") + 
  ggtitle("Taxonomic composition of samples based on MAGs") +
  xlab("Sample ID") + 
  ylab("Normalized relative abundance") +
  coord_flip()

ggsave("abundanceVis.pdf",plot=plot, height = 8, width = 12)