# Prepare roary input script

library(dplyr)

# Load in classification just as in taxonomyVis.R
classification = read.delim("classification.stats",stringsAsFactors = FALSE,header = FALSE)
colnames(classification)=c("bin","NCBI","taxonomy","motu","detect","map","percent","cog")
classification$percent=as.numeric(classification$percent)
classification$taxonomy=substr(classification$taxonomy,1,40)
classification$percent[is.na(classification$percent)]  <- 0
classification$bin=gsub(" $","",classification$bin)

# Load in refined+reassembled consensus bins just as in binningVis.R
reassembledCheckm = read.delim("reassembled.checkm",stringsAsFactors = FALSE,header = FALSE)
colnames(reassembledCheckm) = c("bin","completeness","contamination","GC","lineage","N50","size")
reassembledBins= read.delim("reassembled_bins.stats",stringsAsFactors = FALSE,header = FALSE, sep = " ")
colnames(reassembledBins) = c("bin","contigs","length")
reassembled = left_join(reassembledCheckm,reassembledBins,by="bin")
reassembled$bin=gsub(" $","",reassembled$bin)
reassembled$bin=gsub("permissive","p",reassembled$bin)
reassembled$bin=gsub("orig","o",reassembled$bin)
reassembled$bin=gsub("strict","s",reassembled$bin)

# Join dataframes by bin and filter out low completeness bins
bins = left_join(reassembled,classification,by="bin") %>% filter(completeness >= 90)
bins$taxonomy = gsub("^ ","",bins$taxonomy)
bins$taxonomy = gsub(" $","",bins$taxonomy)

# Identify which species are represented by at least 10 high quality bins
bins %>% group_by(taxonomy) %>% count() %>% filter(n>=10) -> species

# Run for loop to generate text file with bin IDs for each identified species
dir.create(gsub("$","/speciesBinIDs",getwd()))

for (i in species$taxonomy) {
  name = gsub(" ","_",i)
  name = gsub("/","_",name)
  name = gsub("\\[","_",name)
  name = gsub("]","_",name)
  name = gsub("\\(","_",name)
  name = gsub(")","_",name)
  write.table(bins %>% 
                filter(taxonomy == i) %>%
                select(bin),paste0(paste0("speciesBinIDs/",name),".txt"),sep="\n",row.names=FALSE,col.names = FALSE,quote = FALSE)
}