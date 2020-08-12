# Prepare roary input script

library(dplyr)

# Load in classification just as in taxonomyVis.R
classification = read.delim("GTDBtk.stats",stringsAsFactors = FALSE,header = TRUE)
classification$bin = classification$user_genome

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
reassembled$bin=gsub("\\.bin","_bin",reassembled$bin)

# Join dataframes by bin and filter out low completeness bins
bins = left_join(reassembled,classification,by="bin") %>% filter(completeness >= 90) %>% filter(contamination <= 5)
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