# Prepare roary input script

library(dplyr)
library(tidyr)

# Load in classification just as in taxonomyVis.R
classification = read.delim("GTDBtk.stats",stringsAsFactors = FALSE,header = TRUE)
classification$bin = classification$user_genome
gtdbtk_class = classification[,c(2,20)]
gtdbtk_class$classification = gsub("*.__","",gtdbtk_class$classification)
gtdbtk_class %>% separate(classification,c("domain","phylum","class","order","family","genus","species"),sep = ";") -> gtdbtk_class
gtdbtk_class[gtdbtk_class==""]<-'NA'
gtdbtk_class$sample = gsub("_.*$","",gtdbtk_class$bin)

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
bins = left_join(reassembled,gtdbtk_class,by="bin") %>% filter(completeness >= 90) %>% filter(contamination <= 5)

# Identify which species are represented by at least 10 high quality bins
bins %>% group_by(species) %>% count() %>% filter(n>=10) -> species

# Run for loop to generate text file with bin IDs for each identified species
dir.create(gsub("$","/speciesBinIDs",getwd()))

for (i in species$species) {
  #Remove any forbidden characters if present:
  #(spaces, forward slashes, square brackets, parentheses)
  name = gsub(" ","_",i)
  name = gsub("/","_",name)
  name = gsub("\\[","_",name)
  name = gsub("]","_",name)
  name = gsub("\\(","_",name)
  name = gsub(")","_",name)
  write.table(bins %>% 
                filter(species == i) %>%
                select(bin),paste0(paste0("speciesBinIDs/",name),".txt"),sep="\n",row.names=FALSE,col.names = FALSE,quote = FALSE)
}