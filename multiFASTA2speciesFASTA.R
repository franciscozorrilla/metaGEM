# Author: Francisco Zorrilla
# Contact: zorrilla@chalmers.se
# Chalmers University of Technology, Sweden
# 10/09/2018

# This script is part of a larger metagenomics pipeline developed for my thesis project.
# Please take a look at the Snakefile (particularly the parseFASTA rule) and config file to understand script usage.
# Script will extract .fa and .faa files for each identified cluster (species) from a multifasta file based on CONCOCT output.

#Load libraries
source("https://bioconductor.org/biocLite.R")
library(BiocGenerics, lib.loc = "/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/R/")
library(S4Vectors, lib.loc = "/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/R/")
library(IRanges ,lib="/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/R")
library(XVector, lib.loc = "/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/R/")
library(Biostrings, lib.loc = "/c3se/NOBACKUP/groups/c3-c3se605-17-8/projects_francisco/R/")

#Define x containing multispecies protein FASTA with clean IDs
x = as.data.frame(readAAStringSet("cleanID.faa"))
x = cbind(rownames(x),x)
names(x)= NULL

#Define y containing IDs/nodes for each species/cluster
y = read.table("allspecies.txt") #IDs

#Define z containing multispecies DNA FASTA with clean IDs
z = as.data.frame(readDNAStringSet("cleanID.fa"))
z = cbind(rownames(z),z)
names(x)= NULL

#Extract protein FASTA sequences for each species
for (i in 1:length(rownames(y)))  {          
  df <- as.data.frame(x[grep(paste(y[i,], sep="|",collapse ="T"),rownames(x)),])
  assign(paste("species", i, sep = "") , df , envir = globalenv()) 
  write.table(df,paste("speciesProt", i, ".txt",sep = ""),sep="\t",row.names=FALSE)
}    

df <- 0

#Extract DNA FASTA sequences for each species
for (i in 1:length(rownames(y)))  {          
  df <- as.data.frame(z[grep(paste(y[i,], sep="|",collapse ="T"),rownames(z)),])
  assign(paste("species", i, sep = "") , df , envir = globalenv()) 
  write.table(df,paste("speciesDNA", i, ".txt",sep = ""),sep="\t",row.names=FALSE)
}    
