#!/usr/bin/Rscript
#J. G. Harrison
#April 27, 2018

#this script takes an otu table and a taxonomy file and
#outputs an OTU table that does not have non-target taxa (e.g. no plants if you
#are doing fungi). It specifically searches the taxonomy file for the words
#chloroplast, plantae, and mitochondria. It combines plants into one row, mitochondria into another, the ISD
#into another, and stuff that doesn't match any known phyla into another. 

#This script takes several inputs:
#1. a taxonomy file as output from mergeTaxonomyFiles.sh
#2. An otu table
#3. The output of a usearch/vsearch search of OTU sequences to find the ISD in blast6 format
#To make the isd input one can use a command like this:
#vsearch --usearch_global ${f} --db  /project/microbiome/ref_db/synthgene.fasta 
#--strand both --blast6out ${f}_isd --id 0.85 --threads 32

#IMPORTANT: If you do not check that your input files are correctly formatted this script will not work.

#For examples of input data see the dir Example_data/ in the git repo

#Usage: Rscript cleanOTUtable.R taxonomyExample.txt OTUtableExample.txt ISDhitsExample


# Code for development, users should ignore 
# tax <- read.delim(file="processedData/its_90.sintax", header=F,stringsAsFactors = F)
# otus=read.delim(file="processedData/otuTables/its90_otuTable", header=T,stringsAsFactors = F)
# isd=read.delim(file="~/git_repo/bioinformaticsTools/Example_data/ISDhitsExample", header=T,stringsAsFactors = F)

main <- function() {
  inargs <- commandArgs(trailingOnly = TRUE)
  print(inargs)

  #input otu table, sample-well key, and taxonomy file
  tax=read.delim(file=inargs[1], header=F, stringsAsFactors = F)
  otus=read.delim(file=inargs[2], header=T, stringsAsFactors = F)
  isd=read.delim(file=inargs[3], header=F,stringsAsFactors = F)
  
  print(paste("Started with ", length(tax[,1]), " taxa", sep=""))

  plantIndices <- c(grep("[pP]lant", tax[,4]),
    grep("[cC]hloroplast", tax[,4]))

  if( length(plantIndices) > 0 ){
    plantOtus <- tax[plantIndices, 1]
    plantsum <- colSums(otus[otus[,1] %in% plantOtus,2:length(otus)])
    otus <- otus[!(otus[,1] %in% plantOtus),]
    otus[1+length(otus[,1]),] <- c("plant",plantsum)
  }
  
  duds <- which(nchar(tax[,4]) == 0) #taxa not assigned a taxonomic hypothesis
  if( length(duds) > 0 ){
    dudOtus <- tax[duds, 1]
    dudsum <- colSums(otus[otus[,1] %in% dudOtus,2:length(otus)])
    otus <- otus[!(otus[,1] %in% dudOtus),]
    otus[1+length(otus[,1]),] <- c("duds",dudsum)
  }
  
  mtdna <- c(grep("[mM]itochondria", tax[,4]))
  
  if( length(mtdna) > 0 ){
    mtOtu <- tax[mtdna, 1]
    mtsum <- colSums(otus[otus[,1] %in% mtOtu,2:length(otus)])
    otus <- otus[!(otus[,1] %in% mtOtu),]
    otus[1+length(otus[,1]),] <- c("mtDNA",mtsum)
  }
  
  #remove ISD
  isdIndices <- which(otus$OTUID %in% isd$V1)
  if( length(isdIndices) > 0 ){
    isdsum <- colSums(otus[otus$OTUID %in% isd$V1,2:length(otus)])
    otus <- otus[!(otus$OTUID %in% isd$V1),]
    otus[1+length(otus[,1]),] <- c("ISD",isdsum)
  }
  
  ######################################################################################################
  #NOTE: I strongly recommend doing some spot checks of the OTUs that made it through. Find some that were
  #not id'd and go to the NCBI web interface and paste the sequences in there. Could be that some of the
  #stuff removed is a target taxon.
  
  #
  #Also, CHECK that this script worked right. Do your own wrangling to check that the ISD summing worked
  #for example. If you do not understand what this script is doing, then it is probably best to not 
  #use it until you can get clarification. Talk to Josh if you want something here explained. This sort
  # of wrangling is the main way that errors happen (not in analysis), so it is worth being really careful.
  ######################################################################################################
  
  write.csv(otus,file="cleanOTUtable_youshouldrename.csv", row.names = F)
}

main()
