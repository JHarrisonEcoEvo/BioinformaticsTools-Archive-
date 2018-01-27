#!/usr/bin/Rscript
#J. G. Harrison
#Jan. 26, 2018

#this script takes an otu table a taxonomy file and a sample to well key and 
#outputs an OTU table that does not have non-target taxa (e.g. no plants if you are doing fungi)

#This script takes several inputs:
#1. a taxonomy file as output from mergeTaxonomyFiles.sh
#2. An otu table
#3. key linking sample to well

#NOTE: it does not remove cpDNA or mtDNA explicitly. 

#IMPORTANT: This script has very limited utility and very likely won't work for you unless you
#tweak it a bit, or are using my pipeline, so carefully read through it!

#For examples of input data see the dir Example_data/ in the git repo with files named as below

#Usage: Rscript useTax_to_cleanOTUtable.R taxonomyFileExample.txt OTUtableExample.txt Well_to_sampleKeyExample.csv

main <- function() {
  inargs <- commandArgs(trailingOnly = TRUE) 
  print(inargs)

  #input otu table, sample-well key, and taxonomy file
  tax=read.delim(file=inargs[1], header=F)
  otus=read.delim(file=inargs[2], header=T)
  sampNames=read.csv(file=inargs[3])
  sampNames$pl_well=paste(sampNames$plate, sampNames$well, sep=".")
  
  #now we need to pick which OTUs to keep
  #I am first removing any OTUs that the UNITE database said were plants.
  #then I will subset the data to just keep only those OTUs that one or the other database
  #had >80% confidence in placement to phylum. In previous work I have found that <80% placement to
  #phylum often match nontarget taxa on NCBI...so I remove those too.

  fungiOnly= tax[tax[,4]!="d:Plantae",]

  #this keeps any row that for either database there is more than 7 characters in the cell
  #describing the taxonomic hypothesis.
  #this gets rid of anything that was not identified to phylum by one or the other db

  fungiOnly = fungiOnly[which(nchar(as.character(fungiOnly[,3])) > 10 | nchar(as.character(fungiOnly[,4])) > 10),]

  ######################################################################################################
  #NOTE: I strongly recommend doing some spot checks of the OTUs that made it through. Find some that were
  #not id'd as fungi by both databases and go to the NCBI web interface and paste them in there.
  #if they show up as some other eukaryote then it may be worth scrubbing the data more
  ######################################################################################################

  #overwrite our OTU table so that it includes only the taxon of interest
  otus=otus[otus$X.OTU.ID %in% fungiOnly$V1,]

  #remove empty columns (no taxa of interest in a sample)
  otus=otus[,which(colSums(otus[,2:length(otus)])!=0)]

  otus=t(otus)
  colnames(otus)=as.character(otus[1,])
  otus=otus[-1,]
  row.names(otus)=gsub("X(\\d+\\.[a-h]\\d+).*","\\1",row.names(otus))

  #check to see if anything in the row.names is not in the key
  which(!(row.names(otus) %in% sampNames$pl_well))

  otus=data.frame(otus)

  #replace well based row names with the sample names
  otus$samps=row.names(otus)
  newotus=merge(otus, sampNames, by.x="samps", by.y="pl_well")

  #make sure to edit the range in the call to the newotus object
  write.csv(data.frame(newotus$sampleName, newotus[,grep("Otu",names(newotus))]),file="cleanOTUtable_youshouldrename.csv", row.names = F)
  
}

main()  
