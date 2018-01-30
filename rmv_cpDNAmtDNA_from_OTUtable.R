#!/usr/bin/Rscript
#J. G. Harrison
#Jan. 26, 2018

#this script takes an otu table a taxonomy file and a sample to well key and
#outputs an OTU table that does not have non-target taxa (e.g. no plants if you
#are doing fungi). It specifically searches the taxonomy file for the words
#chloroplast and mitochondria

#This script takes several inputs:
#1. a taxonomy file as output from mergeTaxonomyFiles.sh
#2. An otu table
#3. key linking sample to well

#IMPORTANT: This script has very limited utility and very likely won't work for
#you unless you tweak it a bit, or are using my pipeline, so carefully read it!

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
  fungiOnly= tax[-(grep("Chloroplast", tax[,4])),]

  ######################################################################################################
  #NOTE: I strongly recommend doing some spot checks of the OTUs that made it through. Find some that were
  #not id'd as fungi by both databases and go to the NCBI web interface and paste them in there.
  #if they show up as some other eukaryote then it may be worth scrubbing the data more
  ######################################################################################################

  #make new OTU table so that it includes only the taxon of interest
  otus2=otus[otus$X.OTU.ID %in% fungiOnly$V1,]

  #remove empty columns (no taxa of interest in a sample)
  otus2=otus2[,which(colSums(otus2[,2:length(otus2)])!=0)]

  otus2=t(otus2)
  colnames(otus2)=as.character(otus2[1,])
  otus2=otus2[-1,]
  row.names(otus2)=gsub("X(\\d+\\.[a-h]\\d+).*","\\1",row.names(otus2))

  #check to see if anything in the row.names is not in the key
  if(length(which(!(row.names(otus2) %in% sampNames$pl_well)))>1){
    print("Error: some row names of the output otu file are not in the key")
    print("This means that you have some unidentified samples or more likely")
    print("something got mismatched somewhere. Run script piecemeal to see what happened!")
  }

  otus2=data.frame(otus2)

  #replace well based row names with the sample names
  otus2$samps = row.names(otus2)
  newotus = merge(otus2, sampNames, by.x="samps", by.y="pl_well")

  #make sure to edit the range in the call to the newotus object
  outputfile = data.frame(newotus$sampleName, newotus[,grep("Otu",names(newotus))])

  write.csv(outputfile,file="cleanOTUtable_youshouldrename.csv", row.names = F)

}

main()
