#!/usr/bin/Rscript
#J. G. Harrison
#April 27, 2018

#this script takes an otu table and a taxonomy file and
#outputs an OTU table that does not have non-target taxa (e.g. no plants if you
#are doing fungi). It specifically searches the taxonomy file for the words
#chloroplast, plantae, and mitochondria

#This script takes several inputs:
#1. a taxonomy file as output from mergeTaxonomyFiles.sh
#2. An otu table

#IMPORTANT: This script has very limited utility and very likely won't work for
#you unless you tweak it a bit, or are using my pipeline, so carefully read it!

#For examples of input data see the dir Example_data/ in the git repo

#Usage: Rscript cleanOTUtable.R combinedTaxonomyfile.txt OTUtableExample.txt clean

#where clean is either "yes" or "no" and determines if you want
#taxa removed that were not assigned to a phylum

main <- function() {
  inargs <- commandArgs(trailingOnly = TRUE)
  print(inargs)

  #input otu table, sample-well key, and taxonomy file
  tax=read.delim(file=inargs[1], header=F)
  otus=read.delim(file=inargs[2], header=T)
  clean=inargs[3]

  #now we need to pick which OTUs to keep
  #I am first removing any OTUs that either database said were plants.
  #then I will subset the data to just keep only those OTUs that one or the other database
  #had >80% confidence in placement to phylum. In previous work I have found that <80% placement to
  #phylum often match nontarget taxa on NCBI...so I remove those too.

  #note these are just operating on the 4th field because the UNITE database out has a different format
  print(paste("Started with ", length(tax[,1]), " taxa", sep=""))

  fungiOnly <- tax[tax[,3]!="d:Plantae",]
  print(paste("Now have ", length(fungiOnly[,1]), " taxa after removing plant otus found by taxonomy db 1", sep=""))

  fungiOnly <- fungiOnly[fungiOnly[,4]!="d:Plantae",]
  print(paste("Now have ", length(fungiOnly[,1]), " taxa after removing plant otus found by taxonomy db 2", sep=""))
  
  
  if(length(grep("Chloroplast", fungiOnly[,3]))>0){
    fungiOnly= fungiOnly[-(grep("Chloroplast", fungiOnly[,3])),]
    print(paste("Now have ", length(fungiOnly[,1]), " taxa after removing cp otus from db 1", sep=""))
  }
  if(length(grep("Chloroplast", fungiOnly[,4]))>0){
    fungiOnly= fungiOnly[-(grep("Chloroplast", fungiOnly[,4])),]
    print(paste("Now have ", length(fungiOnly[,1]), " taxa after removing cp otus from db 2", sep=""))
  }

  #this keeps any row that for either database there is more than 10 characters
  #describing the taxonomic hypothesis.
  #this gets rid of anything that was not identified to phylum by
  #one or the other db.
  #IMPORTANT: this may very well get rid of target taxa if you are working in
  #a remote system. Use with caution.

  if(clean=="yes"){
    fungiOnly = fungiOnly[which(nchar(as.character(fungiOnly[,3])) > 10 | nchar(as.character(fungiOnly[,4])) > 10),]
  }

  ######################################################################################################
  #NOTE: I strongly recommend doing some spot checks of the OTUs that made it through. Find some that were
  #not id'd as fungi by both databases and go to the NCBI web interface and paste them in there.
  #if they show up as some other eukaryote then it may be worth scrubbing the data more
  ######################################################################################################

  #make new OTU table so that it includes only the taxa of interest
  otus2 <- otus[otus$OTU_ID %in% fungiOnly[,1],]

  #remove empty columns (no taxa of interest in a sample)
  otus2 <- otus2[,which(colSums(otus2[,2:length(otus2)])!=0)]


  otus2 <- t(otus2)
  colnames(otus2)  <- as.character(otus2[1,])
  otus2 <- otus2[-1,]
  row.names(otus2) <- gsub("X(\\d+\\.[a-h]\\d+).*","\\1",row.names(otus2))

  otus2=data.frame(otus2)

  write.csv(otus2,file="cleanOTUtable_youshouldrename.csv", row.names = F)

}

main()
