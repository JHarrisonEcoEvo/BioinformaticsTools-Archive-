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

#For examples of input data see the dir Example_data/ in the git repo

#Usage: Rscript useTax_to_cleanOTUtable.R combinedTaxonomyfile.txt OTUtableExample.txt Well_to_sampleKeyExample.csv clean otusToRemove.txt

#where clean is either "yes" or "no" and determines if you want
#taxa removed that were not assigned to a phylum
#otusToRemove.txt - is a txt file with otus that were identified as mt or cpdna
#that you want removed (this is the output of the check16s* script).

main <- function() {
  inargs <- commandArgs(trailingOnly = TRUE)
  print(inargs)

  #input otu table, sample-well key, and taxonomy file
  tax=read.delim(file=inargs[1], header=F)
  otus=read.delim(file=inargs[2], header=T)
  sampNames=read.csv(file=inargs[3])
  sampNames$pl_well=paste(sampNames$plate, sampNames$well, sep=".")
  clean=inargs[4]
  toRmv=read.delim(inargs[5], header=F)

  #now we need to pick which OTUs to keep
  #note the use of fungi in the object names below doesn't matter. This
  #script will work for bacteria taxonomy files too.
  #I am first removing any OTUs that the UNITE database said were plants.
  #then I will subset the data to just keep only those OTUs that one or the other database
  #had >80% confidence in placement to phylum. In previous work I have found that <80% placement to
  #phylum often match nontarget taxa on NCBI...so I remove those too.

#note these are just operating on the 4th field because the UNITE database out has a different format
  print(paste("Started with ", length(tax[,1]), " taxa", sep=""))

  fungiOnly= tax[tax[,4]!="d:Plantae",]
  print(paste("Now have ", length(tax[,1]), " taxa after removing plant otus", sep=""))

  fungiOnly= tax[-(grep("Chloroplast", tax[,4])),]
  print(paste("Now have ", length(tax[,1]), " taxa after removing cp otus", sep=""))

  fungiOnly= tax[tax[,1] != toRmv[,1]]
  print(paste("Now have ", length(tax[,1]), " taxa after removing input cp and mt OTUs", sep=""))

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
  otus2=otus[otus$OTU_ID %in% fungiOnly$V1,]

  #remove empty columns (no taxa of interest in a sample)
  otus2=otus2[,which(colSums(otus2[,2:length(otus2)])!=0)]

  otus2=t(otus2)
  colnames(otus2)=as.character(otus2[1,])
  otus2=otus2[-1,]
  row.names(otus2)=gsub("X(\\d+\\.[a-h]\\d+).*","\\1",row.names(otus2))

  #check to see if anything in the row.names is not in the key
  if(length(which(!(row.names(otus2) %in% sampNames$pl_well)))>1){
    print("Error: the following row names of the output otu file are not in the key")
    print("This means they had a non-parsable name and so the merge failed")
    print("between the sample names and the well location. This is ok, just")
    print("rename by hand, or tweak the call to gsub.")
    print(row.names(otus2)[which(!(row.names(otus2) %in% sampNames$pl_well))])
    print("These will be at the bottom of the OTU table")
  }

  otus2=data.frame(otus2)

  #replace well based row names with the sample names
  otus2$samps = row.names(otus2)
  newotus = merge(otus2, sampNames, by.x="samps", by.y="pl_well", all.x = TRUE)
  outputfile = data.frame(newotus$sampleName, newotus[,grep("Otu",names(newotus))])

  write.csv(outputfile,file="cleanOTUtable_youshouldrename.csv", row.names = F)

}

main()
