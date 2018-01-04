#!/usr/bin/Rscript
#J. G. Harrison
library("edgeR")
#This script takes an OTU table (csv) as input where rows are otus, 
#columns are samples, and cells are filled with raw read counts. Note that if you get weird results
#you probably have your input matrix transposed!
#The output will include a TMM normalized OTU table (as per Mcurmdie and Holmes 2012)
#and several rarified tables (1k, 5k, 10k reads), note in many cases you may
#want to try analyses without samples with low coverage, these won't be rarefied of course.

#Usage: 
#Rscript ProcessOTUtable.R  OTUtable.csv OTUtable2.csv ...

#for more info on using R outside of an IDE see https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/

#read stdin and argument to rarefy. Note that this allows for piping, and for multiple files
main <- function() {
  inargs <- commandArgs(trailingOnly = TRUE) #cuts out the calls to R that Rscript includes automatically, eg --slave --no-restore, etc
  print(inargs)
  if (length(inargs) == 0) {
    process(file("stdin"))
  } else {
    for (f in inargs) {
      process(f)
    }
  }
}

#Perform work on the input file(s)
process <- function(filename) {
  dat <- read.csv(file = filename, header = T) #here is where you could swap to tab separated.
  
  #Normalize via the TMM method, RLE often fails with sparse data
  fungi = DGEList(counts=as.matrix(t(dat[,2:length(dat)])), group=as.matrix(dat[,1]))
  
  y=calcNormFactors(fungi, method="TMM")
  
  #get normalized counts
  nc=cpm(y, normalized.lib.sizes=T)
  colnames(nc) = dat[,1]
  
  #make a filename that only has the file, not its whole path, this aids naming the out file
  abrv_filename = basename(filename)
  write.csv(t(nc), file=paste("./TMMnormalized_", abrv_filename, sep=""))
  
  rr  = vegan::rrarefy(t(dat[,2:length(dat)]), sample=1000)
  colnames(rr) = dat[,1]
  write.csv(rr, file=paste("./1krarefied_", abrv_filename, sep=""))
  
  rr  = vegan::rrarefy(t(dat[,2:length(dat)]), sample=5000)
  colnames(rr) = dat[,1]
  write.csv(rr, file=paste("./5krarefied_", abrv_filename, sep=""))
  
  rr  = vegan::rrarefy(t(dat[,2:length(dat)]), sample=10000)
  colnames(rr) = dat[,1]
  write.csv(rr, file=paste("./10k_rarefied", abrv_filename, sep=""))
  
  rr  = vegan::rrarefy(t(dat[,2:length(dat)]), sample=20000)
  colnames(rr) = dat[,1]
  write.csv(rr, file=paste("./20k_rarefied", abrv_filename, sep=""))
}

main()  


