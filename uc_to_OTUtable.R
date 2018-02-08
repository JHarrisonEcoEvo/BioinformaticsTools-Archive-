#!/usr/bin/Rscript
#J. G. Harrison
#Jan. 16, 2018

#Usage (input file followed by its length, which is some number): 
#Rscript ProcessOTUtable.R  allucs.uc

#to see what the uc format means click on this:
#http://drive5.com/usearch/manual/opt_uc.html

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

process <- function(filename) {
  otus <- read.table(file = filename, header = F)
  otuTable = t(table(otus$V10,otus$V9))
  write.csv(otuTable, file=paste("./out/OTUtable_", filename, sep=""), row.names=T)
}

main()  

