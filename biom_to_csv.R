#!/usr/bin/Rscript
#J. G. Harrison

#script to convert from the heniously hard to work with biom format to something more memory intensive
#but that works with more analytical methods
#Do you detect some frustration with the biom format??

#usage
#Rscript biom_to_csv.R yourfile

library(biomformat)

#read stdin. Note that this allows for piping, and for multiple files
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
  datbiom <- read_biom(filename)
  x <- biom_data(datbiom)
  yuck <- as.matrix(x)
  write.csv(yuck, file=paste(filename, ".csv", sep=""))
  write.csv(observation_metadata(datbiom), file=paste(filename, "taxonomy.csv", sep=""))
}

main()  

