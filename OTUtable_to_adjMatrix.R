#!/usr/bin/Rscript
#J. G. Harrison

#Take an OTU table and output an adjacency matrix 

#This script takes an OTU table (csv) as input where columns are otus, 
#rows are samples, and cells are filled with raw read counts.
#Note that if you get weird results
#you probably have your input matrix transposed! 
#There should not be an indexing column as the first column in the input table

#Usage (can pipe and take multiple inputs): 
#Rscript OTUtable_to_adjMatrix.R  OTUtable.csv OTUtable2.csv ...

#for more info on using R outside of an IDE see https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/

####################
#Describe functions#
####################

main <- function() {
  inargs <- commandArgs(trailingOnly = TRUE) 
  #above cuts out the calls to R that Rscript includes automatically, 
  #eg --slave --no-restore, etc
  
  print(inargs)
  if (length(inargs) == 0) {
    process(file("stdin"))
  } else {
    for (f in inargs) {
      process(f)
    }
  }
}

#perform work on input file
process <- function(filename) {
  dat <- read.csv(file = filename, header = T) 
  
#get rid of taxa that are not in the data, and samples with no taxa
dat = dat[, colSums(dat) > 0]
dat = dat[rowSums(dat) > 0, ]

#make empty adjacency matrix, note that it is a square matrix of dimensions defined by
#number of taxa (columns in input otu table)

adjmat = matrix(0, nrow = ncol(dat),
                ncol = ncol(dat))

#find those columns that are positive for some taxon i
#these are the taxa that cooccur with i
#and turn these to ones.

for (i in 1:ncol(dat)) {
  adjmat[which(colSums(dat[dat[, i] > 0, ]) > 0), i] = 1
}

names(adjmat) = names(dat)
row.names(adjmat) = names(dat)

#make a filename that only has the file, not its whole path, this helps with naming the out file
abrv_filename = basename(filename)

write.csv(adjmat, file=paste("./Adj_matrix_", abrv_filename, sep=""), row.names=F)
}

#call wrapper function
main()  

