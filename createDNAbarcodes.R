# source("https://bioconductor.org/biocLite.R")
# biocLite("DNABarcodes")
# for vignette see https://www.bioconductor.org/packages/devel/bioc/vignettes/DNABarcodes/inst/doc/DNABarcodes.html
library("DNABarcodes")

#use the modified Sequence Levenshtein distance so we can find indels and up to one substittuion

mySet <-create.dnabarcodes(8, metric="seqlev", heuristic="ashlock", cores=4)
length(mySet)

#here you can use Hamming dist.
#mySet <- create.dnabarcodes(12, dist=5) 
#dist 5 means we can correct up to two substitutions
#length(mySet)

#dont forget to include the illumina adaptor to oligo if that is needed for your application