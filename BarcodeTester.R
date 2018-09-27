# source("https://bioconductor.org/biocLite.R")
# biocLite("DNABarcodes")
# for vignette see https://www.bioconductor.org/packages/devel/bioc/vignettes/DNABarcodes/inst/doc/DNABarcodes.html
library("DNABarcodes")
library("Biostrings")
#read possible barcodes
bcs <- read.csv("~/Desktop/barcodes.csv", header=F)

barcodes <- c(as.character(bcs[,1]), as.character(bcs[,3]))
#function to test for compliments

#sanitycheck code
barcodes[193]="tcgtacga"

for(i in 1:length(barcodes)){
  seq <- strsplit(as.character(barcodes[i]),"")
  
  complementer <- function(x){
    if(x == "a"){return("t")}
    if(x == "t"){return("a")}
    if(x == "g"){return("c")}
    if(x == "c"){return("g")}
  }  
  Compseq <- paste(lapply(seq[[1]], FUN=complementer), collapse="")
  
  out <- grep(Compseq,barcodes, value=T)
  if(length(out)>0){
    print("---")
    print(paste("barcode ", i))
    print(barcodes[i])
    print("has a match with ")
    print(out)
  }
}






len <- NA

for(i in 1:length(as.character(bcs$V1))){
  len[i] <- nchar(as.character(bcs$V3)[i])
}
  
bcs_analysis8 <- analyse.barcodes(as.character(bcs$V3)[which(len==8)],
                                 metric = c("hamming", "seqlev", "levenshtein"), 
                                 cores=2)
bcs_analysis9 <- analyse.barcodes(as.character(bcs$V1)[which(len==9)],
                                 metric = c("hamming", "seqlev", "levenshtein"), 
                                 cores=2)
bcs_analysis10 <- analyse.barcodes(as.character(bcs$V1)[which(len==10)],
                                 metric = c("hamming", "seqlev", "levenshtein"), 
                                 cores=2)

#test GC content
gc <- NA
for(i in 1:length(as.character(bcs$V1))){
gc[i] <- (sum(letterFrequency(DNAString(as.character(bcs$V1[i])), letters=c("C","G")))/nchar(as.character(bcs$V1[i])))
}
hist(gc)

#test for unwatnted primer binding 
#using compliments
#f 16s GTGYCAGCMGCCGCGGTA
for(i in 1:length(as.character(bcs$V1))){
print(grep(pattern= as.character(bcs$V1)[i], x = "CACCGTCGACGGCGCCAT", value=T))
}
for(i in 1:length(as.character(bcs$V1))){
  print(grep(pattern= as.character(bcs$V1)[i], x = "CACTGTCGACGGCGCCAT", value=T))
}
for(i in 1:length(as.character(bcs$V1))){
  print(grep(pattern= as.character(bcs$V1)[i], x = "CACCGTCGCCGGCGCCAT", value=T))
}
for(i in 1:length(as.character(bcs$V1))){
  print(grep(pattern= as.character(bcs$V1)[i], x = "CACTGTCGCCGGCGCCAT", value=T))
}

#reverse 16s
#GGACTACH[act]V[acg]GGGTW[at]TCTAAT
#m is a or c
#y is c or t

for(i in 1:length(as.character(bcs$V1))){
  print(grep(x= as.character(bcs$V1)[i], pattern = "CCTGAT", value=T))
}
for(i in 1:length(as.character(bcs$V1))){
  print(grep(x= as.character(bcs$V1)[i], pattern = "CCCA", value=T))
}
for(i in 1:length(as.character(bcs$V1))){
  print(grep(x= as.character(bcs$V1)[i], pattern = "AGATTA", value=T))
}


#f ITSs 
for(i in 1:length(as.character(bcs$V1))){
  print(grep(pattern= as.character(bcs$V1)[i], x = "GAACCAGTAAATCTCCTTCATT", value=T))
}

#reverse ITS
#GCTGCGTTCTTCATCGATGC
#m is a or c
#y is c or t

for(i in 1:length(as.character(bcs$V1))){
  print(grep(pattern= as.character(bcs$V1)[i], x = "CGACGCAAGAAGTAGCTACG", value=T))
}
