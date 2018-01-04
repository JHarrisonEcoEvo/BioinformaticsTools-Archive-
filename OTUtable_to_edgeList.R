#!/usr/bin/Rscript
#J. G. Harrison

#Take an OTU table and output an edge list for use in network generation

#This script takes an OTU table (csv) as input where rows are otus, 
#columns are samples, and cells are filled with raw read counts. Note that if you get weird results
#you probably have your input matrix transposed! There should not be an indexing column as the first column in the input table

#Usage: 
#Rscript GenerateEdgeList.R  OTUtable.csv OTUtable2.csv ...

#for more info on using R outside of an IDE see https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/

####################
#Describe functions#
####################


#This is the main wrapper function. It calls other functions
#we will call this function at the end of the script
#read stdin and arguments
#Note that this allows for piping, and for multiple files

main <- function() {
  inargs <- commandArgs(trailingOnly = TRUE) 
  #above cuts out the calls to R that Rscript includes automatically, eg --slave --no-restore, etc
  
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
#This function will generate edge lists and write them to the working dir

process <- function(filename) {
  dat <- read.csv(file = filename, header = T) #here is where you could swap to tab separated.

  #sanity check to see if rows outnumbers columns.
  #This will pretty much never happen if input matrix is oriented properly (OTUs>samples)
  #The script will proceed anyway however, just with a warning.
  if(dim(dat)[1] > dim(dat)[2]){
    cat("The number of rows is more than the number of columns\nAre you sure that you have OTUs as columns, and samples as rows?\nProceeding anyway\n")
  }
  
  #split OTU table up into a list where each row is a list element, this should speed up operations a lot
  dat_rows = split(dat, row.names(dat))
  
  #make a vector where we will put our edges, not the to/from is for ease and does NOT imply directionality
  edges=NA
  #make a vector for the number of times taxa cooccur
  cooccurCounts=NA
  
  #note the <<- operator to write to variables in the environment (the ones we just defined specifically)

  lapply(
    dat_rows,
    FUN = function(x) {
      #first cut to just those taxa that are present, and get their names
      taxa_present = names(x)[which(x > 0)]
      
      for (i in taxa_present) {
        for (j in taxa_present) {
          #the two possible edges
          edgeij = paste(i, j, sep = "")
          edgeji = paste(j, i, sep = "")
          
          #if this cooccurrence has not yet been observed, then add it to the edges vector,
          #and add a 1 to that element of the cooccurCounts vector
          #and go to the next iteration of the for loop
          #if it has been observed, then add to counter
          #this counts everything twice, so outside this loop we will divide out counting vector by two
          
          if (edgeij %in% edges == FALSE & edgeji %in% edges == FALSE) {
            #this inner conditional keeps us from making self loops
            if (edgeij != edgeji) {
              edges <<- append(edges, edgeij)
              cooccurCounts[length(edges)] <<- 1
            } else if (edgeij == edgeji) {
              next
            }
          } else if (edgeij %in% edges == TRUE &
                     edgeji %in% edges == FALSE) {
            cooccurCounts[grep(edgeij, edges)] <<- cooccurCounts[grep(edgeij, edges)] + 1
          } else if (edgeij %in% edges == FALSE &
                     edgeji %in% edges == TRUE) {
            cooccurCounts[grep(edgeji, edges)] <<- cooccurCounts[grep(edgeji, edges)] + 1
          }
        }
      }
    }
  )

  #remove the NAs that started our vectors (note we do this at the end for cooccur)
  edges = edges[-1]

  #here is where we account for double counting cooccurrences.
  cooccurCounts = cooccurCounts/2
  
  #split up our edges vector and turn it into something that could be input into a network pckg
  edgSplit = strsplit(edges, split="X")
  
  #make vectors of our final edges
  to=NA
  from=NA
  
  lapply(
    edgSplit,
    FUN = function(x) {
      to <<- append(to, x[2])
      from <<- append(from, x[3])
    }
  )
  
  #make a filename that only has the file, not its whole path, this helps with naming the out file
  abrv_filename = basename(filename)
  out = na.omit(data.frame(to, from, cooccurCounts))
  write.csv(out, file=paste("./Edgelist_", abrv_filename, sep=""), row.names=F)
}

#call wrapper function
main()  

