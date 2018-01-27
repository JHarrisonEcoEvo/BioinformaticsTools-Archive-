#!/bin/bash

#Combine two files of taxonomy hypotheses output from Sintax
#this can be useful if you use two databases for taxonomy assignment

#Usage bash mergeTaxonomyFiles.sh file1 file2

#make sure input files are in same order. For some reason sintax doesn't do this
sort -k 1 $1 > $1
sort -k 1 $2 > $2

#get the fourth field, which is the taxonomic assignment that was chosen at
#>80% probability, based on kmer matches
cut -f4 $1 > tmp.txt
cut -f4 $2 > tmp2.txt

#get OTU names, these should be the same, but will write both just in case
cut -f1 $1 > namesFile1.txt
cut -f1 $2 > namesFile2.txt

paste namesFile1.txt namesFile2.txt tmp.txt tmp2.txt > combinedTaxonomy.txt

rm -rf tmp*
rm -rf namesFile*
