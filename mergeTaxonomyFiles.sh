#!/bin/bash

#Combine two files of taxonomy hypotheses output from Sintax
#this can be useful if you use two databases for taxonomy assignment

#Usage bash mergeTaxonomyFiles.sh file1 file2

#make sure input files are in same order. For some reason sintax doesn't do this
sort -k 1 $1 > tmp.txt
sort -k 1 $2 > tmp2.txt

#get the fourth field, which is the taxonomic assignment that was chosen at
#>80% probability, based on kmer matches
cut -f4 tmp.txt > tax.txt
cut -f4 tmp2.txt > tax2.txt

#get OTU names, these should be the same, but will write both just in case
cut -f1 tmp.txt > namesFile1.txt
cut -f1 tmp2.txt > namesFile2.txt

paste namesFile1.txt namesFile2.txt tax*txt > combinedTaxonomy.txt

rm -rf tmp*
rm -rf tax*txt
rm -rf namesFile*
