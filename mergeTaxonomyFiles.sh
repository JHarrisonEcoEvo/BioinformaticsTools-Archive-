#!/bin/bash

#Combine two files of taxonomy hypotheses output from Sintax
#this can be useful if you use two databases for taxonomy assignment

#Usage bash mergeTaxonomyFiles.sh file1 file2

cut -f4 $1 > tmp.txt
cut -f4 $2 > tmp2.txt

paste tmp.txt tmp2.txt > combinedTaxonomy.txt

rm -rf tmp*
