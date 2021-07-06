#!/usr/bin/python3
#J. Harrison, Univ. of Wyoming

#Program purpose: take an input fasta file and match a sequence to accessions within that file
#Usage: matcher.py query.fasta target.fasta Y/N 

#the Y/N switch selects if the Biopython pairwise alignment tool is desired or not. 

#IMPORTANT: this is a simple program and cannot handle any variation from the 
#following expected file formats: the query sequence must have a header line
#and then the next line should be the target sequence, with no extraneous line endings. 
#Target.fasta should be a fasta file of any length. 
#For help contact J. Harrison, particularly if you are unfamiliar with fasta
#formats and what line endings are. 
#This script uses exact sequence matching, so it cannot handle any sort of mutation
#in either the target or query sequence. 

#For debugging
#target_genome = open("/Users/joshuaharrison/git_repo/bioinformaticsTools/testgenome", "r")


#****************#
#Needed improvements include:
#Going through the input fastas and selecting just the sequence lines
#E.g., find a > then take all lines until the next >. 
#strip line endings from them and concatenate them. 

##################
# Begin programe #
##################

#load modules
import sys

#First get input data from STDIN
query_seq_fn = sys.argv[1]
target_genome_fn = sys.argv[2]
biopython_switch = sys.argv[3]

if len(sys.argv) <=1:
    print("This program expects a query FASTA file and a target FASTA file. It doesn't seem like the paths to those two files were provided.")
#Parse command line
if query_seq_fn:
    query_seq = open(query_seq_fn, "r")
    print ("Successfully opened query sequence " + query_seq_fn)
else:
    print ("No query specified")
if target_genome_fn:
    target_genome = open(target_genome_fn, "r")
    print ("Successfully opened target genome " + target_genome_fn)
else:
    print ("No target genome specified")
 
#Read query line and remove line endings, then print that sequence to STDOUT for the user
query_seq.readline()
query = str.strip(query_seq.readline())
print("This is the query sequence: " + query)

hits = []
for line in target_genome:
    #print(line)
    if query in line:
        hits.append(line)
        break
    
if len(hits) >= 1:
    print("Yay, the query sequence was in the target genome")
    print("NOTE: this is an exact search and therefore any unexpexted sequence variation or file formatting will cause this program to break")


#If we want to do a more robust alignment we can use tools from Biopython
if biopython_switch == "y":
    from Bio import pairwise2
    # Import format_alignment method
    from Bio.pairwise2 import format_alignment
    
    
    # Do local alignment where we look for best alignment over subsets of the sequences
    #being compared.
    # No parameters. Identical characters have score of 1, else 0.
    # No gap penalties.
    #Modified from: https://towardsdatascience.com/pairwise-sequence-alignment-using-biopython-d1a9d0ba861f
    #That page has some examples of how to employ gap penalties and so on. 
    
    alignments = []
    for line in target_genome:    
        alignments = pairwise2.align.localxx(query, line)
    
    # Display alignments in a nicer format
    for a in alignments:
        print(format_alignment(*a))
    
    #extract top alignment and print to screen
    scores = []
    for a in alignments:
        scores.append(a.score)
    
    #Find the best match score and the index of that score. 
    #We will use the latter to say which line in the input genome 
    #was the best match
    best_match = max(scores)
    best_match_index = scores.index(best_match)
    
    
    print("The following was the best alignment: \n" + str(alignments[best_match_index]))

target_genome.close()
query_seq.close()