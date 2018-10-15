#!/usr/bin/python
#J. Harrison, Univ. of Wyoming
#Oct 11, 2018

######################
# Overview and usage #
######################
#This script parses a fastq file by barcode sequence. It takes both forward and
#reverse sequences. It assumes the input sequences are of the form:
#5' barcode - primer region - template
#You will need to know the expected primer sequence to operate this program.
#Barcodes should be in csv files with the forward barcode in the first field and
#reverse barcode in the second field, sample name in the third field, forward
#primer in the fourth field, and reverse primer in the fifth field. Note that
#the fields with primer sequences don't have to be filled, just having the
#primer in the first cell is fine. IMPORTANT: if you are using primers with
#degenerate bases then spell out the possible true bases using IUPAC codes
#within brackets for each degenerate position. For instance, spell out Y as [CT]
#Finally, unless you want to enter a world of pain don't make your barcode
#file with Excel, use google sheets or libreoffice or some other non-medieval
#program.

#Levenshtein distances are used to differentiate barcodes that have errors in
#base calling.

#USAGE
#python demux.py <forward_reads> <reverse_reads> <barcodes>

#Likely ways this program could fail:

#Won't work if are using fasta files instead of fastq files.
#Will fail if you have a barcode file with weird line endings or Excel-
#imposed prefixes.
#Will fail, if your fastq files do not have 4 lines per sample (they should)
#If you run this script twice in the same place you will double the number of
#reads you have because it appends to out files. Doublecheck the final number
#of reads in out files to make sure expectations are met.

#####################
# How program works #
#####################

#1. Reads one line of input fastq at a time and finds the primer binding region
#   (mismatches are allowed, so long as they are few).
#2. Removes bases that come before the primer region (towards 5')
#3. Checks if those bases match anything in the barcode file
#4. If so, then the following sequence information and quality scores are
#   appended to an outfile for that barcode (called by the sample name )
#5. If the barcode doesn't match, then we assume it has a base calling error and
#   use the Levenshtein module to figure out which barcode is the closest match.
#   We use the L. distance because it works on strings of unequal length, unlike
#   the Hamming distance. We only correct for one error (to edit this see
#    comment in code below). If the match is bad, then we write the sequence to
#   a file corresponding to unknown samples.

#load modules
import sys #system stuff
import re #regex
import Levenshtein

#First get input data from STDIN
fn1 = sys.argv[1]
fn2 = sys.argv[2]
fn3 = sys.argv[3]

#Parse command line
if fn1:
    forwardreads = open(fn1, "r")
    print ("Successfully opened forward reads")
else:
    print ("No forward reads specified")
if fn2:
    reversereads = open(fn2, "r")
    print ("Successfully opened reverse reads")
else:
    print ("No reverse reads specified")
if fn3:
    barcodes = open(fn3, "r")
    print ("Successfully opened barcodes")
else:
    print ("No barcodes specified")

#Extract forward and reverse primers
barcodes = open(fn3, "r")
barcodeLine1 = barcodes.readline()
forwardPrimer = barcodeLine1.split(",")[3] #recall that python starts at 0
reversePrimer = barcodeLine1.split(",")[4]
reversePrimer = reversePrimer.rstrip("\n")

#Build lists of barcodes and samples
key_f = []
key_r = []
keybarcodes = {}


for line in open(fn3, "r"):

    key_f.append(line.split(',')[0])
    key_r.append(line.split(',')[1])
    bothbc = line.split(',')[0] + line.split(',')[1]
    keybarcodes[bothbc] = line.split(',')[2]

#Loop through reads, parse by barcode, and append reads to correct out file
for line in forwardreads:
    #Here we see if the input line is an Illumina header
    #If not, then we don't do anything.
    if re.match('^@.+:.+:', line) is not None:
        header = line
        read = forwardreads.readline()
        line3 = forwardreads.readline()
        quality = forwardreads.readline()

        #Compose regex to find primer region, use that as an anchor, and
        #capture the preceeding characters, which is the barcode.
        pattern = "(.+)" + forwardPrimer

        fb = re.findall(pattern, read)

        #Compute Levenshtein distances between the barcode and all options
        #find best option and if it is a good option, then
        #write to out file. If a poor match, then write to a file of reads
        #of unknown provenance.
        if not fb:
            fb = "unknown"
        else:
            if fb[0] not in key_f:
                levdist = []
                levkeys = []
                for k in key_f:
                    levkeys.append(k)
                    levdist.append(Levenshtein.distance(fb[0], k))
                    mindex = levdist.index(min(levdist))

                    #Find the minimum value in this list of distances
                    #and if this value <1 then write to corresponding out
                    #file. This means that one changes would need to be made in the
                    #query sequence to match an expected barcode.

                    if min(levdist) <= 1:
                        fb = key_f[mindex]
                    else:
                        fb = "unknown"
            else:
                fb = fb[0]

    #Now extract the rear barcode
    header_r = reversereads.readline()
    read_r = reversereads.readline()
    line3_r = reversereads.readline()
    quality_r = reversereads.readline()

    pattern_r = "(.+)" + reversePrimer

    rb = re.findall(pattern_r, read_r)

    if not rb:
        rb = "unknown"
    else:
        if rb[0] not in key_r:
            levdist = []
            levkeys = []
            for k in key_r:
                levkeys.append(k)
                levdist.append(Levenshtein.distance(rb[0], k))
                mindex = levdist.index(min(levdist))

                #Find the minimum value in this list of distances
                #and if this value <1 then write to corresponding out
                #file. This means that one changes would need to be made in the
                #query sequence to match an expected barcode.

                if min(levdist) <= 1:
                    rb = key_r[mindex]
                else:
                    rb = "unknown"
        else:
            rb = rb[0]

    #Now paste the forward and reverse barcodes together and see if they
    #are in our dictionary.
    bcs = fb + rb
    #print(bcs)
    if bcs in keybarcodes.keys():
        samp = keybarcodes[bcs]
        out_file = open(samp + "_forward", "a")
        out_file.write(header)
        #Remove the barcode
        basesTormv = len(fb)
        out_file.write(read[basesTormv:])
        out_file.write(line3)
        out_file.write(quality)
        out_file.close()

        out_file = open(samp + "_reverse", "a")
        out_file.write(header_r)
        #Remove the barcode
        basesTormv = len(rb)
        out_file.write(read_r[basesTormv:])
        out_file.write(line3_r)
        out_file.write(quality_r)
        out_file.close()

    else:
        out_file = open("unplaced_forward", "a")
        out_file.write(header)
        out_file.write(read)
        out_file.write(line3)
        out_file.write(quality)
        out_file.close()

        out_file = open("unplaced_reverse", "a")
        out_file.write(header_r)
        out_file.write(read_r)
        out_file.write(line3_r)
        out_file.write(quality_r)
        out_file.close()
