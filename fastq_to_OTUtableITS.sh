#!/bin/bash
#To cite anything here, see the USEARCH webpage and cite appropriately.
#I copied the pipeline from the USEARCH documentation and heavily edited
#J. Harrison Jan 16, 2018
#
#
# #IMPORTANT PLEASE READ!!!
# #This script will output a summary txt called Processing_Summary in a dir called out/
# #within the working dir. This has important QC info so check it out.
# #This script assumes one starts in the dir that has all unpacked fastqs to process
# #The script requires a dir to be present called db/ that has a comparison database of your chosing
# #You can change the path to this database in Step 3
#
# #You will want to make sure the primer stripping step (step 4) has the correct number of
# #bases specified to be removed
# #this will change depending on what primers you are using
#
# #This script assumes demultiplexed reads that are zipped (can comment out zipped option below)
#
# #YOU WILL NEED an extra script to be in the working dir.

# #rmv_fungal_matches.pl

# #After you run this script you will want to generate taxonomic hypotheses for the OTUs/ZOTUs generated here
# #Then you will want to remove unwanted taxa from the OTU table (e.g. host reads)

##############################
#What this script does
#1. Check and see if we have any samples with a lot of errors (this is set at 2 for now)
#2. merges your paired end reads
#3. checks that reads are oriented in the same way
#4. removes bp where the primer binds, as these are prone to errors in base calling
#5. Sequence filtering
#6. Find unique read sequences and abundances
#7. Make OTUs and ZOTUS
#8. Make an OTU table
#9. Calculate summary statistics and add those to summary file
#10. Clean up working dir
#############################


#here the script removes any directories called "out" in the working directory
#and makes a new directory called out
rm -rf out
mkdir out
touch out/Processing_Summary.txt

gunzip *fastq.gz

#######################################
####Step 1. Examine quality of reads
#######################################

#the main benefit to this step is to see how many reads get cut.
#really it can be skipped because of the filtering done later

mkdir fastq_info

for fq in *fastq
do
  usearch -fastx_info $fq -output fastq_info/${fq}info
done

#make a summary txt file that has number of expected errors per sample
grep "^EE" fastq_info/* > fastq_info/expected_error_summary.txt

#cut out just the floating point estimates of expected errors
grep -o "EE mean [0-9].[0-9]" fastq_info/expected_error_summary.txt | cut -d " " -f3 > EE.txt

counter=0
#initialize array of bad samples
badsamples[0]=""
while read p; do
	counter=$((counter+1))
	if (( $(echo "$p > 2" | bc -l) )); then
	 	toadd=`echo $counter`
	 	 badsamples+=(`echo $counter`)
	 	 echo $toadd
	 fi
done < EE.txt #note the strange way to pass an in file here

rm EE.txt

#get the length of the array. if more than one then we have a problem sample
numbad=`echo ${#distro[@]}`

#if problem samples are present, then tell us which those are
#if [ "$numbad > 1" ]; then
if (( $(echo "$numbad > 1" | bc -l) )); then
	echo "The following samples have more than two expected errors, perhaps keep an eye on them:" >> out/Processing_Summary.txt

	#for each element in this array do...
	for i in `echo "${badsamples[@]}"`
	do
		sed -n ${i}p fastq_info/expected_error_summary.txt >> out/Processing_Summary.txt
	done
else
	echo "No samples had more than 2 expected errors. yay!" >> out/Processing_Summary.txt
fi

#######################################
####Step 2. Merge reads
#######################################

#Merge paired reads
#Add sample name to read label (using the -relabel option).
#Sample name is just the file name
#this defaults to 10 threads, or number of cores, whichever is less
#I doubt you will need to run more cores here very often as this is not slow

for f in *R1*
do
	fname=$(basename $f)
	fname2=${fname/R1/R2}

	usearch -fastq_mergepairs ${fname} -reverse ${fname2} -fastqout ${fname}.merged.fq -sample $fname

done

#######################################
####Step 3. Check reads are oriented the same way
#######################################

#note they should be oriented in same way from MiSeq, but worth doublechecking. This may need to be reconsidered
#if we try other sequencing techniques for which I am unfamiliar

#Note that I am just checking for a few samples, as all samples should be the same
#note the database for comparison may need a new path, or be changed if you have a different db

cat `ls *merged.fq | head -n 5` > testFiles.fq

usearch -orient testFiles.fq -db db/rdp_its_v2.fa -tabbedout orient_out.txt

#this should show most reads are + or ? (the latter being ones that didn't match anything in the db)

echo "Sequence orientation for the first 5 samples." >> out/Processing_Summary.txt

orient=`eval cut -f2 orient_out.txt | sort | uniq -c`

echo "$orient" >> out/Processing_Summary.txt

rm -rf orient_out.txt
rm -rf testFiles.fq

#######################################
####Step 4. Remove primer binding region
#######################################

# Strip primers (V4F is 18, V4R is 20) (ITS1F is 22, ITS2R is 20)
# this will be important to doublecheck, as this can vary by primer pair

for f in *merged.fq
do
	usearch -fastx_truncate $f -stripleft 22 -stripright 20 -fastqout ${f}stripped.fq
done

#Note we don't do any more length trimming because we are using merged reads

#######################################
####Step 5. Sequence filtering
#######################################

#Note that we do this after merging so we can leverage info from both reads when calling bases
#If you are using a sequence technology prone to errors, then best to rethink this step
#or at least look into the output more carefully,
#see https://www.drive5.com/usearch/manual/pipe_readprep_filter.html

for f in *stripped.fq
do
	usearch -fastq_filter $f -fastq_maxee 1.0 -fastaout ${f}.filtered.fa
done

#######################################
####Step 6. Find unique read sequences and abundances
#######################################

cat *filtered.fa> combined_filtered.fa

usearch -fastx_uniques combined_filtered.fa -sizeout -fastaout uniqueSequences.fa

# #######################################
# ####Step 7. make OTUs and ZOTUs
# #######################################

## # Make 97% OTUs and filter chimeras
#note we remove singletons here

usearch -cluster_otus uniqueSequences.fa -otus otus970.fa -relabel Otu -minsize 2

# Denoise: predict biological sequences and filter chimeras
#also called ZOTUs, or ESVs..exact sequence variants
usearch -unoise3 uniqueSequences.fa -zotus zotuspre1.fa

#number of Otus before QC steps
echo "Number of OTUS (97% similarity threshold) before quality checks" >> out/Processing_Summary.txt
echo `grep ">" otus970.fa | wc -l` >> out/Processing_Summary.txt

echo "Number of ZOTUS before quality checks" >> out/Processing_Summary.txt
echo "Many of these will be non-target organism (i.e. host)" >> out/Processing_Summary.txt
echo `grep ">" zotuspre1.fa | wc -l` >> out/Processing_Summary.txt

###Clean up these OTU files in a few ways

#First cut out any OTUs that were very short as these are likely errors that acrued during PCR
#note that we write the before and after sums to know how many spurious OTUs there are

usearch -sortbylength otus970.fa -fastaout otus97pre2.fa -minseqlength 64
usearch -sortbylength  zotuspre1.fa -fastaout  zotuspre2.fa -minseqlength 64

#number of Otus before QC steps
echo "Number of 97OTUS that passed min length requirements: <64 bp long" >> out/Processing_Summary.txt
echo `grep ">" otus97pre2.fa | wc -l` >> out/Processing_Summary.txt

echo "Number of zotus that passed min length requirements: <64 bp long" >> out/Processing_Summary.txt
echo `grep ">" zotuspre2.fa | wc -l` >> out/Processing_Summary.txt

#Second, remove low complexity OTUS (ones filled with repeats)
#note I am not saving hits or anything here.

usearch -filter_lowc otus97pre2.fa -output otus97pre3.fa
usearch -filter_lowc  zotuspre2.fa -output  zotuspre3.fa

#number of Otus before QC steps
echo "Number of 97OTUS that passed min complexity reqs. (OTUS with >25% low complexity removed)" >> out/Processing_Summary.txt
echo `grep ">" otus97pre3.fa | wc -l` >> out/Processing_Summary.txt

echo "Number of zotus that passed min complexity reqs. (OTUS with >25% low complexity removed)" >> out/Processing_Summary.txt
echo `grep ">" zotuspre3.fa | wc -l` >> out/Processing_Summary.txt

#Remove any OTUs that are from PhiX (the Illumina internal standard). Likely won't lose much here for miseq
#may be more problematic with Hiseq...keep an eye on this

usearch -filter_phix otus97pre3.fa -output otus97.fa
usearch -filter_phix  zotuspre3.fa -output  zotus.fa

#number of Otus before QC steps
echo "Number of 97OTUS that passed PhiX removal" >> out/Processing_Summary.txt
echo `grep ">" otus97.fa | wc -l` >> out/Processing_Summary.txt

echo "Number of zotus that passed PhiX removal" >> out/Processing_Summary.txt
echo `grep ">" zotus.fa | wc -l` >> out/Processing_Summary.txt

rm -rf zotuspre[0-9].fa
rm -rf otus97pre[0-9].fa
rm -rf otus970.fa

#check for offset sequences
usearch -cluster_fast otus97.fa -id 0.97 -strand both -alnout otus.aln -show_termgaps -userout offsetcheck97.txt -userfields query+target+qstrand+qlo+qhi+tlo+thi
echo "If there is anything here then you should check the file offsetcheck97.txt and see which reads may be messed up" >> out/Processing_Summary.txt
echo `eval cut -f4 offsetcheck97.txt | grep -v "^1$"` >> out/Processing_Summary.txt
echo `eval cut -f6 offsetcheck97.txt | grep -v "^1$"` >> out/Processing_Summary.txt
mv offsetcheck97.txt out/

#do for zotus
usearch -cluster_fast zotus.fa -id 0.97 -strand both -alnout otus.aln -show_termgaps -userout offsetcheckZ.txt -userfields query+target+qstrand+qlo+qhi+tlo+thi
echo "If there is anything here then you should check the file offsetcheck97.txt and see which reads may be messed up" >> out/Processing_Summary.txt
echo `eval cut -f4 offsetcheckZ.txt | grep -v "^1$"` >> out/Processing_Summary.txt
echo `eval cut -f6 offsetcheckZ.txt | grep -v "^1$"` >> out/Processing_Summary.txt
mv offsetcheckZ.txt out/

rm -rf otus.aln

##################################################
# Step 8. Make OTU table
##################################################
#Note the maxrejects option for the otutab function slows down that function a lot bc it reduces heuristics
#it reduces the probability of errors however.
#if speed becomes an issue for some reason, then reconsider this option.

echo "About to cat this make take a minute, go get some coffee!"

cat *stripped.fq > allstrip.fq

#convert to fasta so that -otutab works right
usearch -fastq_filter allstrip.fq -fastaout allstrip.fa

#make 97% OTU table
usearch -otutab allstrip.fa -otus otus97.fa -maxrejects 1000 -otutabout out/otuTable97otus.txt -notmatched unmapped97.fa

#remove the hash from the first line. For some reason usearch puts this there
#note the -i.bak  flag which makes edits inline

sed 's/#//' out/otuTable97otus.txt > interim.txt
sed 's/OTU\ ID/OTU_ID/' interim.txt > out/otuTable97otus.txt
rm -rf interim.txt

#bug in old usearch needs this hack
sed 's/Zotu/Otu/' zotus.fa > zotus2.fa
rm -rf zotus.fa
mv zotus2.fa zotus.fa

#make ZOTU table
usearch -otutab allstrip.fa -otus zotus.fa -maxrejects 1000 -otutabout out/otuTableZotus.txt -notmatched unmappedZOTUS.fa

#remove the hash from the first line. For some reason usearch puts this there
sed 's/#//' out/otuTableZotus.txt > interim.txt
sed 's/OTU\ ID/OTU_ID/' interim.txt > out/otuTableZotus.txt

#count the number of unassigned sequences and put this in the summary files
echo "Number of reads that didn't get mapped when making the 97 OTU table " >> out/Processing_Summary.txt
echo `eval grep "^>" unmapped97.fa | wc -l` >> out/Processing_Summary.txt

echo "Number of reads that didn't get mapped when making the Z OTU table " >> out/Processing_Summary.txt
echo `eval grep "^>" unmappedZOTUS.fa | wc -l` >> out/Processing_Summary.txt

echo "If these numbers are low, then you some unmapped reads can probably be ignored" >> out/Processing_Summary.txt
echo "However, if a lot of reads didn't map, then you will need to dig a bit to figure out why." >> out/Processing_Summary.txt
echo "For guidance see the USEARCH website, which has some good tips. " >> out/Processing_Summary.txt

rm -rf interim.txt
rm -rf allstrip.fq

#additional quality control step. Make sure all OTUs made it into the table
#check these missing OTUs, if there are any, to see what they are.
#could be an error with the otutable function, could be an offset read that didnt get caught
#could also be a strand duplicate (reverse compliment)

cut -f1 out/otuTable97otus.txt > included97s.txt
grep "^>" otus97.fa | sed "-es/>//" > otuTitles97.txt
sort otuTitles97.txt included97s.txt included97s.txt | uniq -u > missing_labels97.txt
usearch -fastx_getseqs otus97.fa -labels missing_labels97.txt -fastaout out/missing97otus.fa

#here we find the notmissing OTUs then match the missing OTUs against them, this helps check for
#offsets or strand dupes

sort missing_labels97.txt missing_labels97.txt otuTitles97.txt | uniq -u > notmissing_labels97.txt
usearch -fastx_getseqs otus97.fa -labels notmissing_labels97.txt -fastaout notmissing.fa
usearch -usearch_global out/missing97otus.fa -db notmissing.fa -strand both -id 0.97 -uc missingVsNotmissing97.uc -alnout missnot.aln

#do above for zotus

cut -f1 out/otuTableZotus.txt > includedZotus.txt
grep "^>" zotus.fa | sed "-es/>//" > otuTitlesZ.txt
sort otuTitlesZ.txt includedZotus.txt includedZotus.txt | uniq -u > missing_labelsZ.txt
usearch -fastx_getseqs zotus.fa -labels missing_labelsZ.txt -fastaout out/missingZotus.fa

sort missing_labelsZ.txt missing_labelsZ.txt otuTitlesZ.txt | uniq -u > notmissing_labelsZ.txt
usearch -fastx_getseqs zotus.fa -labels notmissing_labelsZ.txt -fastaout notmissingZ.fa
usearch -usearch_global out/missingZotus.fa -db notmissingZ.fa -strand both -id 0.97 -uc missingVsNotmissingZ.uc -alnout missnot.aln

unassignedreads=`wc -l missingVsNotmissing*.uc`

if (( $(echo "$unassignedreads > 1" | bc -l) )); then
  echo "Check out the missingVsNotmissing* files to see which OTUs were left out of the OTU table." >> out/Processing_Summary.txt
fi

rm -rf included*
rm -rf otuTitles*
rm -rf missing_labelsZ.txt
rm -rf missing_labels97.txt
rm -rf notmissing_labels97.txt
rm -rf notmissing_labelsZ.txt
rm -rf notmissing*
mv missingVs* out/
mv unmapped* out/
rm -rf *aln

 ##################################################
 # Step 9. Summary info calculation
 ##################################################

#Count number of original reads
echo "Original number of forward raw reads (unmerged, unfiltered)" >> out/Processing_Summary.txt
echo `grep "@M0" *R1*.fastq | wc -l` >> out/Processing_Summary.txt

echo "Original number of reverse raw reads (unmerged, unfiltered)" >> out/Processing_Summary.txt
echo `grep "@M0" *R2*.fastq | wc -l` >> out/Processing_Summary.txt

#number of reads that merged successfully
echo "Number of reads that merged successfully" >> out/Processing_Summary.txt
echo `grep "^@" *merged.fq | wc -l` >> out/Processing_Summary.txt

#number of reads that passed filtering
echo "Number of merged reads that passed filtering" >> out/Processing_Summary.txt
echo `grep "$>" *combined.filtered.fa | wc -l` >> out/Processing_Summary.txt

#number of unique sequences
echo "Number of unique sequences" >> out/Processing_Summary.txt
echo `grep ">" uniqueSequences.fa | wc -l` >> out/Processing_Summary.txt

#number of Otus
echo "Number of OTUS (97% similarity threshold)that passed QC" >> out/Processing_Summary.txt
echo "Many of these will be non-target organism (i.e. host)" >> out/Processing_Summary.txt
echo `grep ">" otus97.fa | wc -l` >> out/Processing_Summary.txt

echo "Number of ZOTUS that passed QC" >> out/Processing_Summary.txt
echo "Many of these will be non-target organism (i.e. host)" >> out/Processing_Summary.txt
echo `grep ">" zotus.fa | wc -l` >> out/Processing_Summary.txt



# # ##################################################
# # # Step 10. Clean up files
# # ##################################################

 rm -rf *merged*

 gzip *fastq
 mkdir rawreads

 mv *fastq.gz rawreads
 mv otus97.fa out/
 mv zotus.fa out/

 rm -rf all_ucs*uc
 rm -rf uniqueSequences.fa
 rm -rf combined_filtered.fa

 tar -czf fastq_info.tar.gz  fastq_info/
 rm -rf fastq_info/
