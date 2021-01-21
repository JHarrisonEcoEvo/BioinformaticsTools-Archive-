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
# #If you are performing 97% clustering then you will want to make sure the primer stripping step (step 4) has the correct number of
# #bases specified to be removed
# #this will change depending on what primers you are using
#
# #This script assumes demultiplexed reads that are zipped (can comment out zipped option below)
#
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

#vsearch can handle zipped files, but usearch can't
#in a real cursory test usearch did a better job merging, so unzip
gunzip *fastq.gz

#######################################
####Step 1. Examine quality of reads
#######################################

#the main benefit to this step is to see how many reads get cut.
#it can be skipped because of the filtering done later
#
# mkdir fastq_info
#
# for fq in *fastq
# do
#   usearch -fastx_info $fq -output fastq_info/${fq}info
# done
#
# #make a summary txt file that has number of expected errors per sample
# grep "^EE" fastq_info/* > fastq_info/expected_error_summary.txt
#
# #cut out just the floating point estimates of expected errors
# grep -o "EE mean [0-9].[0-9]" fastq_info/expected_error_summary.txt | cut -d " " -f3 > EE.txt
#
# counter=0
# #initialize array of bad samples
# badsamples[0]=""
# while read p; do
# 	counter=$((counter+1))
# 	if (( $(echo "$p > 2" | bc -l) )); then
# 	 	toadd=`echo $counter`
# 	 	 badsamples+=(`echo $counter`)
# 	 	 echo $toadd
# 	 fi
# done < EE.txt #note the strange way to pass an in file here
#
# rm EE.txt
#
# #get the length of the array. if more than one then we have a problem sample
# numbad=`echo ${#distro[@]}`
#
# #if problem samples are present, then tell us which those are
# #if [ "$numbad > 1" ]; then
# if (( $(echo "$numbad > 1" | bc -l) )); then
# 	echo "The following samples have more than two expected errors, perhaps keep an eye on them:" >> out/Processing_Summary.txt
#
# 	#for each element in this array do...
# 	for i in `echo "${badsamples[@]}"`
# 	do
# 		sed -n ${i}p fastq_info/expected_error_summary.txt >> out/Processing_Summary.txt
# 	done
# else
# 	echo "No samples had more than 2 expected errors. yay!" >> out/Processing_Summary.txt
# fi

#######################################
####Step 2. Merge reads
#######################################

#Merge paired reads
#Add sample name to read label (using the -relabel option).
#Sample name is just the file name
#this defaults to 10 threads, or number of cores, whichever is less
#I doubt you will need to run more cores here very often as this is not slow
#
# for f in *forward*
# do
# 	fname=$(basename $f)
# 	fname2=${fname/forward/reverse}
#
# #usearch -fastq_mergepairs ${fname} -reverse ${fname2} -fastqout ${fname}.merged.fq -fastq_maxdiffs 50 -fastq_pctid 60 -sample $fname
#
# #if you are getting poor merging see: https://www.drive5.com/usearch/manual/merge_report.html
# #if you have huge alignment then increase the fastq_maxdiffs and decrease fastq_pct
#
# 	vsearch --fastq_mergepairs ${fname} --reverse ${fname2} --fastqout ${fname}.merged.fq --fastq_allowmergestagger --fastq_maxns 10 --fastq_minovlen 10 --fastq_maxdiffs 10 --label_suffix $fname
#   #--fastq_maxns means number of allowable NAs
#   #--fastq_minovlen minimum overlap
#   #--fastq_maxdiffs max number of differences in overlap region
#   #--label_suffix add file name to header
#
# done

#######################################
####Step 3. Check reads are oriented the same way
#######################################

#note they should be oriented in same way from MiSeq, but worth doublechecking. This may need to be reconsidered
#if we try other sequencing techniques for which I am unfamiliar

#Note that I am just checking for a few samples, as all samples should be the same
#note the database for comparison may need a new path, or be changed if you have a different db
#
# cat `ls *merged.fq | head -n 5` > testFiles.fq
#
# usearch -orient testFiles.fq -db ~/ref_db/rdp_its_v2.udb -tabbedout orient_out.txt
#
# #this should show most reads are + or ? (the latter being ones that didn't match anything in the db)
#
# echo "Sequence orientation for the first 5 samples." >> out/Processing_Summary.txt
#
# orient=`eval cut -f2 orient_out.txt | sort | uniq -c`
#
# echo "$orient" >> out/Processing_Summary.txt
#
# rm -rf orient_out.txt
# rm -rf testFiles.fq

#######################################
####Step 4. Remove primer binding region
#######################################

# Strip primers (V4F is 18, V4R is 20) (ITS1F is 22, ITS2R is 20)
# this can vary by primer pair
# this step is not neccesary imo when using zotus/esvs/asvs, but I am doing it
#just in case there is a sequence call error here that would lead to a
#spurious OTU

# for f in *merged.fq
# do
# 	usearch -fastx_truncate $f -stripleft 0 -stripright 0 -fastqout ${f}stripped.fq
# done
#
# #Note we don't do any more length trimming because we are using merged reads

#######################################
####Step 5. Sequence filtering
#######################################

#Note that we do this after merging so we can leverage info from both reads when calling bases
#If you are using a sequence technology prone to errors, then best to rethink this step
#or at least look into the output more carefully,
#see https://www.drive5.com/usearch/manual/pipe_readprep_filter.html

for f in parsed*R1*
do
	vsearch --fastx_filter $f --fastq_maxee 1.0 --fastq_stripleft 19 --fastq_stripright 0 --fastaout ${f}.filtered.fa
done

#######################################
####Step 6. Find unique read sequences and abundances
#######################################

cat *filtered.fa > combined_filtered.fa

usearch -fastx_uniques combined_filtered.fa -fastaout uniqueSequences.fa --sizeout

#Count number of original reads
echo "Original number of forward raw reads (unmerged, unfiltered)" >> out/Processing_Summary.txt
echo `grep "@M0" *R1*.fastq | wc -l` >> out/Processing_Summary.txt

# echo "Original number of reverse raw reads (unmerged, unfiltered)" >> out/Processing_Summary.txt
# echo `grep "@M0" *R2*.fastq | wc -l` >> out/Processing_Summary.txt

#number of reads that merged successfully
# echo "Number of reads that merged successfully" >> out/Processing_Summary.txt
# echo `grep "^@" *merged.fq | wc -l` >> out/Processing_Summary.txt

#number of reads that passed filtering
echo "Number of merged reads that passed filtering" >> out/Processing_Summary.txt
echo `grep ">" combined_filtered.fa | wc -l` >> out/Processing_Summary.txt

#number of unique sequences
echo "Number of unique sequences" >> out/Processing_Summary.txt
echo `grep ">" uniqueSequences.fa | wc -l` >> out/Processing_Summary.txt

# #######################################
# ####Step 7. make OTUs and ZOTUs
# #######################################

## # Make 97% OTUs and filter chimeras
#note we remove singletons here

usearch -unoise3 uniqueSequences.fa -zotus zotusF.fa


##################################################
# Step 8. Make OTU table
##################################################
#Note the maxrejects option for the otutab function slows down that function a lot bc it reduces heuristics
#it reduces the probability of errors however.
#if speed becomes an issue for some reason, then reconsider this option.
sed 's/Zo/O/' zotusF.fa > zotusF4b.fa


sed 's/-/_/g' stripped_parsed_forward.fasta > typofixed_stripped_parsed_forward.fasta

sed 's/\./_/g' typofixed_stripped_parsed_forward.fasta > typofixed2_stripped_parsed_forward.fasta

usearch -otutab  typofixed2_stripped_parsed_forward.fasta -zotus zotusF4b.fa -otutabout otutabF.txt -notmatched unmappedF.fa -sizeout




rm -f included*
rm -f otuTitles*
rm -f missing_labelsZ.txt
rm -f missing_labels97.txt
rm -f notmissing_labels97.txt
rm -f notmissing_labelsZ.txt
rm -f notmissing*
mv missingVs* out/
mv unmapped* out/
rm -f *aln


# # ##################################################
# # # Step 9. Clean up files
# # ##################################################

#  #rm -rf *merged*
#
#  gzip *fastq
#  #mkdir rawreads
#
#  #mv *fastq.gz rawreads
#  mv otus97.fa out/
#  mv zotus.fa out/
#
#  rm -rf all_ucs*uc
#  rm -rf uniqueSequences.fa
#  rm -rf combined_filtered.fa
# # tar -czf fastq_info.tar.gz  fastq_info/
#  rm -rf fastq_info/
#  rm -rf tmpfile*
#  rm -rf allstrip.fa
