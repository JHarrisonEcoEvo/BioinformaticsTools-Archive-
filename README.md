# BioinformaticsTools

A handful of simple scripts that help me automate common tasks.

Hey internet! Beware! These are mostly made to fit my own workflow, or to share
with collaborators. Feel free to use, but don't use blindly!

Suggested order of operations:

If you need to demultiplex your reads use the demux.py script. Note that this assumes dual-indexed paired end reads, with the barcode immediately preceeding the primer region. The script uses Levenshtein distance to correct for one error in the barcode. If your data are different, then the script should be pretty easy to modify. 

1. Run fastq_to_OTUtable to get OTU tables and various QC info.
2. Run your consensus OTU fasta files through SINTAX. I think it is best to use a
couple databases for fungi (UNITE and Warcup) as you will get very different
results depending on database.
3. Merge your taxonomy info using the mergeTaxonomyFiles script.
4. Run the check16sfor_mt_cpDNA script if you are dealing with 16s. This will
give you an output of OTUs that are probably mitochondrial or plastid dna.
5. Run the cleanOTUtable script. Carefully consider the option for
cutting taxa that were not placed in a phylum if you are dealing with the ITS.
This is because a lot of times I found that OTUs that were not well placed in a
phylum turned out to be plant. This may not be needed depending upon your training databases. 
6. You now should have OTU tables without plastid or mtDNA. You can then pass
this to further scripts to either model the data, or use the rarify and normalize
script. Note that I am no longer a fan of the rarefaction/normalization approach
for differential abundance testing.

Explanation of what each script does:

demux.py - demultiplexes dual-indexed paired end reads, with the barcode immediately preceeding the primer region. The script uses Levenshtein distance to correct for one error in the barcode.

fastq_to_OTUtable – takes input fastqs and spits out a 97% OTU table and an
exact sequence variant table. It also gives you a lot of info on QC and them
number of reads lost at each step. There is a lot to consider in this script
so best to read it carefully if you use it.

mergeTaxonomyFiles - merges taxonomy hypotheses output from multiple databases.

OTUtable_output_rarified_TMM.R - takes an OTU table and outputs a rarified
version and a TMM normalized version. Can tweak to change rarefaction level
or normalization techniques

OTUtable_to_adjMatrix.R - makes an adjacency matrix (incidence based) out of an
OTU table. There is a script that makes a weighted adj. matrix too.

OTUtable_to_edgeList.R - get an edge list from an OTU table.

rmv_matches.pl - remove lines from uc output that correspond to matches to a
 taxonomic databases

uc_to_OTUtable.R - take a big uc file and turn it into an OTU table.
Not really any reason to do this now that usearch has this function.

Example_data/ - has some sample data files that may help you try out a script or
tweak the script to your own purposes
