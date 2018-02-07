# BioinformaticsTools

A handful of simple scripts that help me automate common tasks.

Hey internet! Beware! These are mostly made to fit my own workflow, or to share
with collaborators. Feel free to use, but don't use blindly!


Explanation of what each script does:

OTUtable_output_rarified_TMM.R - takes an OTU table and outputs a rarified
version and a TMM normalized version. Can tweak to change rarefaction level
or normalization techniques

OTUtable_to_adjMatrix.R - makes an adjacency matrix (incidence based) out of an
OTU table. There is a script that makes a weighted adj. matrix too.

OTUtable_to_edgeList.R - get an edge list from an OTU table.

rmv_cpDNAmtDNA_from_OTUtable.R - removes possible non target OTUs.
This script will undoubtedly need tweaked to fit somebody else's needs.

rmv_matches.pl - remove lines from uc output that correspond to matches to a
 taxonomic databases

uc_to_OTUtable.R - take a big uc file and turn it into an OTU table.
Not really any reason to do this now that usearch has this function.

useTax_to_cleanOTUtable.R
