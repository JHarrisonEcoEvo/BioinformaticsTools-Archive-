#!/usr/bin/perl
use strict; use warnings;


#Harrison perl 4
my ($ucs) = shift(@ARGV);

#input and output file
open(IN, "< $ucs") or die;
open(OUT, "> all_matches.uc") or die "couldnt open OUT\n";


my $pattern = '*';
	
while (my $line=<IN>){
	#chomp $line;
	my @reads = split(" ", $line);
	#my $xx= $reads[9];
		 if ($reads[9] eq $pattern){next}else{print OUT "@reads", "\n";}
	}

close (IN);
close (OUT);
