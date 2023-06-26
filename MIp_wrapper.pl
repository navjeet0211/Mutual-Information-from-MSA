#!/usr/bin/perl -w
#call MI.pl, peerZ.pl and graph.pl
use strict;
(my $help = <<HELP_DOC) =~ s/\t+//gm;
	
	This generates Z scores with whole alignment neighbors.
	
	Generates a connectivity graph that can be viewed in 
	Graphviz (http://www.graphviz.org/) an open source
	graph visualization program
	
	The peerZ.pl program called by this requires that the 
	Statistics::Lite Perl module be installed
	
	To enable these programs system wide you must modify the source code
	as indicated below on lines 77-79
		
	Arguments for this program are:
	-i input file name (fasta file format only!)
	     required
	-o output file base name
	     required
	     base name will generate three files with the following
	     appended:
	         _count.txt
	         _MIp.txt 
	         _MIp.dot
	-n residue start position in the reference sequence
	     required
	-g gi number of the reference sequence
	     required
	-e minimum entropy value for the calculations (default 0.0001)
	     This is to stratify the data set. In general it is best
	     to leave this at the default.
	-s sum or product (default sum)
	-a return all positions (default F)
	      returns only those above the entropy cutoff
	-h print this help
	     
HELP_DOC

my $file_name = my $min_entropy = my $SorP = my $out_file = my $firstres = my $gi = "";
$min_entropy = 0.0001;
#get the arguements for the calculations
if (@ARGV){
	my $i = 0;
	foreach my $item(@ARGV){
	$i++;
		if ($item eq "-i"){	#input file name
			$file_name = $ARGV[$i];
		}elsif ($item eq "-e"){	#minimum entropy
			$min_entropy = $ARGV[$i];
		}elsif ($item eq "-s"){ #sum or product
			$SorP = $ARGV[$i];
		}elsif ($item eq "-n"){ #residue start position
			$firstres = $ARGV[$i];
		}elsif ($item eq "-g"){ #gi number of reference sequence
			$gi = $ARGV[$i];
		}elsif ($item eq "-o"){ #output file base name
			$out_file = $ARGV[$i];
		}elsif($item eq "-h"){	#print the help information
			print $help;
			exit;
		}
	}
}else{
	print $help;
	exit;
}
if ($firstres eq "" || $file_name eq "" || $out_file eq ""){
	print "#########\nPlease enter all four required values!\n#########\n";
	print $help;
	exit;
}

####
#To enable these system wide:
#copy to a directory in your binary path
#remove the ./ before the program names on the three lines noted below

#remove the ./ here
system ("./MI.pl -i $file_name -o $out_file -a $firstres -t $gi ");
my $Zin = $out_file . "_MI.txt";
my $dotin = $out_file  . "_MIp.txt";

#remove the ./ here
system ("./peerZ.pl -i $Zin -o $out_file ");
#remove the ./ here
system ("./graphZ.pl -i $dotin -o $out_file");
system ("rm $Zin");