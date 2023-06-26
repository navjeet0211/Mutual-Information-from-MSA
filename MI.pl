#! /usr/bin/perl -w
use strict;

#########
#Get MI values from a single MSA
#########
#

(my $version = <<VER_DOC) =~ s/\t+//gm;
	--------------
	     MI.pl
	--------------
	version 0.2.3
	Oct 331, 2007
	
	Minor code cleanups that do not affect function
	just esthetics
	--------------
	version 0.2.2
	March 29, 2005
	
	Now correctly captures the case where both residues
	have 0 entropy (both are absolutely conserved)
	--------------
	version 0.2.1
	Oct 25, 2004
	
	Added a check for invalid sequence characters in 
	the input.
	
	Added the ability to use GI numbers in addition to taxon IDs.
	
	--------------
	version 0.2
	Sept 17, 2004
	
	Added the reciprocal information to the output
	i.e. added position d info to the output
	so it now gives exactly the same output as
	covary.pl without any real speed or memory penalty
	
	Rounded the outputs to 5 signficant digits
	
	--------------
	version 0.1
	Aug 9-11, 2004
	
	This is the initial write of a program
	to determine MI of a single input sequence.
	Hopefully, I can get this to take arbitrarily
	large input sequences.
	
	The program gives the expected output for the 
	MAP1 alignment.
	
	---end-of-version-info
VER_DOC

(my $help = <<HELP_DOC) =~ s/\t+//gm;
	Arguments for this program are:
	-i first_file (mFasta format)
	-g use GI number or Taxon ID 
	     T - use GI number (default)
	     F - use Taxon ID  
	     NOTE: if this option is chosen only one sequence can 
	     be analyzed. If a second sequence is given, it will be
	     ignored.
	-o freq_output file base name. The tab-delimited freq_output file
	     will be named this with .txt appended. The
	     table of counts will be named this with .table appended.
	     The default name is the reference taxon.
	-t reference taxon ID or GI ID 
	-a optional but recommended! 
	     amino acid residue start point for reference sequence
	     the program tries to read this automatically but is 
	     admittedly not very good at it.
	-h print this help
	-v print version information
	
	Given an aa at position i, what is the probability that it varies with
	any other amino acid at every other position. The output is keyed to a
	reference sequence in the alignment. Only the 20 canonical amino acids
	are accepted.
	
	The tab-delimited freq_output is:
	c d entropy_c entropy_d entropy_cd MI_cd aa1 aa2
	
	Requires as input a fasta formatted file. One can optionally
	put taxonomy information in the second field if using the taxon ID option
	(i.e. gi|gi_num|tx|tx_num|free text description).
HELP_DOC

my $file_name_1 = my $ref_tx = my $outfile = "";
my $GI_or_Tx = "T";
my $MSA_start = 0;
my $base = 21;	#logarithm base
my @MI_pos = my @Ref_pos;
my $ID_column = 3;

######GET THE PARAMETERS######
#step through ARGV and set the parameters
if (@ARGV){
	my $i = 0;
	foreach my $item(@ARGV){
		$i++;
		if ($item eq "-i"){
			$file_name_1 = $ARGV[$i];
		}elsif ($item eq "-t"){
			$ref_tx = $ARGV[$i];
		}elsif ($item eq "-o"){
			$outfile = $ARGV[$i];
		}elsif ($item eq "-g"){
			$GI_or_Tx = $ARGV[$i];
		}elsif ($item eq "-v"){
			print $version;
			exit;
		}elsif ($item eq "-a"){
			$MSA_start = $ARGV[$i];
		}elsif ($item eq "-h"){
			print $help;
			exit;
		}
	}
}else{
	print $help;
	exit;
}

if ($outfile eq ""){ #use the taxon name as the default file name
	$outfile = $ref_tx;
}

if ($GI_or_Tx  eq "F"){	#use Taxon as ID
	$ID_column = 3;
	init_open_TX($file_name_1);
}else{
	#print "Sorry GI sorting is not yet implemented. Please make\nsure that there is a unique ID number in the second information field\n";
	$ID_column = 1;
	init_open_TX($file_name_1);
}



############SUBROUTINES###########

####Open file routine that converts fasta formatted MSA to single lines
####keyed by taxon ID number
####Get the counts and co-counts

sub init_open_TX{
	
	#declare variables
	my $last_ID = "";	
	my @line;
	my @temp;
	my %UID_ID;
	
	open(INPUT, "< $_[0]") || die "Could not open $_[0] $!\n";
	my $n = 0;
	#read each line and remove the trailing return 
	while (defined (my $line = <INPUT>)){	
		chomp($line);
		
		#split on | if it is a definition line, and 
		#create a new empty taxon key to hold the sequence info	
		#and set the key value to the taxon retreived from the split line
		#otherwise add the sequence information to the hash, keyed by the last ID
		if ($line =~ m/^>/) {	
			@line = split(/\|/, $line);
			
				$UID_ID{$line[$ID_column]} = "";	#set up a blank to hold the sequence
				$last_ID = $line[$ID_column];	#keep the ID number
				
				#if we are dealing with the reference taxon  
				#grab the start point of the sequence
				#assumes sequence numbering is of the form "/5-134" and
				#that this is the last piece of information in the definition line
				#only do this if the MSA start point is not set in the command line
				if ($MSA_start == 0){
					if($line[$ID_column] eq $ref_tx){
						@temp = split(/-|\//,$line[-1]);
						$MSA_start = $temp[1]; 
					}
				}
		}else{	#this hash holds the sequence keyed by tx ID 
			if ($line !~ m/[XxBbJjOoUuZz]/){
				$UID_ID{$last_ID} .= $line;
			}else{
				print "The sequence $last_ID contains an invalid character at line $n\n";
				exit;
			}
		}
		
		$n++;
	}
	close INPUT;
	
	#print "Got the sequences!\n";
	#now process the reference sequence, getting only the positions that do
	#not contain gaps
	
	my @ref_line = split(//,$UID_ID{$ref_tx});
	
	#this generates two arrays with numbers that 
	my $i = 0;
	my $j = $MSA_start;
	foreach(@ref_line){
		if($_ ne "-"){	#if it is not a gap character
			push @MI_pos, $i;
			push @Ref_pos, $j++;
			#$j++;
		}
		$i++;
	}
	my %count_hash;
	my %co_count_hash;
	
	#print "got the ungapped reference sequence\n";
	#now get the residue counts and co-counts
 	my @key = keys %UID_ID; #gets the sequence IDs
 	my $num_seq = 0;
 	my $aa1aa2;
 	foreach(@key){	#run thru each sequence in the %UID_ID hash
 		$num_seq++;
		@line = split(//,$UID_ID{$_});
		#The outer loop gets the counts and places them in a hash keyed by
		#postion number in the MSA and then by amino acid name
		#print "getting counts and co-counts for sequence $_\n";
		for (my $pos = 0; $pos < @MI_pos; $pos++){ #check only ungapped postions in the refseq
			if (exists $count_hash{$MI_pos[$pos]}{$line[$MI_pos[$pos]]}){	#initialize if not already
				$count_hash{$MI_pos[$pos]}{$line[$MI_pos[$pos]]}++;
			}else{
				$count_hash{$MI_pos[$pos]}{$line[$MI_pos[$pos]]} = 1;	#otherwise increment
			}
			#This inner loop gets the co-counts and places them in a hash keyed by
			#position1,position2,aa1aa2 (concatene makes the hash smaller but still
			#is a unique name). 
			for (my $pos2 = $pos + 1; $pos2 < @MI_pos; $pos2++){
				$aa1aa2 = "$line[$MI_pos[$pos]]$line[$MI_pos[$pos2]]";
				if (exists $co_count_hash{$MI_pos[$pos]}{$MI_pos[$pos2]}{$aa1aa2}){
					$co_count_hash{$MI_pos[$pos]}{$MI_pos[$pos2]}{$aa1aa2}++;	#increment
				}else{
					$co_count_hash{$MI_pos[$pos]}{$MI_pos[$pos2]}{$aa1aa2} = 1;	#initialize
				}
			}
		}
	}
	#print "got the counts and co-counts\n";
	@key = sort {$a <=> $b}(keys(%count_hash));
	my @amino_acids = qw (A C D E F G H I K L M N P Q R S T V W Y -);
	my $freq;
	$i = 0;
	my $entropy;
	my @pos_entropy;
	
	#initiate the header line for the frequency output file
	my $freq_output = "pos\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\t-\tentropy\n";
	
	#get the counts and entropies for each position
	foreach my $key (@key){	#check each ungaped position in the msa
					#this form is easier to read than the one where I use
					#positions in @MI_pos
		$entropy = 0;
		$freq_output .= "$Ref_pos[$i]\t";
		foreach my $residue (@amino_acids){	#run thru the residues
			if (!exists $count_hash{$key}{'X'}){	#only ungapped positions
				if (exists $count_hash{$key}{$residue}){	#not all residues at all positions
					$freq = $count_hash{$key}{$residue}/$num_seq;
					$entropy += -1 * $freq * (log($freq)/log($base));
					$freq_output .= "$count_hash{$key}{$residue}\t";
				}else{
					$freq_output .= "0\t";
				}
			}else{
				$freq_output .= "gap";
			}
		}
		#create a hash of entropies keyed to the positions in the MSA
		push @pos_entropy, $entropy;
		$freq_output .= "$entropy\n";
		$i++;
	}
	my $countoutfile = $outfile . "_count.txt";
	open (OUT, "> $countoutfile") || die "cannot open  $countoutfile\n";
		print OUT $freq_output;
	close OUT;
	#print $freq_output;
	#foreach (@pos_entropy){print "$_\t";}
	#works
	#print $count_hash{0}{'-'} . "<- count of postion 1, - , number of seqs: $num_seq\n";
	
	#now do the heavy lifting
	#calculate the co-count frequency, joint entropy and MI
	#save it into a string for printing later
	my $joint_entropy = 0;
	my $pair;
	my $Hc = my $Hd = my $Hcd = my $MI;
	my $MI_out = "c\td\tH_c\tH_d\tH_cd\tMI_cd\taa_c\taa_d\n";
	for (my $pos = 0; $pos < @MI_pos; $pos++){ #check only ungapped postions in the refseq
		if (!exists $count_hash{$MI_pos[$pos]}{'X'}){	#check only ungapped alignment positions
			for (my $pos2 = $pos + 1; $pos2 < @MI_pos; $pos2++){
				$joint_entropy = $MI = 0;
				if (!exists $count_hash{$MI_pos[$pos2]}{'X'}){	#check only ungapped alignment positions
					foreach $pair (keys%{$co_count_hash{$MI_pos[$pos]}{$MI_pos[$pos2]}}){	#list of co-counts
						$freq = $co_count_hash{$MI_pos[$pos]}{$MI_pos[$pos2]}{$pair} / $num_seq; #get the co-count frequencies
						$joint_entropy += -1 * $freq * (log($freq) / log($base));	#sum them
					}
				}
				if ($joint_entropy !=0){
					$MI = sprintf("%.5f", ($pos_entropy[$pos] + $pos_entropy[$pos2] - $joint_entropy));
		#### NEW ####
		#added the information for the converse to the output. 
		#this makes the output congruent with covary.pl and useable by the
		#downstream programs without modification.
					$Hc = sprintf("%.5f", $pos_entropy[$pos]);
					$Hd = sprintf("%.5f", $pos_entropy[$pos2]);
					$Hcd = sprintf("%.5f", $joint_entropy);
	
					$MI_out .= "$Ref_pos[$pos]\t$Ref_pos[$pos2]\t$Hc\t$Hd\t$Hcd\t$MI\t$ref_line[$MI_pos[$pos]]\t$ref_line[$MI_pos[$pos2]]\n";
					$MI_out .= "$Ref_pos[$pos2]\t$Ref_pos[$pos]\t$Hd\t$Hc\t$Hcd\t$MI\t$ref_line[$MI_pos[$pos2]]\t$ref_line[$MI_pos[$pos]]\n";
				}elsif($joint_entropy ==0){
		#### NEW ####
		#added the information when both positions are absolutely conserved 
					$Hc = 0;
					$Hd = 0;
					$Hcd = 0;
	
					$MI_out .= "$Ref_pos[$pos]\t$Ref_pos[$pos2]\t$Hc\t$Hd\t$Hcd\t$MI\t$ref_line[$MI_pos[$pos]]\t$ref_line[$MI_pos[$pos2]]\n";
					$MI_out .= "$Ref_pos[$pos2]\t$Ref_pos[$pos]\t$Hd\t$Hc\t$Hcd\t$MI\t$ref_line[$MI_pos[$pos2]]\t$ref_line[$MI_pos[$pos]]\n";
				
				}
			}
		} 
	}	
	my $mioutfile = $outfile . "_MI.txt";
	open (OUT, ">$mioutfile") || die "cannot open $mioutfile\n";
		print OUT $MI_out;
	close OUT;
}
