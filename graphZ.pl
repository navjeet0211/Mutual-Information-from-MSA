#! /usr/bin/perl -w
use strict;

########
#This program outputs a set of nodes in a format
#that can be read by Graphviz to draw a connectivity
#graph between MI-containing residues
########

####declare the variables
my $file_name_1;
my $min_H = 0.3;
my $outfile;
my $Z_range;
my $version;
my $min_Z = 4.5;
my $mid_Z = 5;
my $max_Z = 6;
my @sorted_array;
my $output = "graph A {\n\tgraph [nodesep=\"0.10\", size=\"12,12\", bb=\"0,0,737,832\"];\n\tnode [fontsize=\"12.00\", shape=circle, height = 0.4 width = 0.6 fixedsize=\"true\"];\n";
my $Zcol = 12; 

my $parse = "T";


#this type of formating sets the variable to a multiline format. The
#substitution removes any preceeding tabs. Indents are set with spaces
(my $help = <<HELP_DOC) =~ s/\t+//gm;
	Arguments for this program are:
	-i first_file (peerZ output format)
	-m minimum H value; default 0.3
	-o output file base name. The tab-delimited output file
	     will be named this with .dot appended.
	-r list of numerical ranges; default "6 5 4.5"
	     the first value will be drawn as a solid line
	     the second value will be drawn as a dashed line
	     the third value will be drawn as a dotted line
	     lines will be labeled with their respective Z score
	-h print this help
	-v print version information
HELP_DOC
######GET THE PARAMETERS######
#step through ARGV and set the parameters
if (@ARGV){
	my $i = 0;
	foreach my $item(@ARGV){
		$i++;
		if ($item eq "-i"){
			$file_name_1 = $ARGV[$i];
		}elsif ($item eq "-m"){
			$min_H = $ARGV[$i];
		}elsif ($item eq "-o"){
			$outfile = $ARGV[$i];
		}elsif ($item eq "-r"){
			$Z_range = $ARGV[$i];
		}elsif ($item eq "-v"){
			print $version;
			exit;
		}elsif ($item eq "-h"){
			print $help;
			exit;
		}
	}
}else{
	print $help;
	exit;
}
$parse = "F" if $Zcol == 9;
my @unsort_Z_range = split/\s+/, $Z_range if $Z_range;
if (@unsort_Z_range == 3){
	my @sorted_Z_range = sort  {$b <=> $a} @unsort_Z_range;
	$max_Z = $sorted_Z_range[0];
	$mid_Z = $sorted_Z_range[1];
	$min_Z = $sorted_Z_range[2];
}

#open the file and generate a sorted list of outputs
init_open();

#take the sorted list and generate the graph
#split each value of the sorted array and use the last entry as the Z score
#add style = "dashed_dotted"]\n; to each array element and save as a string
init_make_graph();

########SUBROUTINES#####
sub init_make_graph{
	my @line;
	my $value;
	foreach (@sorted_array){
		@line = split/\s/, $_;
		if ($line[-1] >= $max_Z){
			$output .= "$_ style = \"bold\" ];\n";
		}elsif ($line[-1] >= $mid_Z && $line[-1] < $max_Z){
			$output .= "$_ style = \"solid\"];\n";
		}else{
			$output .= "$_ style = \"dashed\"];\n";
		}
	}	#end foreach
	$output .= "}\n";
	my $graphoutfile = $outfile . "_MIp.dot"; 
	open (OUT, "> $graphoutfile") || die;
		print OUT $output;
	close OUT;
}	#end sub init_make_graph



sub init_open{
	my %found_hash;
	my @idx;
	my @unordered_array;
	my $value;
	open(INPUT, "< $file_name_1") || die "Could not open $file_name_1 $!\n";
	#read each line and remove the trailing return 
	while (defined (my $line = <INPUT>)){	
		chomp($line);
		if ($line =~ m/^\d/) {	#split on the tabs if it is a data line
			my @line = split(/\t/, $line);
			#check to see if both positions have been found together before
			#keep a hash of found positions and check if both exist already
			#if the hash is empty do the next positions
			if (!exists $found_hash{$line[0]}){
				#capture the two positions and the Z score if the 
				#entropy of both is > the minimum entropy and 
				#the Z score is greater than the lowest value in the
				#Z score range
				if ($line[2] >= $min_H && $line[3] >= $min_H && $line[$Zcol] >= $min_Z){
					#generate the Graphviz code for the line
					#keep in an unordered list indexed by Z score
					$value = sprintf ("%.1f", $line[$Zcol]);
					push @idx, $value;
					push @unordered_array, "$line[6]$line[0] -- $line[7]$line[1] [label = $value";
					$found_hash{$line[0]} = 0;
					$found_hash{$line[1]} = 0;
				} #end if
			}else{	#if the hash is full delete the values
				delete($found_hash{$line[0]});
				delete($found_hash{$line[1]})
			}
		}	#end if
	} #end while
	close INPUT;
	#now sort based on the data field saved in the index each line
	@sorted_array = @unordered_array[ sort { $idx[$b]<=>$idx[$a] } 0 .. $#idx ]; 
}	#end sub init_open
