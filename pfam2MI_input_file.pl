#!/usr/bin/perl -w
open(IN, "HSA-MSA-NEW.fa") or die("no such file\n");
open(OUT, ">hsa-align-new.fa");

my $i = 1;
while(<IN>)
 {
 my($line) = $_;
 chomp($line);
 if($line =~ s/^>//){
     print OUT ">gi|",$i,"|",$line,"\n";
     $i++;}
 else{
     print OUT $line,"\n";}
 }
close(IN);
close(OUT);
