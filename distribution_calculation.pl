#!/usr/bin/perl -w
use strict;
while(@ARGV<2){
	print "\nUSAGE:\n\tdistribution_calculation.pl [input_file] [bin] \n\n";die;
}

my $inputfile=shift @ARGV;
open(INP,$inputfile) or die"NOTE: please input a similarity file, a bin!\n";
my $bin=shift @ARGV;
my $maxNum=0;
my @allnum=();
while(<INP>){
	chomp($_);
	if($_){
		my @TheLine=split(/\s+/,$_);
		my $number=pop(@TheLine);
		push(@allnum,$number);
		if($number>$maxNum){
			$maxNum=$number;
		}
	}
}
print $maxNum,"\n";
my $interval=0;
if(($maxNum)/$bin-int(($maxNum)/$bin)>0){
	$interval=int(($maxNum)/$bin)+1;
}else{
	$interval=int(($maxNum)/$bin);
}
my @InterArray=();
for(my $i=0;$i<$interval;$i++){
	push(@InterArray,0);
}


foreach my $number(@allnum){
	for(my $j=0;$j<$interval;$j++){
		if($number>=$j*$bin and $number<($j+1)*$bin){
			$InterArray[$j]++;
			last;
		}
	}
}


open(OUT,">$inputfile.bin$bin.dist");

for(my $k=0;$k<$interval;$k++){
	print OUT $k*$bin,"\t",$InterArray[$k],"\t",($InterArray[$k])/scalar(@allnum),"\n";
}

