#!/usr/bin/perl -w

while(@ARGV<2){
	print "\nUSAGE:\n\tselect_sim_cutoff.pl <sim_file>  <sim_score_cutoff>\n\n";
	die;
}
my $inpf=shift;
my $threshold=shift;
open(INP,$inpf);
open(OUT, ">$inpf.cutoff");
my $a=0;
my @array;
while(<INP>){
	chomp($_);
	@array=split(/\s+/,$_);
	$a=pop(@array);
	if($a>=$threshold){
		print OUT "$_\n";
	}	
}
