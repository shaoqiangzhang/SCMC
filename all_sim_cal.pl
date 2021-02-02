#!/usr/bin/perl -w
use strict;
while(@ARGV<2){
	print "\nUSAGE:\n\tall_cell_sim_calculation.pl <scRNA-seq_expression_file> <cell_number> \n\n";
	die;
}
my $inpfile=shift @ARGV;
my $ncells=shift @ARGV;

for(my $i=0;$i<($ncells-1);$i++){
	system("./cal_sim_cuda $inpfile $i  >>$inpfile.sim");	
}
