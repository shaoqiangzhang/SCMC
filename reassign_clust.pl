#!/usr/bin/perl
use warnings;
use strict;
#use List::Util qw/sum max min/;
#-----------------------------------------------------------
##By Shaoqiang Zhang(zhangshaoqiang@tjnu.edu.cn)
##July-20-2021
#-----------------------------------------------------------
while(@ARGV<1){
        print "\nUSAGE:\n\treassign_clust.pl <mcl_file> <sim_file> <min_size>\n\n";
        die;
}

my $mclfile=shift @ARGV;
my $simfile=shift @ARGV;
my $minsize=shift @ARGV;

open(INPF,$mclfile) or die "ERROR: please input a clustering result file!!\n";
open(CLU,">$mclfile.renew");

my %cell2clust_hash=();
my @cellClustArray=();
my @aray_of_small_clusters=();
my %rarecellshash=();
my $count=0;
while(<INPF>){
	chomp $_;
	my $theline=$_;
	if($theline){
		$count++;
		@cellClustArray=split(/\s+/,$_);
		if(scalar(@cellClustArray)>$minsize){
			foreach my $cell (@cellClustArray){
				$cell2clust_hash{$cell}=$count; ##save cells to large cluster labels
			}
		}else{
			push(@aray_of_small_clusters, $theline);
			foreach my $cell (@cellClustArray){
				$rarecellshash{$cell}=-1;
			}
		}
	}
}##save rare cell cluster labels to hash

open(SIM,$simfile) or die "ERROR: please input a cell similarity file";
my %bestNeighborHash=();
my %bestneighborsim=();
my @cellsimarray=(); 
while(<SIM>){
	chomp $_;
	my @cellsimarray=split(/\t+/,$_);
	if(exists($rarecellshash{$cellsimarray[0]}) and not(exists($rarecellshash{$cellsimarray[1]}))){
		unless(exists($bestNeighborHash{$cellsimarray[0]})){
			$bestNeighborHash{$cellsimarray[0]}=$cellsimarray[1];
			$bestneighborsim{$cellsimarray[0]}=$cellsimarray[2];
		}elsif($bestneighborsim{$cellsimarray[0]}<$cellsimarray[2]){
			$bestNeighborHash{$cellsimarray[0]}=$cellsimarray[1];
			$bestneighborsim{$cellsimarray[0]}=$cellsimarray[2];
		}
	}elsif(exists($rarecellshash{$cellsimarray[1]}) and not(exists($rarecellshash{$cellsimarray[0]}))){
		unless(exists($bestNeighborHash{$cellsimarray[1]})){
			$bestNeighborHash{$cellsimarray[1]}=$cellsimarray[0];
			$bestneighborsim{$cellsimarray[1]}=$cellsimarray[2];
		}elsif($bestneighborsim{$cellsimarray[1]}<$cellsimarray[2]){
			$bestNeighborHash{$cellsimarray[1]}=$cellsimarray[0];
			$bestneighborsim{$cellsimarray[1]}=$cellsimarray[2];
		}
	}
}



foreach my $oneline(@aray_of_small_clusters){
	my @cellarray=split(/\s+/,$oneline);
	my $maxcellsim=0;
	my $maxcell;
	foreach my $cell (@cellarray){
		if($bestneighborsim{$cell}>$maxcellsim){
			$maxcell=$cell;
			$maxcellsim=$bestneighborsim{$cell};
		}
	}
	foreach my $cell (@cellarray){
		$cell2clust_hash{$cell}=$cell2clust_hash{$bestNeighborHash{$maxcell}};
	}##set cells of the small cluster to a big cluster with max 
}

foreach my $key  ( sort { $a <=> $b } keys %cell2clust_hash ) { 
	print CLU $cell2clust_hash{$key},"\n";
}



