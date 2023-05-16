#!/usr/bin/perl;
use warnings;
use strict;
open IN,"$ARGV[0]" or die$!;
my %hash;
while (<IN>){
	chomp;
	my @line=split/\t/,$_;
	my $M1=join"\t",$line[0],$line[1];
	my $M2=join"\t",$line[1],$line[0];
	$hash{$M1}=$line[4];
	$hash{$M2}=$line[4];
	};
close IN;
open IN2,"$ARGV[1]" or die$!;
open OUT,"> $ARGV[2]" or die$!;
print OUT"window_bin";
my @arr;
while (<IN2>){
	if($_=~/^L1/){}
	else{	chomp;
		print OUT"\t$_";
		push @arr,$_;
	}};
print OUT"\n";
close IN2;
open IN3,"$ARGV[1]" or die$!;
while (<IN3>){
	if($_=~/^L1/){}
	else{	chomp;
		print OUT"$_";
		my $name=$_;
		foreach(0..$#arr-1){
			if($arr[$_] eq $name){
				print OUT"\t1"}
			else{	my $Ms=join"\t",$name,$arr[$_];
				print OUT"\t$hash{$Ms}"}
		};
		if($arr[-1] eq $name){
			print OUT"\t1\n"}
		else{	my $Ms=join"\t",$name,$arr[-1];
			print OUT"\t$hash{$Ms}\n"};
		};
	};
