#!usr/bin/perl -w
use strict;
use SVG;


my $svg = SVG->new('width', 3000, 'height', 1800);
my $win;
#$svg->rect('x', 100, 'y', 100, 'width', $len, 'height', 5, 'stroke', "none", 'fill', "black");
my $name =$ARGV[1];
open IN,$ARGV[0] or die $!;
my $start;my $len;
while (<IN>){
	chomp;
	next if($_=~/^#/);
	my @a = split /\t/,$_;
	if ($a[2] eq "mRNA"){
		if ($a[8]=~/ID=(.+?);/){
			my $id = $1;
			if ($id eq $name){
				$start = $a[3];
				$len = $a[4] - $a[3] + 1;
			}
		}
	}
}
$win = 2400/$len;
close IN;
$svg->rect('x', 100, 'y', 100, 'width', $len*$win, 'height', 5, 'stroke', "none", 'fill', "black");
$svg->rect('x', 100, 'y', 150, 'width', $len*$win, 'height', 3, 'stroke', "none", 'fill', "black");

open IN1,$ARGV[0] or die $!;
while (<IN1>){
	chomp;
	next if($_=~/^#/);
        my @a = split /\t/,$_;
	if ($a[2] eq "mRNA"){
		if ($a[8]=~/ID=(.+?);/){
		my $id2 = $1;
		if ($id2 eq $name){
		if ($a[6] eq "+"){
                        $svg->text('x', 90 , 'y', 140, '-cdata', "5'", 'text-anchor', 'middle', 'font-size', 30);
                        $svg->text('x', 90 + $len*$win + 20 , 'y', 140, '-cdata', "3'", 'text-anchor', 'middle', 'font-size', 30);
                }else{
                        $svg->text('x', 90 , 'y', 140, '-cdata', "3'", 'text-anchor', 'middle', 'font-size', 30);
                        $svg->text('x', 90 + $len*$win + 20 , 'y', 140, '-cdata', "5'", 'text-anchor', 'middle', 'font-size', 30);
                }	
		}
		}
	}
	if($a[2] eq "CDS"){
        if ($a[8]=~/Parent=(.+);/){
        my $id = $1;
        if ($id eq $name){
        my $s = $a[3] - $start;
        $svg->rect('x', 100 + $s*$win, 'y', 100 - 15,  'width', ($a[4] - $a[3] + 1)*$win , 'height',  30, 'stroke', "black", 'fill', "#DAA520");
	}
        }elsif($a[8]=~/Parent=(.+)/){
	my $id = $1;
        if ($id eq $name){
        my $s = $a[3] - $start;
        $svg->rect('x', 100 + $s*$win, 'y', 100 - 15,  'width', ($a[4] - $a[3] + 1)*$win , 'height',  30, 'stroke', "black", 'fill', "#DAA520");
        }
	}
	}elsif($a[2] =~/prime_UTR/){
		if ($a[8]=~/ID=(.+?)\.utr/){
                        my $name3 = $1;
                        if ($name3 eq $name){
                                $svg->rect('x', 200 + ($a[3] - $start)*$win, 'y', 85,  'width', ($a[4] - $a[3])*$win , 'height',  30, 'stroke', "black", 'fill', "#bebebe");
                        }
                }
	}
}
close IN1;

my $n = int ($len/500);

for (my $i = 0;$i <= $n;$i++){
	my $table = $i*500;
	$svg->rect('x', 100 + ($i*500)*$win, 'y', 145, 'width', 3, 'height', 8, 'stroke', "none", 'fill', "black");
	$svg->text('x', 100 + ($i*500)*$win, 'y', 185, '-cdata', "$table bp", 'text-anchor', 'middle', 'font-size', 30);
}


$svg->rect('x', 300 , 'y',40,   'width', 100 , 'height',  30, 'stroke', "black", 'fill', "#DAA520");
$svg->rect('x', 100 , 'y', 40,  'width', 100 , 'height',  30, 'stroke', "black", 'fill', "#A020F0");
$svg->rect('x', 500, 'y', 50, 'width', 100, 'height', 5, 'stroke', "none", 'fill', "black");


$svg->text('x', 250 , 'y', 65, '-cdata', "UTR", 'text-anchor', 'middle', 'font-size', 30);
$svg->text('x', 450 , 'y', 65, '-cdata', "CDS", 'text-anchor', 'middle', 'font-size', 30);
$svg->text('x', 650 , 'y', 65, '-cdata', "Intron", 'text-anchor', 'middle', 'font-size', 30);

##############draw SNP##############

if ($ARGV[2]=~/\.gz/){
open IN3,"<:gzip",$ARGV[2] or die $!;
}else{
open IN3,$ARGV[2] or die $!;
}
my $e = 0;
while (<IN3>){
	chomp;
	my @a= split /\t/,$_;
	if ($_=~/#/){
		next;
	}else{
		my $pos = $a[1] - $start;
		$e++;			
	}
}
close IN3;

my $w = 2200/$e;
if ($ARGV[2]=~/\.gz/){
open IN4,"<:gzip",$ARGV[2] or die $!;
}else{
open IN4,$ARGV[2] or die $!;
}
my $m = 0;
while (<IN4>){
        chomp;
        my @a= split /\t/,$_;
        if ($_=~/#/){
                next;
        }else{
                my $pos = $a[1] - $start;
		$svg->line('x1', 100+$pos*$win , 'y1',150,  'x2', 200+$m*$w+$w/2 , 'y2',300, 'stroke', "black", 'stroke-dasharray', "1 2");
        	$m++;
	}
}
close IN4;


##############draw SNP_type##############
open IN5,$ARGV[3] or die $!;
my %hash;
my $sample_num = 0;
$/=">";<IN5>;
while (<IN5>){
        chomp;
        my ($id,$seq) = split /\n/,$_,2;
        $seq=~s/\n//g;
        $hash{$id}=$seq;
	$sample_num++;
}
close IN5;

my $win2=1000/$sample_num;
open IN6,$ARGV[4] or die $!;
$/="\n";
my $nn = 0;
while (<IN6>){
        chomp;
        $nn+=$win2;
        my @a = split /\t/,$_;
        my @b= split //,$hash{$a[0]};
        my $u=0;
        foreach my $i(@b){
                if ($i eq "A"){
                        $svg->rect('x', 200 + $u*$w, 'y', 300 +$nn, 'width', $w , 'height',  $win2, 'stroke', "none", 'fill', "#6AAB9C");
                }elsif($i eq "T"){
                        $svg->rect('x', 200 + $u*$w, 'y', 300 +$nn, 'width', $w , 'height',  $win2, 'stroke', "none", 'fill', "#FA9284");
                }elsif($i eq "C"){
                        $svg->rect('x', 200 + $u*$w, 'y', 300 +$nn, 'width', $w , 'height',  $win2, 'stroke', "none", 'fill', "#E06C78");
                }elsif($i eq "G"){
                        $svg->rect('x', 200 + $u*$w, 'y', 300 +$nn, 'width', $w , 'height',  $win2, 'stroke', "none", 'fill', "#5874DC");
                }elsif($i eq "-"){
                        $svg->rect('x', 200 + $u*$w, 'y', 300 +$nn, 'width', $w , 'height',  $win2, 'stroke', "none", 'fill', "#384E78");
                }elsif($i eq "H"){
                        $svg->rect('x', 200 + $u*$w, 'y', 300 +$nn, 'width', $w , 'height',  $win2, 'stroke', "none", 'fill', "#ababa9");
                }
		$u++;
        }
        #$svg->text('x', 50 , 'y', 300 +$nn, '-cdata', "$a[0]", 'text-anchor', 'middle', 'font-size', 20);
}
close IN6;

$svg->rect('x', 300 , 'y',1350,   'width', 50 , 'height',  30, 'stroke', "black", 'fill', "#6AAB9C");
$svg->rect('x', 400 , 'y',1350,   'width', 50 , 'height',  30, 'stroke', "black", 'fill', "#FA9284");
$svg->rect('x', 500 , 'y',1350,   'width', 50 , 'height',  30, 'stroke', "black", 'fill', "#E06C78");
$svg->rect('x', 600 , 'y',1350,   'width', 50 , 'height',  30, 'stroke', "black", 'fill', "#5874DC");
$svg->rect('x', 700 , 'y',1350,   'width', 50 , 'height',  30, 'stroke', "black", 'fill', "#ababa9");
$svg->rect('x', 800 , 'y',1350,   'width', 50 , 'height',  30, 'stroke', "black", 'fill', "#384E78");


$svg->text('x', 370 , 'y', 1375, '-cdata', "A", 'text-anchor', 'middle', 'font-size', 30);
$svg->text('x', 470 , 'y', 1375, '-cdata', "T", 'text-anchor', 'middle', 'font-size', 30);
$svg->text('x', 570 , 'y', 1375, '-cdata', "C", 'text-anchor', 'middle', 'font-size', 30);
$svg->text('x', 670 , 'y', 1375, '-cdata', "G", 'text-anchor', 'middle', 'font-size', 30);
$svg->text('x', 770 , 'y', 1375, '-cdata', "H", 'text-anchor', 'middle', 'font-size', 30);
$svg->text('x', 870 , 'y', 1375, '-cdata', "-", 'text-anchor', 'middle', 'font-size', 30);



print $svg->xmlify;

