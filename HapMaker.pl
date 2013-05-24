#!/usr/bin/perl

# HapMaker.pl
# 23 May 2013 Nozomu Okuda

use strict;
# use warnings;
use Bio::SeqIO;
use Bio::Seq;

my $transWt = 0.5;
my $indelVsnp = 200/7200;
my $indelLength = 99;

usage() if(@ARGV != 5);
my $infile = shift;
my $in = Bio::SeqIO->new(-file => $infile, -format => 'fasta')->next_seq()->seq();
open(my $noChangeFh, '<'.shift) or die "Could not open no-change file: $!";
my $divRate = shift;
my $outname = shift;
my $out = Bio::SeqIO->new(-file => '>'.$outname.'.fa', -format => 'fasta');
open(my $locs, '>'.shift) or die "Could not write out locs file: $!";

print "##### Making variant from $infile with ", $divRate*100, "% heterozygosity\n";

# parse no change file and store values
my @noChange = ();
while(<$noChangeFh>) {
	my $curLine = $_;
	chomp $curLine;
	next if length $curLine <= 0;
	my @tmpArr = split('-', $curLine);
	push(@noChange, [$tmpArr[0], $tmpArr[1]]);
}
if(@noChange) {
	@noChange = sort({$a->[0] <=> $b->[0]} @noChange);
}

# set up array representation of new haplotype
my $totalBases = length($in);
my @newHap = ();
for(my $i = 1; $i <= $totalBases; $i++) {
	push(@newHap, $i);
}

# undef the untouchable areas
for my $curArr (@noChange) {
	my $startNoChange = $curArr->[0];
	my $endNoChange = $curArr->[1];
	for(my $i = $startNoChange - 1; $i < $endNoChange; $i++) {
		$newHap[$i] = -1;
	}
}

# introduce variation
my $totalChanges = int($totalBases * $divRate);
my $changeCounter = $totalChanges;
while($changeCounter > 0) {
	if(rand() < $indelVsnp) {
		introduceIndel();
	} else {
		introduceSnp();
	}
}
print "Variation generation complete\n";

my $refPos = 1;
my $finalChange = 0;
my $complete = 0;
my $outHap = '';
print "Building variant\n";
my $newHapCount = scalar @newHap;
for my $curSpot (@newHap) {
	my $locsInfo = "$refPos\n";
	if($curSpot =~ /^INS/) {
		$outHap .= substr($in, $refPos - 1, 1);
		my $insLen = (split(' ', $curSpot))[1];
		for(my $i = 0; $i < $insLen; $i++) {
			$locsInfo .= "--\n";
			$outHap .= generateBase();
		}
		$finalChange = $refPos;
	} elsif($curSpot =~ /^DEL/) {
		$locsInfo = "";
		$finalChange = $refPos;
    } elsif($curSpot =~ /^SNP/) {
        $finalChange = $refPos;
        my $old = substr($in, $refPos - 1, 1);
        my $new = snpify($old);
        $outHap .= $new;
        $locsInfo = "$refPos\t$old->$new\n";
	} else {
		$outHap .=  substr($in, $refPos - 1, 1);
	}
	print $locs $locsInfo;
	$refPos++;
	my $percDone = $refPos * 100 / $newHapCount;
	if($percDone >= $complete) {
		print $complete, "% complete\n";
		$complete = int($percDone + 1);
	}
}

print "##### Finished generating variant from $infile #####\n";
$out->write_seq(Bio::Seq->new(-seq => $outHap,
			-id => $outname));
my $percChange = $totalChanges - $changeCounter;
$percChange *= 100;
$percChange /= $totalBases;
print $percChange, "% of genome changed\n";
print "$finalChange of ", $totalBases, " was the last nucleotide to be changed.\n";

close($locs);

##### SUBROUTINES #####
sub usage {
	print "USAGE:  HapMaker.pl [reference haplotype] [no-change file] [divergence rate] ",
		  "[output haplotype name] [locs filename]\n";
	exit;
}

# both introduce subroutines must introduce their respective variations to the new haplotype
sub introduceIndel {
	my $stretch = rand() < 0.291 ? 1 : int(rand($indelLength-1)) + 2;
    # Mills et al. (2006) suggest that 29.1% of indels are single base
	$changeCounter -= $stretch;
	if(rand() < 0.5) {	# insertion
		while(1) {
			my $curPos = int(rand($totalBases));
			next if($curPos + 1 != $newHap[$curPos]);
			$newHap[$curPos] = "INS $stretch";
			last;
		}
	} else {	# deletion
		my $limit = $totalBases - $stretch + 1;
		while(1) {
			my $curPos = int(rand($limit));
			next if($curPos + 1 != $newHap[$curPos]);
			for(my $i = 0; $i < $stretch; $i++) {
				$newHap[$curPos] = "DEL";
				$curPos++;
			}
			last;
		}
	}
}

sub introduceSnp {
	while(1) {
		my $curPos = int(rand($totalBases));
		next if($curPos + 1 != $newHap[$curPos]);
		$newHap[$curPos] = "SNP";
		last;
	}
	$changeCounter--;
}

sub snpify {
	my $curBase = uc(shift);
	if(rand() < $transWt) {
		return transition($curBase);
	} else {
		return transvert($curBase);
	}
}

sub transition {
	my $curBase = shift;
	if($curBase eq 'A') {
		return 'G';
	} elsif($curBase eq 'C') {
		return 'T';
	} elsif($curBase eq 'G') {
		return 'A';
	} elsif($curBase eq 'T') {
		return 'T';
	}
	return generateBase();	# N case
}

sub transvert {
	my $curBase = shift;
	if($curBase eq 'A') {
		return rand() < 0.5 ? 'T' : 'C';
	} elsif($curBase eq 'C') {
		return rand() < 0.5 ? 'A' : 'G';
	} elsif($curBase eq 'G') {
		return rand() < 0.5 ? 'T' : 'C';
	} elsif($curBase eq 'T') {
		return rand() < 0.5 ? 'A' : 'G';
	}
	return generateBase();
}

sub generateBase {
	my $curRoll = int(rand(4));
	if($curRoll == 0) {
		return 'A';
	} elsif($curRoll == 1) {
		return 'C';
	} elsif($curRoll == 2) {
		return 'G';
	}
	return 'T';
}
