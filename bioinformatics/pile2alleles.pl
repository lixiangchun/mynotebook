#!/usr/bin/perl
# Markus Rasmussen (mostly borrowed @ Seqanswers)
# samtools pileup -vcf refgenome.fasta reads.bam | perl THIS.pl > detailed (unannotated) snv data

use warnings;
use strict;
use PerlIO::gzip;
@ARGV == 1 or die "Need two input arguments.";

my %snps;
my %iupac;
$iupac{'A'}{'R'} = 'G';
$iupac{'A'}{'W'} = 'T';
$iupac{'A'}{'M'} = 'C';
$iupac{'A'}{'Y'} = 'C|T';
$iupac{'A'}{'K'} = 'G|T';
$iupac{'A'}{'S'} = 'G|C';
$iupac{'G'}{'R'} = 'A';
$iupac{'G'}{'S'} = 'C';
$iupac{'G'}{'K'} = 'T';
$iupac{'G'}{'Y'} = 'C|T';
$iupac{'G'}{'W'} = 'A|T';
$iupac{'G'}{'M'} = 'A|C';
$iupac{'C'}{'Y'} = 'T';
$iupac{'C'}{'S'} = 'G';
$iupac{'C'}{'M'} = 'A';
$iupac{'C'}{'R'} = 'A|G';
$iupac{'C'}{'W'} = 'A|T';
$iupac{'C'}{'K'} = 'G|T';
$iupac{'T'}{'Y'} = 'C';
$iupac{'T'}{'W'} = 'A';
$iupac{'T'}{'K'} = 'G';
$iupac{'T'}{'R'} = 'A|G';
$iupac{'T'}{'S'} = 'G|C';
$iupac{'T'}{'M'} = 'A|C';

#my $filename = <STDIN>;
#$filename = chomp($filename);
open IN,"<:gzip(autopop)",$ARGV[0];
while (<IN>) 
{
	chomp;
	#if ($_ == ""){next;}
	#            chr1   10884  c      G      6     44   60    3      GGG   ??>
	if ($_ =~ /(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t(\S+)\t(\S+)\t\S+/)
	{
		my ($chr, $pos, $ref, $cons, $consQual, $depth, $baseString) = ($1,$2,$3,$4,$5,$6,$7); #(split /\t/, $_)[0,1,2,3,4,7,8];
		$ref = uc $ref;
		if ($ref eq '*'){next;}
		if ($cons eq 'N'){next;}
		my $snp = ($cons =~ m/[AGCT]/) ? $cons : $iupac{$ref}{$cons};
		$snps{$chr}{$pos}{'ref'} = $ref;
		$snps{$chr}{$pos}{'snp'} = $snp;
		$snps{$chr}{$pos}{'qual'} = $consQual;
		my $snpCount;
		if ($snps{$chr}{$pos}{'snp'} eq 'A') {
			$snpCount = $baseString =~ tr/Aa/Aa/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'C') {
			$snpCount = $baseString =~ tr/Cc/Cc/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'G') {
			$snpCount = $baseString =~ tr/Gg/Gg/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'T') {
			$snpCount = $baseString =~ tr/Tt/Tt/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'A|C') {
			$snpCount = $baseString =~ tr/AaCc/AaCc/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'A|G') {
			$snpCount = $baseString =~ tr/AaGg/AaGg/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'A|T') {
			$snpCount = $baseString =~ tr/AaTt/AaTt/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'G|C') {
			$snpCount = $baseString =~ tr/GgCc/GgCc/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'G|T') {
			$snpCount = $baseString =~ tr/GgTt/GgTt/;
		} elsif ($snps{$chr}{$pos}{'snp'} eq 'C|T') {
			$snpCount = $baseString =~ tr/CcTt/CcTt/;
		} else {
			$snpCount = 0;
		}
		$snps{$chr}{$pos}{'depth'} = $depth;
		$snps{$chr}{$pos}{'freq'} = $snpCount;
		$snps{$chr}{$pos}{'pct'} = sprintf '%.2f', ($snpCount/$depth);
		}
	else
	{
		#print STDERR "In file $filename, line incompatible: $_ \n";
		print STDERR "Line $. incompatible: $_ \n";
	}
}
close IN;

foreach my $chr (sort keys %snps) {
	foreach my $pos (sort {$a <=> $b} keys %{ $snps{$chr} }) {
		print "$chr\t$pos\t${snps{$chr}{$pos}{'ref'}} > ${snps{$chr}{$pos}{'snp'}}\t${snps{$chr}{$pos}{'qual'}}\t${snps{$chr}{$pos}{'depth'}}\t${snps{$chr}{$pos}{'freq'}}\t${snps{$chr}{$pos}{'pct'}}\n";
	}
}
