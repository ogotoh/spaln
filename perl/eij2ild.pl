#!/usr/bin/perl

#############################################################################
#
#	eij2ild.pl
#
#	make ild (intron-legth-distribution) file from *.eij
#
#
#	Osamu Gotoh, ph.D.      (-2001)
#	Saitama Cancer Center Research Institute
#	818 Komuro, Ina-machi, Saitama 362-0806, Japan
#
#	Osamu Gotoh, Ph.D.      (2001-2023)
#	National Institute of Advanced Industrial Science and Technology
#	Computational Biology Research Center (CBRC)
#	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
#
#	Osamu Gotoh, Ph.D.      (2003-)
#	Department of Intelligence Science and Technology
#	Graduate School of Informatics, Kyoto University
#	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
#
#	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
#
#	Usage: eij2ild.pl xxx.eij ----> stdout
#
#############################################################################

use strict;

my %ilen;
$_ = $ARGV[0];

die "$_ not found or empty !\n" unless (-s $_);

while (<>) {
	next if (/^\#/);
	my @a = split;
	my $l = $a[-2];
	++$ilen{$l};
}

print "ilen\tfreq\n";
foreach (sort {$a <=> $b;} keys(%ilen)) {
	print "$_\t$ilen{$_}\n";
}
