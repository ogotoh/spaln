#!/usr/bin/perl

#############################################################################
#
#	canonical.pl
#
#	count the number of canonical intron boundaries
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
# Usage: canonical.pl *.eij
#
#############################################################################

use strict;

my %atac;
my %gcag;
my %gtag;
my %other;

foreach my $fn (@ARGV) {
	my $dot = rindex($fn, '.');
	my $eij = substr($fn, 0, $dot);
	open(CON, "eij4mer.pl $fn |") or die "Can't run eij4mer !\n";
	while (<CON>) {
	    my ($didi, $frq) = split;
	    my $acgt = ($didi =~ tr/ACGT//);
	    next unless ($acgt == 4);
	    if    ($didi eq "ATAC") {$atac{$eij} = $frq;}
	    elsif ($didi eq "GCAG") {$gcag{$eij} = $frq;}
	    elsif ($didi eq "GTAG") {$gtag{$eij} = $frq;}
	    else	{$other{$eij} += $frq;}
	}
	close(CON);
}

foreach my $eij (sort keys(%gtag)) {
	my $sum = ($atac{$eij} + $gcag{$eij} + $gtag{$eij} + $other{$eij}) / 100;
	print "#genespec\t   ATAC\t   GCAG\t   GTAG\t Others\t";
	print "   atac\t   gcag\t   gtag\t others\n";
	printf("%s\t%7d\t%7d\t%7d\t%7d\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n", $eij,
	    $atac{$eij}, $gcag{$eij}, $gtag{$eij}, $other{$eij},
	    $atac{$eij}/$sum, $gcag{$eij}/$sum,  $gtag{$eij}/$sum, $other{$eij}/$sum);
}

