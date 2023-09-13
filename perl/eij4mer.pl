#!/usr/bin/perl

#############################################################################
#
#	eij4mer.pl
#
#	count the number of teramers flanking introns
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
# Usage: eij4mer.pl *.eij
#
#############################################################################

use strict;

my %tetra;

while (<>) {
	next if (/^#/);
	my @a = split;
	my ($don, $acc) = split(/\.+/, pop(@a));
	$tetra{$don . $acc}++;
}

foreach my $q (sort keys(%tetra)) {
	printf("%s\t%d\n", $q, $tetra{$q});
}

