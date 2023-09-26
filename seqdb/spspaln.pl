#!/usr/bin/perl

##############################################################################
#
#	spspaln.pl: Run spaln with the speied-specific table inferred
#			from each chromosome/contig id (cid)
#
#
#	Osamu Gotoh, ph.D.      (-2001)
#	Saitama Cancer Center Research Institute
#	818 Komuro, Ina-machi, Saitama 362-0806, Japan
#
#	Osamu Gotoh, Ph.D.      (2001-)
#	National Institute of Advanced Industrial Science and Technology
#	Computational Biology Research Center (CBRC)
#	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
#
#	Osamu Gotoh, Ph.D.      (2003-2012)
#	Department of Intelligence Science and Technology
#	Graduate School of Informatics, Kyoto University
#	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
#
#	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
#____________________________________________________________________________
#
#	Usage: spspaln.pl -dgenome_gf [other spaln options] queries
#
#	The genomic fasta file must follow the following rules
#
#	The first eight charactors of cid specifies the species
#	The first letter immediately followin '>' must be in upper case
#	Other 7 letters must be lower case
#	The letters after the 8 charactors must be uniqe to each other
#
#	Example of caeneleg_g.gf:
#
#	>CaenelegMtDNA dna chromosome chromosome WBcel235 MtDNA 1 13794 1
#	CAGTAAATAGTTTAATAAAAATATAGCATTTGGGTTGCTAAGATATTATTACTGATAGAA
#	....
#	>CaenelegIII dna chromosome chromosome WBcel235 III 1 13783801 1
#	CCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAG
#	....
#
#############################################################################

use lib "$ENV{HOME}/perl";      # perl module
use strict;
use SspTab;

sub usage {
	print STDERR "spspaln.pl: species-sorted spaln\n";
	print STDERR "Usage: spspaln.pl [spaln-options] queries.fasta\n";
	exit (1);
}

my $spidn = 8;
my $dbs;
my $opt;
my $vbs = 2;

while ($_ = $ARGV[0], /^-/) {
	shift;
	if (/^-V(\d*)/) {$vbs = $1? $1: 3;}
	elsif (/^-d(\S+)/) {$dbs = $1; $opt .= "$_ ";}
	elsif (/^-D(\S+)/) {$dbs = $1; $opt .= "$_ ";}
	else {$opt .= "$_ ";}
}

my $spc = $dbs? substr($dbs, 0, $spidn): undef;

&opengnmdb(".");  # initialize gnm2tab

##### given genome db #####

if ($spc) {
	my $ssp = &sspopt($dbs);
	my $args = join(' ', @ARGV);
	my $cmd = "spaln -Q7 $ssp$opt $args";
	print $cmd, "\n" if ($vbs & 1);
	system($cmd) if ($vbs & 2);
	exit (0);
}

##### genome db specified by query ####

my $tmpf = "spspaln$$";

while (<>) {
	if (/^>/) {
	    my $sp = substr($_, 1, $spidn);
	    if ($spc ne $sp) {
		&exec_spaln() if ($spc);
		$spc = $sp;
		open(FD, "> $tmpf") or die "Can't write to $tmpf \n";
	    }
	}
	print FD;
}
&exec_spaln();

sub exec_spaln {
	close(FD);
	my $sp = lc($spc);
	my $gnm = $sp . "_g";
	my $ssp = &sspopt($gnm);
	return unless ($ssp);
	return unless (&fseqdb("$gnm.idx"));
	my $cmd = "spaln -Q7 -d$gnm $ssp$opt $tmpf";
	print $cmd, "\n" if ($vbs & 1);
	system($cmd) if ($vbs & 2);
	unlink($tmpf);
}
