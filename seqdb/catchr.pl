#!/usr/bin/perl

# Catenate a single or separate chromosomal sequences in a single multi-fasta file
# Copyright(c) Osamu Gotoh <<o.gotoh@i.aist.go.jp>>
#
# Input: FASTA, Genbank/DDBJ, or EMBL format sequence files
# Output: to STDOUT in FASTA format, in which
# from the top of sequence ID to the 'Keyword' (inclusive) will be replaced by 'Prefix'
#
# Example:
# catch,pl -wAADE -hDrospseu AAED.gbk 
# will converte 
#
# Input: AAED.gbk
#--------------------------------------------------------------------------------
# LOCUS       AADE02000001           28560 bp    DNA     linear   INV 26-JUN-2018
# DEFINITION  Drosophila pseudoobscura pseudoobscura strain MV2-25 contig_1,
#            whole genome shotgun sequence.
# annotations
# ORIGIN
# sequence
# //
# LOCUS       AADE02000002           22723 bp    DNA     linear   INV 26-JUN-2018
# ...
#--------------------------------------------------------------------------------
#
# to
#
#--------------------------------------------------------------------------------
# >Drospseu02000001
# sequence
# >Drospseu02000002
# ...
#--------------------------------------------------------------------------------
#
# <2019-03-28>

use strict;

my $nbr;
my $oriid = 1;
my $addnbr;
my $idword;
my $prefix;
my $omit_old_id;
my $Un = "Un";
my $field = -1;
my $addchr;
my $gbk;
my $embl;
my $anno;

sub usage {
    print "Usage:\n";
    print "\tcatchr.pl [-wKeyword] [-hPrefix] [-n] [-o] input_files > target\n";
    exit (1);
}

while ($_ = $ARGV[0], /^-/) {
    shift;
    /^-c(\S*)/	&& ($addchr = $1);
    /^-f(\d*)/	&& ($field = $1);
    /^-h$/	&& &usage;
    /^-h(\S+)/	&& ($prefix = $1);
    /^-i/	&& ($oriid = 0);
    /^-n/	&& ($addnbr = 1);
    /^-o/	&& ($omit_old_id = 1);
    /^-u(\S+)/	&& ($Un = $1);
    /^-w(\S*)/	&& ($idword = $1);
    /^-\?/	&& &usage;
}

while (<>) {
    if (/^>/) {
	my ($gid, $rest) = split(' ', $_, 2);
	tr/|:/  /;
	if ($field >= 0) {
	    $gid =~ tr/|:/  /;
	    my @a = split(' ', $gid);
	    $gid = $field? $a[$field]: substr($a[0], 1);
	    print ">$prefix$gid $rest";
	} elsif ($gid =~ /^>(\d+)$/) {
	    print ">$prefix", $addchr, $1, " $rest";
	} elsif ($idword && (($oriid && /^>(.*)$idword\s*(\S+)(.*)/o) ||
		(!$oriid && $gid =~ /^>(.*)$idword\s*(\S+)(.*)/o))) {
	    ++$nbr;
	    my $id = $2;
	    chop($id) if (substr($id, -1, 1) eq ',');
	    print ">$prefix$id";
	    print "_$nbr" if ($addnbr);
	    my $oldid = $1;
	    if ($omit_old_id) {
		(undef, $oldid) = split(' ', $1, 2);
	    }
	    print " $oldid$3\n";
	} elsif ($addnbr) {
	    ++$nbr;
	    print ">$prefix$Un";
	    print "_$nbr ";
	    print substr($_, 1);
	} else {
	    print ">$prefix", substr($_,1);
	}
    } elsif (/^LOCUS\s+(\S+)/) {
	my $id = $1;
	if ($idword && $id =~ /^(.*)$idword\s*(\S+)(.*)/o) {
	    $id = $2;
	}
	print ">$prefix$id\n";
	$gbk = $anno = 1;
    } elsif (/^ID\s+(\S+)/) {
	my $id = $1;
	if ($idword && $id =~ /^(.*)$idword\s*(\S+)(.*)/o) {
	    $id = $2;
	}
	print ">$prefix$id\n";
	$embl = $anno = 1;
    } elsif (/^\//) {
	$anno = 1;
    } elsif ($anno) {
	$anno = 0 if (/^ORIGIN/ || /^SQ/);
    } elsif ($gbk) {
	$_ = substr($_, 10);
	tr / //d;
	print;
    } elsif ($embl) {
	chomp;
	$_ = substr($_, 5, 65);
	tr / //d;
	print $_, "\n";
    } else {
	print;
    }
    $nbr = 0 if eof;
}
