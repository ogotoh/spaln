#!/usr/bin/perl

# Catenate separate chromosomal sequences in a single multi-fasta file
# Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
# <2007-07-23>

my $nbr;
my $oriid = 1;
my $addnbr;
my $idword;
my $prefix;
my $omit_old_id;
my $Un = "Un";
my $field = -1;
my $addchr;

sub usage {
    print STDERR "catchr -hPrefix [-n] [-o][-wWord] files > target\n";
    exit (1);
}

while ($_ = $ARGV[0], /^-/) {
    shift;
    /^-c(\S*)/	&& ($addchr = $1);
    /^-f(\d*)/	&& ($field = $1);
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
		(my $dummy, $oldid) = split(' ', $1, 2);
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
    } else {
	print;
    }
    $nbr = 0 if eof;
}
