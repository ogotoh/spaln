#!/usr/bin/perl

#############################################################################
#
#	makeidx.pl
#
#	format genomic sequence or protein sequence database
#	to be used by spaln
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
#
#___________________________________________________________________________
#
#	Usage: makeidx.pl [-ianp] [-s src_dir] [-d dest_dir] [-g] fasta_file[.gz]
#	-i: Make index
#	-a: Make block info of aa sequences
#	-n: Make block info of genomic sequence for DNA queries
#	-p: Make block info of genomic sequence for aa queries
#
#	Examples:
#
#	% makeidx.pl -inp -g dictdisc_g.gf.gz
#	% makeidx.pl -ia dictdisc.faa.gz
#
###############################################################################

use strict;

my ($obj, $gz, $opt, $compress) = ('i');
my $makefile = "./Makefile";
my $src_dir = "./";

sub usage {
	print "Usage: makeidx.pl [-ianp] [-s src_dir] [-d seqdb/] [-g] fasta_file[.gz]\n";
	exit(1);
}

while ($_ = $ARGV[0], /^-/) {
	shift;
	if (/^-h/) {&usage;}
	elsif (/^-s(\S+)/) {$src_dir = $1;}	# define src_dir
	elsif (/^-s/) {$src_dir = shift;}
	elsif (/^-g/) {$opt .= '-g '; $compress = 1;}
	elsif (/^(-y\S+)/ || /^(-X\S+)/ || /^(-R\S*)/ || /^(-D\S*)/)
	    {$opt .= "$1 ";}
	elsif (/^-(\w+)/)  {$obj .= lc($1);}	# [a|n|p|np]
}

my $fsrc = shift;
unless (-s $fsrc) {
	$src_dir .= '/' if (substr($src_dir, -1, 1) ne '/');
	$fsrc = "$src_dir$fsrc";
	exit(1) unless (-s $fsrc);
}
my $sl = rindex($fsrc, '/');
my $fa = substr($fsrc, $sl + 1);

#### examine free disk space

open(DF, "df |") or die "Can't run df !\n";
while (<DF>) {
	chomp;
	my @a = split;
	my $mp = pop(@a);
	next if ($mp ne "/data");
	my $pc = pop(@a);
	chop($pc);
	die "$mp about to be full !\n" if ($pc >= 99);
}
close(DF);

#### start of main ####

print STDERR "$fa $obj\n";
my $fnm = $fa;
my @ind = ('i', 'a', 'n', 'p');
my @ext = ('idx', 'bka', 'bkn', 'bkp');

my $dot = rindex($fa, '.');
&usage unless ($dot > 0);
my $src_ext = substr($fa, $dot + 1);
my $fnm = substr($fa, 0, $dot);

if (($gz = $src_ext eq "gz")) {
	$dot = rindex($fnm, '.');
	&usage unless ($dot > 0);
	$src_ext = substr($fnm, $dot + 1);
}
my $bdy = substr($fnm, 0, $dot);

my @exp = ();
for (my $i = 0; $i < @ind; ++$i) {
	next if (index($obj, $ind[$i]) < 0);
	my $ffn = $bdy . '.' . $ext[$i];
	if ($i == 0) {
	    my $odr = $bdy . '.odr';
	    unlink($odr) if (-e $odr && -M $odr > -M $fsrc);
	}
	if (-s $ffn) {
	    if (-M $ffn > -M $fsrc) {
		unlink($ffn);
		++$exp[$i];
	    }
	} else {
	    $ffn .= '.gz';
	    if (-s $ffn) {
		if (-M $ffn > -M $fsrc) {
		    unlink($ffn);
		    ++$exp[$i];
		}
	    } else {
		++$exp[$i];
	    }
	}
}
exit(0) unless (@exp);

for (my $i = 0; $i < @ind; ++$i) {
	next unless ($exp[$i]);
	if ($i == 0)  {&makeidx($fa);}
	else {
	    my $ffn = $bdy . '.' . $ext[$i];
	    &makeblk($fa, $ffn);
	}
}

#### end of main ####

sub makeidx {
	my $fa = shift;
	my $cmd = "makdbs -K";
	$cmd .= ($src_ext eq 'fa' || $src_ext eq 'faa')? 'A': 'D';
	$cmd .= ' -g' if ($compress);
	$cmd .= " $fa";
	system($cmd);
}

sub makeblk {
	my ($fa, $blk) = @_;
	my $cmd = "makblk.pl -W$blk $opt $fa";
	print STDERR "$cmd\n";
	system($cmd);
}
