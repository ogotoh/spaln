#!/usr/bin/perl

#############################################################################
#
#	make_eij.pl
#
#	make *.eij by mapping nucleotie transcripts onto genome
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
#
#
#	Usage: make_eij.pl -d G_g [-Tssp other_spaln_option] C_c.cf
#
#############################################################################

use lib "$ENV{HOME}/perl";      # perl module
use SspTab;
use Util;

use strict;

sub usage {
	print STDERR "Usage: $0 -d G_g [spaln options] C_c.cf\n";
	print STDERR "\tNote: G_g must have been formatted with 'spaln -W -KD G_g.fna'\n";
	print STDERR "\t\tDon't put space between option key and arguement in spaln options\n";
	print STDERR "EXAMPLE:\n";
	print STDERR "\tOK: make_eij.pl -d gallgall_g -TTetrapod -t8 gallgall_c.cf\n";
	print STDERR "\tNG: make_eij.pl -d gallgall_g -T Tetrapod -t8 gallgall_c.cf\n";
	exit(1);
}

my $debug = 3;
my $sp_opt = "-Q7 -O12 -yX0 -LS -pq";
my $gnmd;
my $tgt;
my $alp = "AlnParam";

while ($_ = $ARGV[0], /^-/) {
	shift;
	if (/^-D/) {$debug = 2;}
	elsif (/^-h/) {&usage;}
	elsif (/^-d(\S*)/) {&Util::getoptarg(\$gnmd, $1);}
	elsif (/^-o(\S*)/) {&Util::getoptarg(\$tgt, $1);}
	elsif (/^-q/) {$debug = 1;}
	else {$sp_opt .= " $_";}
}
my $qry = join(' ', @ARGV);
&usage unless ($gnmd && $qry);

if ($sp_opt =~ /(-C\d+)/) {	# non-standard genetic code
	open(ALP, "> $alp") or die "Can't write to $alp !\n";
	print ALP $1, "\n";
	close (ALP);
}
my $cmd = "spaln -d$gnmd $sp_opt $qry";
if ($debug & 2) {print "$cmd\n";}
if ($debug & 1) {system($cmd);}
my $dot = rindex($qry, '.');
$qry = substr($qry, 0, $dot) if ($dot > 0);
my $grcd = $qry . '.grd';
$grcd .= '.gz' unless (-s $grcd);
die "no spaln output!" unless ($debug == 2 || -s $grcd);
$tgt = $qry . '.eij' unless ($tgt);
$cmd = "sortgrcd -F2 -O15 $grcd > $tgt";
if ($debug & 2) {print "$cmd\n";}
if ($debug & 1) {system($cmd);}

