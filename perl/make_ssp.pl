#!/usr/bin/perl

#############################################################################
#
#	make_ssp.pl
#
#	make species-specific paramter files
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
# Usage: make_spp.pl -d G [-eN] [-uN] *.eij
# or	 make_spp.pl -c CDS.fna
#
# swich (level) {
#   case 1: .ild		# intron length distribution
#   case 2: .cano		# canonical sp junction
#   case 3: .SP3, .SP5		# MSA around sp junction
#   case 4: .tri3, .tri5	# position-specific 3-mer frequence
#   case 5: .pwm		# position specific weight matrix
#   case 6: 			# resevered
#   case 7: .0m3, 0m5		# 0-th Markov model around sp junction
#   case 8: .1m3, 1m5		# 1-th Markov model around sp junction
#   case 9: .2m3, 2m5		# 2-th Markov model around sp junction
#   case 10: .gci		# G+C content of intron
#   case 11: .SPb		# 5' + 3' fused MSA around sp junction
#   case 12: .iwdfq		# k-mer frequency in intron
#   case 13: .ipt		# intron potential 4-th Markov model
#   case 14: .dcf		# dicodon frequency
#   case 15: .cdp		# coding potential
#   case 16: .ildp		# fitted ild parameters
# }
#
###############################################################################

use lib "$ENV{HOME}/perl";      # perl module
use Util;

use strict;

#######	constants #######

my @lim = (16, 1024, 4096, 16384);
my @seqdb = (".", $ENV{"ALN_DBS"}, "$ENV{HOME}/seqdb");

sub usage {
	print STDERR "Usage:\n";
	print STDERR "\tmake_ssp.pl [-dG] [-eN] [-S] [-oS] [-uN] [-I] [-M] [-q] X.eij\n";
	print STDERR "or\tmake_ssp.pl -dG -S -cC\n";
	print STDERR "\nOptions:\n";
	print STDERR "\t-cC:	CDS (def =  X[0..7]_c.cf)\n";
	print STDERR "\t-dG:	genome db (def = X[0..7]_g)\n";
	print STDERR "\t-eN: make (N < 10)? [1-N]th: Nth param file\n";
	print STDERR "\t  swich (N) {\n";
	print STDERR "\t    case 1: .ild           # intron length distribution\n";
	print STDERR "\t    case 2: .cano          # canonical sp junction\n";
	print STDERR "\t    case 3: .SP3, .SP5     # MSA around sp junction\n";
	print STDERR "\t    case 4: .tri3, .tri5   # position-specific 3-mer frequence\n";
	print STDERR "\t    case 5: .pwm           # position specific weight matrix\n";
	print STDERR "\t    case 7: .0m3, 0m5      # 0-th odr Markov model around sp junction\n";
	print STDERR "\t    case 8: .1m3, 1m5      # 1-th odr Markov model around sp junction\n";
	print STDERR "\t    case 9: .2m3, 2m5      # 2-th odr Markov model around sp junction\n";
	print STDERR "\t    case 10: .gci          # G+C content of intron\n";
	print STDERR "\t    case 11: .SPb          # 5' + 3' fused MSA around sp junction\n";
	print STDERR "\t    case 12: .iwdfq        # k-mer frequency in intron\n";
	print STDERR "\t    case 13: .ipt          # intron potential 4-th odr Markov model\n";
	print STDERR "\t    case 14: .dcf          # dicodon freqency\n";
	print STDERR "\t    case 15: .cdp          # coding potential 5-th odr Markov model\n";
	print STDERR "\t    case 16: .ildp         # fitted ild parameters\n";
	print STDERR "\t  }\n";
	print STDERR "\t-oS: output label (def = G[0..7])\n";
	print STDERR "\t-q:  quiet mode\n";
	print STDERR "\t-uN: remove N % similar sequences\n";
	print STDERR "\t-I:  reserve intermediate files except MSA\n";
	print STDERR "\t-JN: minimum ORF length (90)\n";
	print STDERR "\t-M:  reserve MSA\n";
	print STDERR "\t-S:  generate spaln-compatible files\n";
	print STDERR "\nExamples:\n";
	print STDERR "\tmake_ssp.pl -d X_g -S -e9,16 X.eij\n";
	print STDERR "\tmake_ssp.pl -d X_g -e9,13,15 -J120 X.eij\n";
	print STDERR "\tmake_ssp.pl -d X_g -S -J120 -c Y_c.cf\n";
	print STDERR "\nComments:\n";
	print STDERR "\tAssuming that the genomic sequence X_g.gf(.gz) have been\n";
	print STDERR "\tformatted for DNA queries in the 'Seqdb' directory,\n";
	print STDERR "\tand Y_c.cf(.gz) exists in the current directory,\n";
	print STDERR "\tthis command makes X.ild, X.cano, etc.\n";
	print STDERR "\tFor -e16 option to work, 'fitild' must have been installed,\n";
	print STDERR "\tand 'IldModel.txt' must exist in the 'Table' directory.\n";
	print STDERR "\t-S option generates parameter files used by Spaln.\n";
	exit(1);
}

#######	global variables #######

my $debug = 3;
my $gnm;
my $cds;
my $level = 9;
my $pcntid = 0;
my $reserve_msa = 0;
my $reserve_imf = 0;
my $label;
my $min_orf = 90;
my $spaln_f = 0;
my $alp = "AlnParam";

#################################################################
#
#	main
#
#################################################################

while ($_ = $ARGV[0], /^-/) {
	shift;
	/^-h/	&& &usage();		# help
	/^-I/	&& ($reserve_imf = 1);	# retain intermediate files
	/^-M/	&& ($reserve_msa = 1);	# retain MSA
	/^-q/	&& ($debug &= 1);	# quiet
	/^-S/	&& ($spaln_f = 1);	# files used by spaln
	if (/^-D(\S*)/)	{&Util::getoptarg(\$debug, $1, 2);}	# debug
	if (/^-c(\S*)/) {&Util::getoptarg(\$cds, $1);}		# CDS
	if (/^-d(\S*)/) {&Util::getoptarg(\$gnm, $1);}		# genome
	if (/^-e(\S*)/) {&Util::getoptarg(\$level, $1);}	# job
	if (/^-J(\d*)/) {&Util::getoptarg(\$min_orf, $1);}	# minimum ORF
	if (/^-o(\S*)/) {&Util::getoptarg(\$label, $1);}	# output name
	if (/^-u(\d*)/) {&Util::getoptarg(\$pcntid, $1);}	# remove % similar
}

my $eij = shift;
&usage() unless (($eij && -s $eij) ^ ($cds && -s $cds));	# exclusive or
if ($eij) {&test_eij($eij) or die "$eij is not proper input !\n";}

my $eijgz;
$eij = $cds if ($cds);			# temporal
if (substr($eij, -3, 3) eq '.gz') {
	$eijgz = $eij;
	&System("gunzip $eij");
	$eij = substr($eijgz, 0, -3);
}

my $sl = rindex($eij, '/') + 1;
my $dot = rindex($eij, '.');
$dot = length($eij) if ($dot < $sl);
my $path = substr($eij, 0, $sl);
my $tbdy = substr($eij, $sl, $dot - $sl);
my $tid = (substr($tbdy, -2, 2) eq '_c')? substr($tbdy, 0, -2): $tbdy;
if ($cds) {
	$cds = $eij;
	$eij = undef;
}

$gnm = $tid . '_g' unless ($gnm);
my $genspc = (substr($gnm, -2, 2) eq '_g')? substr($gnm, 0, -2): $gnm;
my $cf = $cds? $cds: "$tbdy.cf";
my $cfgz = "$cf.gz";
my $outpfx = $path . ($label? $label: $genspc);

if ($cds) {			# coding potential
	&make_cdp();
	exit (0);
}

my ($num) = split(' ', `wc $eij`);
my $ndate = -M $eij;
my $idx;
my $seqdb;
foreach (@seqdb) {
	$_ .= '/' unless (substr($_, -1, 1) eq '/');
	$idx = $_ . "$gnm.idx";
	$idx .= '.gz' unless (-s $idx);
	next unless (-s $idx);
	$seqdb = $_ . $gnm;
	last;
}
die "$idx not found !\n" unless ($seqdb);
my $gdate = -M $idx;
$ndate = $gdate if ($ndate > $gdate);
my @jobs = split(/,/, $level);
my $sp5;
my $sp3;
my $tri5;
my $tri3;

foreach $level (@jobs) {
    if ($level < 10) {
	next if ($level < 1);	# eij -> ild
	&make("eij2ild.pl", $eij, "$outpfx.ild");

	next if ($level < 2);	
	die "Too few samples !\n" if ($num < $lim[0]);	# eij -> cano
	&make("canonical.pl", $eij, "$outpfx.cano");

	next if ($level < 3);	# eij -> SPx
	my $pwm = "$outpfx.pwm";
	$tri5 = "$outpfx.tri5";
	$tri3 = "$outpfx.tri3";
	my $spx = $level == 3 || $level > 6;
	if ($spx || !-s $tri5 || -M $tri5 > -M $eij) {
	    &makeSPc($outpfx, '5');
	}
	if ($spx || !-s $tri3 || -M $tri3 > -M $eij) {
	    &makeSPc($outpfx, '3');
	}

	next if ($level < 4);	# SPx -> trix
	next if ($num < $lim[1]);
	$sp5 = &gz_or_not("$outpfx.SP5");
	if (!-s $tri5 || (-s $sp5 && -M $tri5 > -M $sp5)) {
	    my $cmd = "npssm -t $tri5 $sp5";
	    &System($cmd);
	}

	$sp3 = &gz_or_not("$outpfx.SP3");
	if (!-s $tri3 || (-s $sp3 && -M $tri3 > -M $sp3)) {
	    my $cmd = "npssm -t $tri3 $sp3";
	    &System($cmd);
	}

	next if ($level < 5);	# trix -> pwm
	if (!$spaln_f && &make("npssm -m0 -l15 -u31 -f", $tri3, $pwm)) {
	    my $cmd = "npssm -m0 -l1 -u16 -f $tri5 >> $pwm";
	    &System($cmd);
	}

	next if ($level < 6);	# idx -> wdfq
	my $wdfq = &make_wdfq;

	if ($level > 8 && $num >= $lim[3]) {
	    &mksig(2, $wdfq, $tri5, $sp5, $spaln_f? "Splice5": "$outpfx.2m5");
	    &mksig(2, $wdfq, $tri3, $sp3, $spaln_f? "Splice3": "$outpfx.2m3");
	} elsif ($level > 7 && $num >= $lim[2]) {
	    &mksig(1, $wdfq, $tri5, $sp5, $spaln_f? "Splice5": "$outpfx.1m5");
	    &mksig(1, $wdfq, $tri3, $sp3, $spaln_f? "Splice3": "$outpfx.1m3");
	} elsif ($level > 6 && $num >= $lim[1]) {
	    &mksig(0, $wdfq, $tri5, $sp5, $spaln_f? "Splice5": "$outpfx.0m5");
	    &mksig(0, $wdfq, $tri3, $sp3, $spaln_f? "Splice3": "$outpfx.0m3");
	}
    } elsif ($level == 10) {
	&make("eijnc.pl -i -d$gnm", $eij, "$outpfx.gci");
    } elsif ($level == 11) {
	&makeSPc($outpfx, 'b');
    } elsif ($level == 12 || $level == 13) {
	my $iwdfq = "$outpfx.iwdfq";
	&make("kmers -d$gnm -w6 -l8 -r16 -e ", $eij, $iwdfq);
	return if ($level == 12);
	my $wdfq = &make_wdfq;	# idx -> wdfq
	my $ipt = "$outpfx.ipt";	# intron potential
	if (!-s $ipt || -M $ipt > -M $iwdfq) {
	    my $cmd = "exinpot -d$gnm -i $iwdfq -g $wdfq -e $eij -b $ipt";
	    &System($cmd);
	    if ($spaln_f && -s $ipt) {rename($ipt, "IntronPotTab.dat");} 
	}
	unlink($iwdfq) unless ($reserve_imf);
    } elsif ($level == 14 || $level == 15) {
	next unless (-s $cf || -s $cfgz);
	&System("gunzip $cfgz") if (!-s $cf);
	if ($level == 14) {
	    my $dcf = $genspc . ".dcf";
	    next if (-s $dcf && -M $dcf < -M $cf);
	    my $cmd = "exinpot -m5 -c -J$min_orf -O4 -b $dcf $cf";
	    &System($cmd);
	} else {	# $level == 15
	    make_cdp();
	}
    } elsif ($level == 16) {
	next unless (-s "$outpfx.ild");
	my $ildp = "$outpfx.ildp";
	next if (-s $ildp && -M $ildp < -M "$outpfx.ild");
	die "Install fitild !\n" unless (`which fitild`);
	my $alprm = `fitild -d IldModel.txt $outpfx.ild`;
	open(ALP, "> $ildp") or die "Can't write to $ildp !\n";
	print ALP $alprm, "\n";
	close(ALP);
	next unless ($spaln_f);
# make AlnParam
	my @a = split(' ', $alprm);
	splice(@a, 0, 3); splice(@a, 2, 1); splice(@a, -3);
	$alprm = sprintf("-yI\"%d %d %.4f ", $a[0], $a[1], $a[2]);
	if (@a == 6) {
	    $alprm .= sprintf("%.2f %.2f %.4f\"", $a[3], $a[4], $a[5]);
	} elsif (@a == 10) {
	    $alprm .= sprintf("%.4f %.2f %.2f %.4f ", $a[3], $a[4], $a[5], $a[6]);
	    $alprm .= sprintf("%.2f %.2f %.4f\"", $a[7], $a[8], $a[9]);
	} elsif (@a == 14) {
	    $alprm .= sprintf("%.4f %.2f %.2f %.4f ", $a[3], $a[4], $a[5], $a[6]);
	    $alprm .= sprintf("%.4f %.2f %.2f %.4f ", $a[7], $a[8], $a[9], $a[10]);
	    $alprm .= sprintf("%.2f %.2f %.4f\"", $a[11], $a[12], $a[13]);
	}
	if (open(ALP, "+< $alp")) {	# read-write
	    my @options;
	    while (<ALP>) {
		push(@options, $_) unless (/^-yI/);
	    }
	    seek(ALP, 0, 0);		# rewind
	    foreach (@options) {
		print ALP;
	    }
	    print ALP $alprm, "\n";
	    close(ALP);
	} elsif (open(ALP, "> $alp")) {
	    print ALP $alprm, "\n";
	    close(ALP);
	} else {
	    print STDERR "Can't write to $alp !\n";
	}
    }
} continue {
    if ($level < 10) {
	unless ($reserve_imf) {
	    unlink($tri5);
	    unlink($tri3);
	}
	&save_msa($sp5);
	&save_msa($sp3);
	if ($spaln_f) {
	    if (-s "Splice5.psm") {rename("Splice5.psm", "Splice5.dat");}
	    if (-s "Splice3.psm") {rename("Splice3.psm", "Splice3.dat");}
	}
   }
}

&System("gzip $eij") if ($eijgz);
exit (0);

#################################################################
#
#	subroutines
#
#################################################################

sub System {
	my $cmd = shift;
	if ($debug & 2) {print STDERR $cmd, "\n";}
	if ($debug & 1) {
	    print STDERR "$cmd has failed !\n" if (system($cmd));
	}
}

sub test_eij {
	my $eij = shift;
	open(EIJ, $eij) or return (0);
	my $good = 0;
	while (<EIJ>) {
	    next if (/^#/);
	    my @a = split;
	    $good = (($a[1] eq '+' || $a[1] eq '-') &&
		$a[-2] =~ /^\d+$/ && $a[-1] =~ /^\w\w\.\.\w\w/);
	    last;
	}
	close(EIJ);
	return ($good);
}

sub save_msa {
	my $msa = shift;
	return unless ($msa && -e $msa);
	if ($reserve_msa) {
	    &System("gzip $msa") unless ($msa =~ /\.gz/);
	} else {
	    unlink($msa);
	}
}

sub make {
	my ($cmd, $src, $tgt) = @_;
	if (!-s $tgt || -M $src < -M $tgt) {
	    $cmd .= " $src > $tgt";
	    &System($cmd);
	    return (1);
	}
	return (0);
}

sub mksig {
	my ($mm, $wdfq, $tri, $sp, $tgt) = @_;
	return (0) unless (!-s $tgt || -M $tgt > -M $sp);
	my $cmd = "npssm -m$mm -r $wdfq -f $tri -b $tgt $sp";
	&System($cmd);
}

sub gz_or_not {
	my $fn = shift;
	my $fngz = $fn . '.gz';
	if (-s $fn && -s $fngz) {
	    if (-M $fngz > -M $fn) {
		&Ststem("rm $fngz");
	    } else {
		&Ststem("rm $fn");
		$fn = $fngz;
	    }
	} elsif (-s $fngz) {
	    $fn = $fngz;
	}
	return ($fn);
}

sub makeSPc {
	my ($bdy, $c) = @_;
	my $tgt = &gz_or_not("$bdy.SP" . uc($c));
	return unless (!-s $tgt || (-M $tgt > $ndate));
	$tgt = substr($tgt, 0, -3) if (substr($tgt, -3, 3) eq '.gz');
	&System("eijnc.pl -$c -d$gnm $eij > $tgt");
	exit(1) unless (-s $tgt);	# empty
	return if ($pcntid == 0);

# remove closely similar members

	my $dvn = "$bdy.dvn" . $c;
	my $spc = "$bdy.sp" . uc($c);
	rename $tgt, $spc;
	&make("dvn -A1 -O0 -H$pcntid", $spc, $dvn);
	if (-s $dvn) {
	    my $rmv = "$bdy.rm$c";
	    my $cmd = "clade.pl -M $dvn > $rmv";
	    &System($cmd);
	    if (-s $rmv) {
		$cmd = "rdn -ce -l100 -f$rmv -L $spc > $tgt";
		&System($cmd);
	    }
	    &System("rm $dvn");
	    &System("rm $spc");
	    &System("rm $rmv");
	} else {
	    rename $spc, $tgt;
	    &System("rm $dvn") if (-e $dvn);
	}
}

sub make_wdfq {
	my $wdfq = "$outpfx.wdfq";
	my $cmd = "kmers -KD -w6 -d $gnm > $wdfq";
	&System($cmd) if (!-s $wdfq || -M $wdfq > $gdate);
	return ($wdfq);
}

sub make_cdp {
	my $cdp = $genspc . ".cdp";
	next if (-s $cdp && -M $cdp < -M $cf);
	my $wdfq = &make_wdfq;
	my $cmd = "exinpot -m5 -c -g $wdfq -J$min_orf -O4 -b $cdp $cf";
	&System($cmd);
	if ($spaln_f && -s $cdp) {rename($cdp, "CodePotTab.dat");}
}
