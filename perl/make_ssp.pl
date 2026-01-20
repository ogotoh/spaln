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
# Usage: make_spp.pl -d G [-eN] [-uN] *.eij[.gz]
# or	 make_spp.pl -c CDS_fna[.gz]
#
# swich (level) {
#   case 1: .ild		# intron length distribution
#   case 2: .cano		# canonical sp junction
#   case 3: .SP3, .SP5		# MSA around sp junction
#   case 4: .tri3, .tri5	# position-specific 3-mer frequence
#   case 5: .pwm.dat		# position specific weight matrix
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
my @gnm_ext = (".ent", ".idx", ".odr", ".seq");
my $ext_gz = ".gz";
my $ext_dgz = ".dgz";
my $avail_mem = 6;	# GB
my $uncompress_genome = 0;

sub usage {
	print STDERR "Usage:\n";
	print STDERR "\tmake_ssp.pl [-dG] [-eN] [-S] [-oS] [-uN] [-I] [-M] [-q] X.eij\n";
	print STDERR "or\tmake_ssp.pl -dG -S -cC\n";
	print STDERR "\nOptions:\n";
	print STDERR "\t-cC:	CDS (def =  X_m.cf[.gz])\n";
	print STDERR "\t-dG:	genome db (def = X_g)\n";
	print STDERR "\t-eN: make (N < 10)? [1-N]th: Nth param file\n";
	print STDERR "\t  swich (N) {\n";
	print STDERR "\t    case 1: .ild           # intron length distribution\n";
	print STDERR "\t    case 2: .cano          # canonical sp junction\n";
	print STDERR "\t    case 3: .SP3, .SP5     # MSA around sp junction\n";
	print STDERR "\t    case 4: .tri3, .tri5   # position-specific 3-mer frequence\n";
	print STDERR "\t    case 5: .pwm	   # position specific weight matrix\n";
	print STDERR "\t    case 7: .0m3, 0m5      # 0-th odr Markov model around sp junction\n";
	print STDERR "\t    case 8: .1m3, 1m5      # 1-th odr Markov model around sp junction\n";
	print STDERR "\t    case 9: .2m3, 2m5      # 2-th odr Markov model around sp junction\n";
	print STDERR "\t    case 10: .gci          # G+C content of intron\n";
	print STDERR "\t    case 11: .SPb          # 5' + 3' fused MSA around sp junction\n";
	print STDERR "\t    case 12: .iwdfq        # k-mer frequency in intron\n";
	print STDERR "\t    case 13: .ipt          # intron potential 4-th odr Markov model\n";
	print STDERR "\t    case 14: .dcf          # dicodon freqency\n";
	print STDERR "\t    case 15: .cdp,.cus     # coding potential 5-th odr Markov model\n";
	print STDERR "\t    case 16: .ildp         # fitted ild parameters\n";
	print STDERR "\t  }\n";
	print STDERR "\t-g:  gzipped output\n";
	print STDERR "\t-oS: output label (def = G[0..7])\n";
	print STDERR "\t-q:  quiet mode\n";
	print STDERR "\t-uN: remove N % similar sequences\n";
	print STDERR "\t-z:  memory in GB to allocate gzipped whole genomic seq (6)\n";
	print STDERR "\t-I:  reserve intermediate files except MSA\n";
	print STDERR "\t-JN: minimum ORF length (90)\n";
	print STDERR "\t-M:  reserve MSA\n";
	print STDERR "\t-S:  generate spaln-compatible files\n";
	print STDERR "\t-SC: generate spaln-compatible files removing intermediate files\n";
	print STDERR "\nExamples:\n" ;
	print STDERR "\tmake_ssp.pl -d X_g -S -g X.eij\n";
	print STDERR "\tmake_ssp.pl -d X_g -e9,13,15 -J120 -c Y_c.cf X.eij\n";
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
my $gzipped = 0;

#################################################################
#
#	main
#
#################################################################

while ($_ = $ARGV[0], /^-/) {
	shift;
	/^-h/	&& &usage();		# help
	/^-g/	&& ($gzipped = 1);	# gzipped
	/^-I/	&& ($reserve_imf = 1);	# retain intermediate files
	/^-q/	&& ($debug &= 1);	# quiet
	/^-M/	&& ($reserve_msa = 1);	# retain MSA
	/^-S/	&& ($spaln_f = 1);	# used by spaln/aln
	if (/^-D(\S*)/)	{&Util::getoptarg(\$debug, $1, 2);}	# debug
	if (/^-c(\S*)/) {&Util::getoptarg(\$cds, $1);}		# CDS
	if (/^-d(\S*)/) {&Util::getoptarg(\$gnm, $1);}		# genome
	if (/^-e(\S*)/) {&Util::getoptarg(\$level, $1);}	# job
	if (/^-J(\d*)/) {&Util::getoptarg(\$min_orf, $1);}	# minimum ORF
	if (/^-o(\S*)/) {&Util::getoptarg(\$label, $1);}	# output name
	if (/^-u(\d*)/) {&Util::getoptarg(\$pcntid, $1);}	# remove % similar
	if (/^-z(\d*)/) {&Util::getoptarg(\$avail_mem, $1);}	# available memory 
}

my $eij = shift;
my $eijgz;
my $eijdate;
my $cdsgz;
my $cdsdate;

&usage unless ($gnm);
my $sl = rindex($gnm, '/') + 1;
my $wkdir = ($sl >= 0)? substr($gnm, 0, $sl): "";
my $genspc = substr($gnm, $sl);
if (substr($genspc, -2, 2) eq '_g') {$genspc = substr($genspc, 0, -2);}
else	{$gnm .= '_g';}

if ($eij) {
	if ($wkdir) {
	    $sl = rindex($eij, '/') + 1;
	    my $path = substr($eij, 0, $sl);
	    $eij = $wkdir . $eij if ($wkdir ne $path);
	}
	my $isgz = substr($eij, -3, 3) eq $ext_gz;
	if (-s $eij) {
	    $eijdate = -M $eij;
	    if ($isgz) {
		$eijgz = $eij;
		$eij = substr($eijgz, 0, -3);
	    } elsif (!&test_eij()) {
		$eij = undef;
	    }
	} elsif ($isgz) {
	    $eij = substr($eij, 0, -3);
	    $eijdate = -M $eij if (-s $eij);
	} else {
	    $eijgz = $eij . $ext_gz;
	    $eijdate = -M $eijgz if (-s $eijgz);
	}
}

if ($cds) {
	if ($wkdir) {
	    $sl = rindex($cds, '/') + 1;
	    my $path = substr($cds, 0, $sl);
	    $cds = $wkdir . $cds if ($wkdir ne $path);
	}
	my $isgz = substr($cds, -3, 3) eq $ext_gz;
	if (-s $cds) {
	    $cdsdate = -M $cds;
	    if ($isgz) {
		$cdsgz = $cds;
		$cds = substr($cdsgz, 0, -3);
	    }
	} elsif ($isgz) {
	    $cds = substr($cds, 0, -3);
	    $cdsdate = -M $cds if (-s $cds);
	} else {
	    $cdsgz = $cds . $ext_gz;
	    $cdsdate = -M $cdsgz if (-s $cdsgz);
	}
}

&usage unless  ($eijdate || $cdsdate);

my $outpfx = $wkdir . ($label? $label: $genspc);
my $alp = $wkdir . "AlnParam";
my $spsig5 = $wkdir . "Splice5";
my $spsig3 = $wkdir . "Splice3";
my $ipt = $wkdir . "IntronPotTab";
my $cdp = $wkdir . "CodePotTab";
my $sps5 = $outpfx . ".SP5";
my $sps3 = $outpfx . ".SP3";
my $tri5;
my $tri3;
unless ($spaln_f) {
	$tri5 = &gz_or_not($outpfx . ".tri5", $ext_dgz);
	$tri3 = &gz_or_not($outpfx . ".tri3", $ext_dgz);
	$ipt = $outpfx . ".ipt";
	$cdp = $outpfx . ".cdp";
}
my $iptgz = $ipt . $ext_gz;
my $iptdgz = $ipt . $ext_dgz;
my $cdpgz = $cdp . $ext_gz;
my $cdpdgz = $cdp . $ext_dgz;

unless ($eijdate) {
	&make_cdp();
	&compress($eij);
	exit (0);
}

my $ndate = $eijdate;
my $idx;
my $seqdb;
foreach (@seqdb) {
	$_ .= '/' unless (substr($_, -1, 1) eq '/');
	$idx = $_ . $gnm . ".idx";
	$idx .= $ext_gz unless (-s $idx);
	next unless (-s $idx);
	$seqdb = $_ . $gnm;
	last;
}
die "$idx not found !\n" unless ($seqdb);
my $grp = $seqdb . ".grp";
if (open(GRP, $grp)) {
    while (<GRP>) {
	my @grp = split;
	if ($grp[$#grp] eq "E_O_F") {
	    $uncompress_genome = $grp[0] > $avail_mem * 1.e9;
	    last;
	}
    }
    close(GRP);
}
my $gdate = -M $idx;
$ndate = $gdate if ($ndate > $gdate);
my $wdate = $ndate;
my @jobs = $spaln_f? (9, 13, 16): split(/,/, $level);
push(@jobs, 15) if ($cds);
my %jobs = map { $_ => 1 } (@jobs);

foreach $level (@jobs) {
    if ($level < 10) {
	next if ($level < 1);
# eij -> ild
	my $ild = $outpfx . ".ild";
	my $ilddat = $ild . ".dat";
	my $ilddgz = $ild . $ext_dgz;
	my $cmd = "eij2ild ";
	$cmd .= "-g -b $ild ";
	if (-s $ild && -M $ild < $ndate) {
	    $cmd .= $ild;
	    &System($cmd);
	    unlink($ild);
	} elsif (-s $ilddat && -M $ilddat < $ndate) {
	    $cmd .= $ilddat;
	    &System($cmd);
	    unlink($ilddat);
	} else {
	    unlink($ild) if (-e $ild);
	    unlink($ilddat) if (-e $ilddat);
	    unlink($ilddgz) if (-e $ilddgz);
	    &expandeij() if (-s $eijgz);
	    next unless (-s $eij);
	    $cmd .= $eij;
	    &System($cmd);
	}
	next if ($level < 2);

# eij -> cano
	my $cano = $outpfx . ".cano";
	my $num = &renew_cano($cano) if (-s $cano && -M $cano < $eijdate);
	unless ($num) {
	    &expandeij() if (-s $eijgz);
	    next unless ($eij && -s $eij);
	    $num = &cannonical($eij, $cano);
	    next unless ($num);
	}
	if ($num < $lim[0]) {
	    print STDERR "$eij: Too few samples !\n";
	    last;
	}
	next if ($level < 3);

	unless ($spaln_f) {
# SPx -> trix
	    my $pwm = $outpfx . ".pwm";	# text
	    my $pwmdt = $pwm  . ".dat";	# binary
	    unless (-s $tri5 && -M $tri5 < $ndate) {
		unlink($tri5) if (-e $tri5);
		my $sp5 = &makeSPc($sps5);
		my $cmd = "npssm -t $tri5 $sp5";
		$tri5 .= $ext_dgz if (substr($tri5, -4, 4) ne $ext_dgz);
		&System($cmd);
	    }

	    unless (-s $tri3 && -M $tri3 < $ndate) {
		unlink($tri3) if (-e $tri3);
		my $sp3 = &makeSPc($sps3);
		my $cmd = "npssm -t $tri3 $sp3";
		$tri3 .= $ext_dgz if (substr($tri3, -4, 4) ne $ext_dgz);
		&System($cmd);
	    }

	    next if ($level < 5);
# trix -> pwmdt
	    if (-e $pwm) {
		if (-s $pwm && -M $pwm < $ndate) {
		    my $cmd = "npssm -n$num -p $pwm -b $pwmdt";
		    &System($cmd);
		}
		unlink($pwm);
	    }
	    unless (-s $pwmdt && -M $pwmdt < -M $tri3) {
		my $tmp3 = "tmp3$$.dat";
		my $cmd = "npssm -m0 -p -l15 -u31 -i $tri3 -b $tmp3";
		&System($cmd);
		my $tmp5 = "tmp5$$.dat";
		$cmd = "npssm -m0 -p -l1 -u16 -i $tri5 -b $tmp5";
		&System($cmd);
		$cmd = "npssm -p '$tmp3<$tmp5' -b $pwmdt";
		&System($cmd);
		unlink($tmp3); unlink($tmp5); 
	    }
	}	# spaln_f
	next if ($level < 6);

# idx -> wdfq
	my $wdfq = &make_wdfq(".wdfq");

	if ($level > 8 && $num >= $lim[3]) {
	    &mksig(2, $wdfq, $sps5, $spaln_f? $spsig5: $outpfx . ".2m5", $tri5);
	    &mksig(2, $wdfq, $sps3, $spaln_f? $spsig3: $outpfx . ".2m3", $tri3);
	} elsif ($level > 7 && $num >= $lim[2]) {
	    &mksig(1, $wdfq, $sps5, $spaln_f? $spsig5: $outpfx . ".1m5", $tri5);
	    &mksig(1, $wdfq, $sps3, $spaln_f? $spsig3: $outpfx . ".1m3", $tri3);
	} elsif ($level > 6 && $num >= $lim[1]) {
	    &mksig(0, $wdfq, $sps5, $spaln_f? $spsig5: $outpfx . ".0m5", $tri5);
	    &mksig(0, $wdfq, $sps3, $spaln_f? $spsig3: $outpfx . ".0m3", $tri3);
	}
    } elsif ($level == 10) {
	expandeij() if (-s $eijgz);
	next unless ($eij);
	&uncompress_gnm(".seq");
	my $outf = $outpfx .".igc.dat";
	&make("eijnc.pl -b -c -d$gnm ", $eij, $outf); # binary output
    } elsif ($level == 11) {
	&makeSPc($outpfx . "SPb");
    } elsif (($level == 12 || $level == 13)) {
	expandeij() if (-s $eijgz);
	next unless (-s $eij);
	my $iwdfq = &make_wdfq(".iwdfq", "-l8 -r16 -e $eij");
	return if ($level == 12);
	my $wdfq = &make_wdfq(".wdfq");	# idx -> wdfq
	my $iwdate = -M $iwdfq;
	$iwdate = $wdate if ($wdate < $iwdate);
	next if ((-s $ipt && -M $ipt < $iwdate) ||
		 (-s $iptgz && -M $iptgz < $iwdate) ||
		 (-s $iptdgz && -M $iptdgz < $iwdate));
	&uncompress_gnm(".seq");
	my $cmd = "exinpot -d$gnm -i $iwdfq -r $wdfq -e $eij -b $ipt";
	$cmd .= " -g" if ($gzipped);
	$cmd .= " -O5 >> iptscore.txt";
	&System($cmd);
	&rm_ext($ipt, ".ipt") if ($spaln_f);
    } elsif (($level == 14 || $level == 15)) {
	if ($level == 14) {
	    my $dcf = $outpfx .".dcf";
	    next if (-s $dcf && -M $dcf < $cdsdate);
	    expandcds() if (-s $cdsgz);
	    next unless (-s $cds);
	    my $cmd = "exinpot -m5 -c -J$min_orf -O5 -b $dcf";
	    $cmd .= " -g" if ($gzipped);
	    $cmd .= " $cds";
	    &System($cmd);
	    &rm_ext($dcf, ".dcf") if ($spaln_f);
	} else {	# $level == 15
	    &make_cdp();
	}
    } elsif ($level == 16) {
	my $ild = $outpfx . ".ild";
	unless (-s $ild) {
	    my $ildext = $ild . ".dat" ;
	    $ildext .= $ext_gz unless (-s $ildext);
	    $ildext = $ild . $ext_dgz unless (-s $ildext);
	    next unless (-s $ildext);
	    $ild = $ildext;
	}
	my $ildp = $outpfx . ".ildp";
	my $alprm;
	if (-s $ildp && -M $ildp < -M $ild) {
	    $alprm = `cat $ildp`;
	} else {
	    die "Install fitild !\n" unless (`which fitild`);
	    $alprm = `fitild -d IldModel.txt $ild`;
	    open(ALP, "> $ildp") or die "Can't write to $ildp !\n";
	    print ALP $alprm;
	    close(ALP);
	}
	next unless ($spaln_f);
# make AlnParam
	my @a = split(' ', $alprm);
	splice(@a, 0, 3); splice(@a, 1, 1); splice(@a, -3);
	$alprm = sprintf("-yI\"%d %d %.4f ", $a[0], $a[1], $a[2]);
	if (@a == 7) {
	    $alprm .= sprintf("%.2f %.2f %.4f\"", $a[4], $a[5], $a[6]);
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
	&save_msa($sps5);
	&save_msa($sps3);
    }
}

if ($cds) {
    if (length($cds) > 5 && substr($cds, -5, 2) eq "_m") {
	my $ccf = $cds;
	substr($ccf, -4, 1) = "c";
	&compress($ccf) if (-s $ccf);
    }
    &compress($cds);
    &compress($cdp) unless ($spaln_f);
}
&compress($eij);
&compress_gnm();

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
	print "$eij is empty !\n" unless ($good);
	return ($good);
}

sub expandeij()
{
	&System("gunzip $eijgz");
	$eij = undef unless (test_eij());
	$eijgz = undef;
}

sub expandcds()
{
	&System("gunzip $cdsgz");
	$cdsgz = undef;
} 

sub cannonical {
	my($eij, $cano) = @_;
	my %tetra;

	open(EIJ, $eij) or return (0);
	my $num = 0;
	while (<EIJ>) {
	    next if (/^#/);
	    my @a = split;
	    my ($don, $acc) = split(/\.+/, pop(@a));
	    ++$tetra{$don . $acc};
	    ++$num;
	}
	close(EIJ);
	return (0) unless ($num);
	my $sl = rindex($outpfx, '/');
	my $id = substr($outpfx, $sl + 1);
	my $pct = $num / 100;
	my $cano4mer = $tetra{"GTAG"} + $tetra{"GCAG"} + $tetra{"ATAC"};
	open(CANO, "> $cano") or die "Can't write to $cano !\n";
	printf(CANO "%s\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n", 
	    $id, $num, $tetra{"GTAG"}/$pct, $tetra{"GCAG"}/$pct, 
	    $tetra{"ATAC"}/$pct, ($num - $cano4mer)/$pct);
	close(CANO);
	return ($num);
}

sub renew_cano()
{
	my $cano = shift;
	open(CANO, $cano) or return (0);
	my $old_style;
	my $id;
	my $num = 0;
	my $sls;
	my @a;
	while (<CANO>) {
	    if (/^#genespec/) {
		$old_style = 1;
		next;
	    } elsif (!$old_style) {		# new stype
		($id, @a) = split;
		$sls = rindex($id, '/');
		my $ub = substr($id, -2);
		$ub = $ub eq '_c' || $ub eq '_m';
		my $sum = shift(@a);
		return ($sum) if ($sls <= 0 && !$ub);	# as is
		$id = substr($id, $sls + 1);	# new but ealier style
		$id = substr($id, 0, -2) if ($ub);
		unlink($cano);
		open(CANO, "> $cano") or die "Can't write to $cano !\n";
		printf(CANO "%s\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n", 
		    $id, $sum, $a[0], $a[1], $a[2], $a[3]);
		close(CANO);
		return($sum);
	    } else {
		($id, @a) = split;
		last;
	    }
	}
	close(CANO);
	my $sls = rindex($id, '/') + 1;
	$id = substr($id, $sls) if ($sls > 0);
	$id = substr($id, 0, -2) if (substr($id, -2, 1) eq '_');
	my $sum = $a[0] + $a[1] + $a[2] + $a[3];
	return (0) if ($sum == 0);
	my $pct = $sum / 100;
	unlink($cano);
	open(CANO, "> $cano") or die "Can't write to $cano !\n";
	printf(CANO "%s\t%d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n", 
	    $id, $sum, $a[2]/$pct, $a[1]/$pct, $a[0]/$pct,$a[3]/$pct);
	close(CANO);
	return ($sum);
}

sub save_msa {
	my $msa = shift;
	return unless ($msa && -e $msa);
	if ($reserve_msa == 0) {
	    unlink($msa);
	} else {
	    &System("gzip $msa") unless ($msa =~ /\.gz/);
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
	my ($mm, $wdfq, $sps, $tgt, $trix) = @_;
	my $sig = $tgt;
	my $sig2 = $tgt . ".psm";		# old style
	if (-s $sig2 && -M $sig2 < $wdate) {
	    my $cmd = "npssm -f $sig2 -b $sig2";
	    $cmd .= " -g" if ($gzipped);
	    &System($cmd);
	    unlink($sig2);
	    return;
	}
	unlink($sig2) if (-e $sig2);
	$sig .= $ext_dgz;
	return if (-s $sig && -M $sig < $wdate);
	$sps = &makeSPc($sps);
	my $cmd = "npssm -m$mm -r $wdfq -b $tgt";
	$cmd .= " -i $trix" if ($trix && -s $trix);
	$cmd .= " -g" if ($gzipped);
	$cmd .= " $sps";
	&System($cmd);
}

sub gz_or_not {
	my ($fn, $ext) = @_;
	my $fngz = $fn . $ext;
	if (-s $fn && -s $fngz) {
	    if (-M $fngz > -M $fn) {
		unlink($fngz);
	    } else {
		unlink($fn);
	    }
	} elsif (-s $fngz) {
	    $fn = $fngz;
	}
	return ($fn);
}

sub makeSPc {
	my $msa = shift;
	my $c = substr($msa, -1, 1);
	my $msagz = $msa . $ext_gz;
	return ($msagz) if (-s $msagz && (-M $msagz < $ndate));
	return ($msa) if (-s $msa && (-M $msa < $ndate));

	&expandeij() if (-s $eijgz);
	die "$eij not defined" unless ($eij && -s $eij);
	&uncompress_gnm(".seq");
	&System("eijnc.pl -$c -d$gnm $eij > $msa");
	die "$msa is empty !\n" unless (-s $msa);
	return ($msa) if ($pcntid == 0);

# remove closely similar members

	my $bdy = substr($msa, 0, -1);
	my $dvn = "$bdy.dvn" . $c;
	&make("dvn -A1 -O0 -H$pcntid", $msa, $dvn);
	if (-s $dvn) {
	    my $spc = "$bdy.sp" . uc($c);
	    my $rmv = "$bdy.rm$c";
	    my $cmd = "clade.pl -M $dvn > $rmv";
	    &System($cmd);
	    if (-s $rmv) {
		$cmd = "rdn -ce -l100 -f$rmv -L $msa > $spc";
		&System($cmd);
		if (-s $spc) {
		    unlink($msa);
		    rename($spc, $msa);
		} elsif (-e $spc) {
		    unlink($spc);
		}
	    }
	    unlink($dvn);
	    unlink($rmv);
	} elsif (-e $dvn) {
	    unlink($dvn);
	}
	return ($msa);
}

sub make_wdfq {
	my ($ext, $opt) = @_;
	my $wfq = $outpfx . $ext;
	my $wdfq = $wfq . $ext_dgz;
	my $wfqdz = $wfq . ".dat.gz";
	if (-s $wfqdz && -M $wfqdz < $gdate) {
	    rename($wfqdz, $wdfq);
	} elsif (-s $wfq && -M $wfq < $gdate) {
	    my $cmd = "kmers -g -b $wfq -c $wfq";
	    &System($cmd);
	    unlink($wfq);
	} elsif (!(-s $wdfq && -M $wdfq < $gdate)) {
	    my $cmd = "kmers -KD -w6 -d $gnm $opt -g -b $wfq";
	    &System($cmd);
	}
	$wdate = -M $wdfq if ($ext eq ".wdfq");
	return ($wdfq);
}

sub make_cdp {
	my $wdfq = &make_wdfq(".wdfq");
	return if ($cdsdate < $wdate && (
		(-s $cdp && -M $cdp < $cdsdate) ||
		(-s $cdpgz && -M $cdpgz < $cdsdate) ||
		(-s $cdpdgz && -M $cdpdgz < $cdsdate)));
	&expandcds() if (-s $cdsgz);
	my $cmd = "exinpot -m5 -O5 -c -r $wdfq -J$min_orf";
	$cmd .= " -g" if ($gzipped);
	$cmd .= " -b $cdp $cds >> cdpscore.txt";
	&System($cmd);
	&rm_ext($cdp, ".cdp") if ($spaln_f);
}

sub compress {
	my $fn = shift;
	my $fngz = $fn . $ext_gz;
	if (-s $fn) {
	    if (-s $fngz) {
		if (-M $fn < -M $fngz) {
		    unlink($fngz);
		    &System("gzip $fn");
		} else {
		    unlink($fn);
		}
	    } else {
		&System("gzip $fn");
	    }
	}
}

sub compress_gnm {
	foreach my $ext (@gnm_ext) {
	    my $file = "$gnm$ext";
	    next if (-e "$file$ext_gz");
	    &System("gzip $file") if (-s $file);
	}
}

sub uncompress_gnm {
	return unless ($uncompress_genome);
	my $ext = shift;
	if ($ext) {
	    my $file = "$gnm$ext$ext_gz";
	    &System("gunzip $file") if (-s $file);
	    return;
	}
	foreach my $ext (@gnm_ext) {
	    my $file = "$gnm$ext$ext_gz";
	    &System("gunzip $file") if (-s $file);
	}
}

sub rm_ext {
	my ($bdy, $pext) = @_;
	return if (-s $bdy);
	my	$fn = $bdy . $pext;
	return if (-s $fn);
	my $dext = $gzipped? ".dgz": ".dat";
	my $sspf = $bdy . $dext;
	return if (-s $sspf);
	$fn .= $dext;
	rename($fn, $sspf) if (-s $fn);
}

