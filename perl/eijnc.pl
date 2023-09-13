#!/usr/bin/perl

#############################################################################
#
#	eijnc.pl
#
#	generate gapless MSA surrounding eij(s)
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
# Usage: eijnc.pl -[3|5|j|b|i|c] [-dG_g] [-a] xxx.eij
#
#############################################################################

use strict;

my ($Don, $Acc, $Join, $Intron, $Branch, $Compos) = (1, 2, 3, 4, 5, 6);
my @code = ('', '5', '3', 'J', 'I', 'B', 'C');
my @lwing = (0, 3, 29, 30, 0, 60, 0);
my @rwing = (0, 21, 3, 30, 0, 0, 0);
my @nlen = (0, 24, 32,120, 0, 60, 0);
my $Bndry = $Don;
my $Debug;
my @data;
my $uniq = 1;
my %count;
my ($Dbs, $nlen, $wing, $stepup, $chrnbr);
my $test_seq_len;

sub usage {
	print STDERR "Usage:\teijnc.pl -[3|5|j|b|i|c] [-dG_g] [-a] xxx.eij\n";
	print STDERR "\t3: 3', 5: 5', j: 5'+3', b: bp, i: intron, c: composition\n";
	print STDERR "\ta: don't remove redundancy\n";
	exit (1);
}

while ($_ = $ARGV[0], /^-/) {
        shift;
        last if /^-$/;
	/^-a/ && ($uniq = 0);
	/^-d(\S+)/ && ($Dbs = $1);
	/^-d$/ && ($Dbs = shift);
	/^-C/ && ($stepup = 1);
	/^-n/ && ($chrnbr = 1);
	/^-3/ && ($Bndry = $Acc);
	/^-5/ && ($Bndry = $Don);
	/^-j/ && ($Bndry = $Join);
	/^-i/ && ($Bndry = $Intron);
	/^-b/ && ($Bndry = $Branch);
	/^-c/ && ($Bndry = $Compos);
	/^-h/ && (&usage());
	/^-w(\d+)/ && ($wing = $1);
	/^-w$/ && ($wing = shift);
	/^-D(\S*)/ && ($Debug = $1? $1: 1);
}

my $wms1 = $wing - 1;
my $wdth = ($Bndry == $Intron)? 60: $wing + $wing;

if ($Bndry == $Branch) {
	$wms1 = $wms1 + $wms1 + 1;
	$wing = 0;
}

my $input = $ARGV[0];
my $genspc = substr($input, 0, 8);
my $dot = rindex($input, '.');
$input = substr($input, 0, $dot) if ($dot > 0);
my @entries;
my $chr;

$Dbs = $genspc if (!$Dbs && $genspc =~ /^[a-z]/);
$Dbs .= '_g' if ($Dbs && substr($Dbs, 8, 2) ne '_g');
my $dbs = $Dbs;
my $pgs;

while (<>) {
    next if (/^#/);
    my ($junc, $rvs, $l, $r) = split;
    unless ($Dbs) {
        my $gs = lc(substr($junc, 0, 8));
	if ($pgs ne $gs) {
	    &output() if ($pgs);
	    $pgs = $gs;
	    $dbs = $gs . "_g";
	}
    }
    $rvs = $rvs eq '-';
    if ($l > $r) {
	my $t = $l;
	$l = $r;
	$r = $t;
    }
    $chr = $junc if ($junc ne $chr);

    my ($dl, $dr, $al, $ar);
    if ($rvs && $Bndry != $Intron) {
	$al = $l - $rwing[$Bndry];
	$ar = $l + $lwing[$Bndry] - 1;
	$dl = $r - $rwing[$Bndry] + 1;
	$dr = $r + $lwing[$Bndry];
        next if ($al < 0);
    } else {
	$dl = $l - $lwing[$Bndry];
	$dr = $l + $rwing[$Bndry] - 1;
	$al = $r - $lwing[$Bndry] + 1;
	$ar = $r + $rwing[$Bndry];
        next if ($dl < 0);
    }
    my $gene;
    my $sid;
    my $sgn = $rvs? '<': '';
    if ($Bndry == $Intron) {
	$gene = sprintf("\$%s %d %d %s", $chr, $dl, $ar, $sgn);
	$sid  = sprintf("%s%d.%d", $chr, $l, $r);
    } elsif ($Bndry == $Don) {
        $gene = sprintf("\$%s %d %d %s", $chr, $dl, $dr, $sgn);
	$sid  = sprintf("%s%d", $chr, $l);
    } elsif ($Bndry == $Acc || $Bndry == $Branch) {
	$gene = sprintf("\$%s %d %d %s", $chr, $al, $ar, $sgn);
	$sid  = sprintf("%s%d", $chr, $r);
    } elsif ($Bndry == $Join) {
	if ($rvs) {
	    $gene = sprintf("\$%s %d %d %d %d <", $chr, $al, $ar, $dl, $dr);
	} else {
	    $gene = sprintf("\$%s %d %d %d %d", $chr, $dl, $dr, $al, $ar);
	}
	$sid  = sprintf("%s%d.%d", $chr, $l, $r);
    } elsif ($Bndry == $Compos) {
	$gene = sprintf("\$%s %d %d %s", $chr, $l, $r, $sgn);
	$sid  = sprintf("%s%d.%d", $chr, $l, $r);
    }
    next unless (++$count{$sid} == 1);	# fist time only
    push(@entries, $gene);
}
&output();

sub output {
    my $catalog = "utncat$$";

    open(CAT, "> $catalog") or die "Can't write to $catalog !\n";
    foreach (@entries) {
	print CAT $_, "\n";
    }
    close(CAT);

    my $cmd = "utn ";
    if ($Bndry == $Compos)      {$cmd .= "-c -O0";}
    else                        {$cmd .= "-l$nlen";}
    $cmd .= " -d$dbs -i:$catalog";

    if ($Debug) {
	print "$cmd\n";
	exit (1) if ($Debug == 1);
    }

    my $nrcd = 0;
    open(UTN, "$cmd |") or die "Can't run $cmd !\n";
    while (<UTN>) {
	if ($Bndry == $Intron) {
	    if (/^(>\S+)/) {
		printf "%s.%d\n", $1, ++$nrcd;
	    } else {
		print;
	    }
	} elsif ($Bndry == $Compos) {
	    my @a = split;
	    my $sid = sprintf("%s.%d", $a[0], ++$nrcd);
	    printf "%-23s\t%7d\t%7.2f\t%7.2f\n", $sid, $a[3], $a[4], $a[5];
	} elsif (/^>(\S+)/) {
	    ++$nrcd;
	    $chr = $1;
	} elsif (/N/ || /-/ || length() <= $nlen) {	# unfinished
	    next;
	} else {
	    chomp;
	    push(@data, sprintf("%s| %s.%d", $_, $chr, $nrcd));
	}
    }
    close(UTN);
    unlink($catalog);
    @entries = ();
}

if (@data) {
    my $nseq;
    my $maxr;
    if ($uniq) {
	my @dt = ();
	my %ht = ();
	foreach (@data) {
	    ++$nseq;
	    my ($sq) = split;
	    push(@dt, $_) unless ($ht{$sq}++);
	    $maxr = $ht{$sq} if ($ht{$sq} > $maxr);
	}
	@data = @dt;
    }
    my $unseq = @data;
    printf ">%s [%d]", $input, $unseq;
    printf "\t%d %d", $nseq, $maxr if (defined($nseq));
    print "\n\n";
    foreach (@data) {
	print "       1 $_\n";
    }
}

