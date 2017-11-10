#!/usr/bin/perl

use strict;

sub usage {
	print "Usage: makblk.pl [-Wx.[bkn|bkp|bka] -K[D|P|S]] [options] seq.1 seq.2 ...\n";
	print "\t or  makblk.pl [-Wx.[bkn|bkp|bka] -K[D|P|S]] [options] x.grp\n";
	print "\nOptions:\n";
	print "\t-r:\tuse reduced alphabet\n";
	print "\t-X...\ttransferred to spaln\n";
	print "-XCc:\tc=0: conti seed; c=1: spaced seed; c>1: conti + (c-1) spaced seeds\n";
	exit(0);
}

my $Debug;
my $opt;
my $MinMaxGene = 16384;
my $genlencoeff = 36;	# assume gene_size o< sqrt(genome_size)
my ($mol, $out, $ver, $grp, $ext);
my ($ktuple, $blksz, $maxgene, $xa, $xb, $xc, $gkmg, $bkmg);
my @Bitpat = ("1011,10011", "101011,1000111", "10100111,100101101", "1010011011,1010100111",
	"1001110111,100011011011", "100110110111,1010010111011",
	"1001110110111,10100101011111", "100111001101111,1010011010101111",
	"1000111101111011,1001110101111011", "101001101011111011,100011100011111111",
	"10101001110100111111,101100011011010011111", "100011110001111110111,1010110010101011111011",
	"101011001001011101101111,1001100001100111101111011",
	"100101101101111011010111,1011010001011110101111011");
my @Alpha = (20, 18, 12, 10);
my @GnmAlpha = (10, 11, 12,13);
my @Nbitpat = (1, 4, 4, 5);
my $reduced_alphabet;

while ($_ = $ARGV[0], /^-/) {
	shift;
	if (/^-D/) {$Debug = 1;}
	elsif (/^-R/)	{$reduced_alphabet = 1;}
	elsif (/^-W(\S+)/)	{$out = $1;}
	elsif (/^-K(\S+)/)	{$mol = uc($1);}
	elsif (/^-v(\S+)/)	{$ver = $1;}
	elsif (/^-b(\d+\S*)/ || /^-Xb(\d+\S*)/) {$blksz = $1;}
	elsif (/^-G(\d+\S*)/ || /^-XG(\d+\S*)/)	{$maxgene = $1;}
	elsif (/^-k(\d+)/ || /^-Xk(\d+)/)	{$ktuple = $1;}
	elsif (/^-A(\d+)/ || /^-XA(\d+)/)	{$xa = $1;}
	elsif (/^-C(\d+)/ || /^-XC(\d+)/)	{$xc = $1;}
	else 			{$opt .= "$_ ";}
}

if ($out) {
	my @a = split(/\./, $out);
	$ext = pop(@a);
}
if (!$ext && $mol) {
	my %exts = ('D', 'bkn', 'P', 'bkp', 'A', 'bka');
	$ext = $exts{$mol};
}
&usage unless ($ext);

sub kmg {
	my $val = uc(shift);
	if ($val =~ /(\d+)K/) {return ($1 * 1024);}
	if ($val =~ /(\d+)M/) {return ($1 * 1024 * 1024);}
	if ($val =~ /(\d+)G/) {return ($1 * 1024 * 1024 * 1024);}
	return (int($val));
}

sub estimate {
	my ($gsz, $seqs) = @_;
	my $lout = $out;
	unless ($lout) {
	    my @a = split('.', $seqs);
	    pop(@a);
	    $lout = join('.', @a) . '.'. $ext;
	}
	$gsz *= 2 if ($ext eq 'bkn' || $ext eq 'bkp');
	if ($blksz) {
	    $blksz = &kmg($blksz);
	} else {
	    $blksz = sqrt($gsz);
	    $blksz = int($blksz / 1024 + 1) * 1024;
	}
	if ($maxgene) {
	    $maxgene = &kmg($maxgene);
	} else {
	    $maxgene = $genlencoeff * sqrt($gsz);
	    $maxgene = int($maxgene / 1024 + 1) * 1024;
	    $maxgene = $MinMaxGene if ($maxgene < $MinMaxGene);
	}
	my $maxblk = $maxgene / $blksz;
	$maxblk = int($maxblk / 16 + 1) * 16;
	if ($ext eq 'bkn') {
	    my $xopt = '-KD';
	    $ktuple = int(0.59 * log($gsz)) unless ($ktuple);
	    $ktuple = 16 if ($ktuple > 16);
	    $ktuple = 3 if ($ktuple < 3);
	    my $km3 = $ktuple - 3;
	    my $sht = $km3? $ktuple: 2;
	    if ($xc) {
		$xb = $Bitpat[$km3] unless ($xb);
		$xopt .= " -XC$xc -XB$xb";
	    } elsif ($xb) {
	        $xopt .= " -XB$xb";
	    }
	    $opt = "$xopt -Xs$sht " . $opt;
	} elsif ($ext eq 'bkp') {
	    my $xopt = '-KP';
	    $ktuple = int(0.35 * log($gsz)) unless ($ktuple);
	    $ktuple = 6 if ($ktuple > 6);
	    $ktuple = 3 if ($ktuple < 3);
	    ++$ktuple if ($reduced_alphabet);
	    my $km3 = $ktuple - 3;
	    my $sht = $km3? $ktuple: 2;
	    if ($xc) {
		$xb = $Bitpat[$km3] unless ($xb);
		$xa = $reduced_alphabet? $GnmAlpha[$km3]: 20 unless ($xa);
		$xopt .= " -XA$xa -XB$xb -XC$xc";
	    } else {
		$xopt .= " -XA$xa" if ($xa);
		$xopt .= " -XB$xb" if ($xb);
	    }
	    $opt = "$xopt -Xs$sht " . $opt;
	} elsif ($ext eq 'bka') {
	    $opt = '-KA ' . $opt;
	    $ktuple = int(0.3 * log($gsz)) unless ($ktuple);
	    $ktuple = 6 if ($ktuple > 6);
	    $ktuple = 3 if ($ktuple < 3);
	    my $km3 = $ktuple - 3;
	    my $sht = $km3? $ktuple: 2;
	    $xa = $Alpha[$km3] unless ($xa);
	    $xb = $Bitpat[$km3] unless ($xb);
	    $xc = $Nbitpat[$km3] unless ($xc);
	    $opt = "-XC$xc -XA$xa -XB$xb -Xs$sht " . $opt;
	} else {
	    die "Extension must be .bka, .bkn, or .bkp !\n";
	}
	my $cmd = "spaln";
	$cmd .= $ver if ($ver);
	$cmd .= " -W$lout -Xk$ktuple";
	if ($ext eq 'bkn' || $ext eq 'bkp') {
	    $cmd .= " -Xb$blksz -XG$maxgene -Xg$maxblk";
	}
	$cmd .= " $opt$seqs";
	print $cmd, "\n";
	exit(1) if (!$Debug && system($cmd));
}

sub fromgrp {
	my $gfn = shift;
	my $gsz = 0;
	my $seqs;
	open (GRP, $gfn) or die "Can't open $gfn !\n";
	while (<GRP>) {
	    my @a = split;
	    if ($a[2] eq "E_O_F") {
		$gsz = $a[0];
	    } else {
		$seqs .= ' ' . $a[2];
	    }
	}
	&estimate($gsz, $seqs);
}

while (my $g = $ARGV[0]) {
	shift;
	my @a = split(/\./, $g);
	my $b = $g;
	my $e;
	if (@a > 1) {
	    $e = pop(@a);
	    $b = join('\.', @a);
	}
	my $gfn = $b . ".grp";
	if ($e eq 'grp') {&fromgrp($gfn);}
	elsif (-s $gfn) {
	    my ($gsz, $chr, $dmy) = split(' ', `tail -1 $gfn`);
	    &estimate($gsz, $g, 1);
	} else {
	    printf STDERR "$gfn not found !\n";
	    next;
	}
}
