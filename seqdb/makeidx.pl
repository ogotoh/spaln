#!/usr/bin/perl

# Usage: makeidx.pl [-ianp] [-ggz_dir] [-k] [fasta_file | gzipped fasta_file]
# -i: Make index
# -a: Make block info of aa sequences
# -n: Make block info of genomic sequence for DNA queries
# -p: Make block info of genomic sequence for aa queries

# Example: foreach.pl makeidx.pl `cat fasta.lst`

use strict;

my ($obj, $keep, $gz, $opt) = ('i');
my $makefile = "./Makefile";
my $gnmgfgz = "./";

while ($_ = $ARGV[0], /^-/) {
	shift;
	if (/^-k/) {$keep = 1;}			# don't delete uncompressed file
	elsif (/^-g(\S+)/) {$gnmgfgz = $1;}	# define gz_dir
	elsif (/^(-y\S+)/ || /^(-X\S+)/ || /^(-R\S*)/ || /^(-D\S*)/) {$opt .= "$1 ";}
	elsif (/^-(\w+)/)  {$obj .= lc($1);}	# [a|n|p|np]
}

my $fsrc = shift;
unless (-s $fsrc) {
	$gnmgfgz .= '/' if (substr($gnmgfgz, -1, 1) ne '/');
	$fsrc = "$gnmgfgz$fsrc";
	exit(1) unless (-s $fsrc);
}
my $sl = rindex($fsrc, '/');
my $fa = substr($fsrc, $sl + 1);

### examine free disk space

open(DF, "df |") or die "Can't run df !\n";
while (<DF>) {
	chomp;
	my @a = split;
	my $mp = pop(@a);
	next if ($mp ne "/data");
	my $pc = pop(@a);
	chop($pc);
	exit(1) if ($pc >= 99);
}
close(DF);

print STDERR "$fa $obj\n";
my $fnm = $fa;
my $bdy = $fa;
my @ind = ('i', 'a', 'n', 'p');
my @ext = ('idx', 'bka', 'bkn', 'bkp');

my $dot = rindex($fa, '.');
if ($dot > 0) {
	$bdy = $fnm = substr($fa, 0, $dot);
	if (($gz = substr($fa, $dot + 1) eq "gz")) {
	    $dot = rindex($fnm, '.');
	    $bdy = ($dot > 0)? substr($fnm, 0, $dot): $fnm;
	    my $exp = 0;
	    for (my $i = 0; $i < @ind; ++$i) {
		next if (index($obj, $ind[$i]) < 0);
		my $ffn = $bdy . '.' . $ext[$i];
		++$exp unless (-s $ffn && -M $ffn < -M $fa);
	    }
	    exit(0) unless ($exp); 
	    system("gzip -cd $fsrc > $fnm") unless (-s $fnm);
	}
}

for (my $i = 0; $i < @ind; ++$i) {
    if (index($obj, $ind[$i]) >= 0) {
	my $ffn = $bdy . '.' . $ext[$i];
	my $odr = $bdy . '.odr';
	unlink($odr) if ($i == 0 && -e $odr);
	my $cmd = "make -f $makefile ";
	$cmd .= "OPT='$opt' " if ($opt);
	$cmd .= $ffn;
	system($cmd);
    }
}

unlink $fnm if ($gz && !$keep);
