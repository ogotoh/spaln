#!/usr/bin/perl

#############################################################################
#
#	clade.pl
#
#	make clades from .dvn/.dvp file
#
#	Input	.dvn/.dvp
#	output	< Thr clades
#	Option	-M (don't include Master) -S (Master only)
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
#############################################################################

my $Thr = 100;
my $Mstr = 1;
my $Slav = 1;
my $Uniq = 1;
my $Sim = 0;
my $avoid;
my %member;
my %master;
my %subsid;

while ($_ = $ARGV[0], /^-/) {
	shift;
	last if /^-$/;
	/^-D/ && ($Debug = 1);
	/^-H(\d+)/ && ($Thr = $1);
	/^-A(\S+)/ && ($avoid = $1);
	/^-M/ && ($Mstr = 0);	# Slave only, Don't output Master
	/^-S/ && ($Slav = 0);	# Master only, Don't output Slave
	/^-U/ && ($Uniq = 0);	# Don't output unique members
	/^-R/ && ($Sim = 1);	# Resemblance
}

while (<>) {
	my (@a) = split;
	my ($scr) = shift(@a);
	my ($b) = pop(@a);
	my ($a) = pop(@a);
	if ($Sim && $scr < $Thr || !$Sim && $scr > $Thr) {
	    $member{$a} = 1;
	    $member{$b} = 1;
	    next;
	}
	if ($avoid && $a =~ /$avoid/ && !($b =~ /$avoid/)) {
	    my($tmp) = $a;
	    $a = $b;
	    $b = $tmp;
	}
	unless ($subsid{$b}) {
	    $subsid{$a} = $a unless $subsid{$a};
	    $mstr = $subsid{$b} = $subsid{$a};
	    $master{$mstr} .= $b . " ";
	}
}

sub bymemb {
	($master{$b} =~ tr/ / /) <=> ($master{$a} =~ tr/ / /);
}

foreach $mstr (sort bymemb keys(%master)) {
	print $mstr," " if ($Mstr);
	print $master{$mstr} if ($Slav);
	print "\n";
}

exit(0) unless ($Uniq);

foreach $mem (sort keys(%member)) {
	print $mem,"\n" unless($subsid{$mem});
}

