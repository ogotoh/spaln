#############################################################################
#
#	subroutines for species specific table
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
###############################################################################

package	SspTab;
require	Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(speccode opengnmdb sspopt gnmdb setgnmdb avrintlen dbentry ftable fseqdb);

my @tabdirs = (".", "$ENV{HOME}/table", $ENV{"ALN_TAB"});
my @sdbdirs = (".", "$ENV{HOME}/seqdb", $ENV{"ALN_DBS"});
my $tabdir;
my %gnm2tab;
my $gnmdbf;
my $cmndbf;
my $gnmidfn = 8;	# number of letters for identifier
my $gnmidxe = "_g";	# extention to indicate genomic seq.
my $qtl3_4 = 100;	# default 3/4 quantile intron length

sub speccode {
	my $gen = ucfirst(shift);
	my $spc = lc(shift);
	my $genus = substr($gen, 0, 4);
	if (my $ul = 4 - length($genus)) {
	    $genus .= '_' x $ul;
	}
	my $spec = substr($spc, 0, 4);
	if (my $ul = 4 - length($spec)) {
	    if ($spec =~ /^sp\./) {$spec = 'spec';}
	    else {$spec .= '_' x $ul;}
	}
	return $genus . $spec;
}

sub opengnmdb {
	$cmndbf = shift;
	return if (%gnm2tab);
	my $Gnm2Tab;
	foreach (@tabdirs) {
	    $Gnm2Tab = "$_/gnm2tab";
	    if (-s $Gnm2Tab) {
		$tabdir = $_;
		last;
	    }
	}
	open(G2T, $Gnm2Tab) or die "Can't open $Gnm2Tab !\n";
	while (<G2T>) {
	    my ($gid, $spc) = split;
	    $gnm2tab{$gid} = $spc;
	}
	close(G2T);
}

sub sspopt {    # Species specific options
	my $gnm = shift;
	my $fst = ($gnm =~ /^\$/)? 1: 0;
	my $spc = lc(substr($gnm, $fst, $gnmidfn));
	$spc = $gnm2tab{$spc};
	return (defined($spc)? "-T$spc ": "");
}

sub ftable {
	my $genspc = shift;
	my $deftab;
	foreach my $td (@tabdirs) {
	    my $dir = "$td/$gnm2tab{$genspc}";
	    $deftab  = $td if (!$deftab && -s "$td/AlnParam");
	    return ($dir) if (-d $dir);
	}
	return ($deftab);
}

sub fseqdb {
	my $fn = shift;
	if ($fn) {
	    foreach my $td (@sdbdirs) {
		my $id = "$td/$fn";
		return ($id) if (-s $id);
	    }
	} else {
	    for (my $i = 1; $i < @sdbdirs; ++$i) {
		my $td = $sdbdirs[$i];
		return "$td/" if (-d $td);
	    }
	}
	return (undef);
}

sub avrintlen {
	my $gnm = shift;
	my $fst = ($gnm =~ /^\$/)? 1: 0;
	my $spc = lc(substr($gnm, $fst, $gnmidfn));
	$spc = $gnm2tab{$spc};
	my $avi = $qtl3_4;
	if (open(AP, "$tabdir/$spc/AlnParam")) {
	    while (<AP>) {
		next unless (/^-yI/);
		(my $dmy, $avi) = split;
		last;
	    }
	    close(AP);
	}
	return ($avi);
}

sub gnmdb {
	my $seq = shift;
	$seq = '$' . $seq if (substr($seq, 0, 1) ne '$');
	return ("-d$gnmdbf ") if ($gnmdbf && `utn -n -pq -d$gnmdbf '$seq 1 1'`);
	my $gdb = lc(substr($seq, 1, $gnmidfn)) . $gnmidxe;
	return ("-d$cmndbf ") if ($cmndbf && `utn -n -pq -d$cmndbf '$seq 1 1'`);
	return ("-d$gdb ") if (`utn -n -pq -d$gdb '$seq 1 1'`);
	return (undef);
}

sub setgnmdb {
	$gnmdbf = shift;
}

sub dbentry {
	my $hline = shift;
	@a = split(' ', $hline);
	if ($a[0] eq '>' || $a[0] eq ';D' || $a[2] eq '+' || $a[2] eq '-')
	    {shift(@a);}
	elsif (/^>/) {$a[0] = substr($a[0], 1);}
	elsif (/^;D/) {$a[0] = substr($a[0], 2);}
	my $st = $a[$#a - 1];
	my $bias = ($st eq 'T' || $st eq 'Y')? 8: 0;
	$st = 'X' unless ($st eq 'N' || $st eq 'Q' || $st eq 'T' || $st eq 'Y');
	my $gene = ($st? '$': "") . $a[$bias];
	return ($gene, $a[$bias + 1], $a[$bias + 4], $a[$bias + 6], $st);
}

