##############################################################################
#
#	Extract information about exon-intron structure from a extended
#	fasta protein/DNA/RNA sequence file
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

package	Util;
require	Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(getoptarg, System);

sub getoptarg {
	my ($opt, $arg, $dflt) = @_;
	return ($$opt = $arg) if ($arg eq '0');
	my $nextarg = $ARGV[0];
	if (!$arg && defined($nextarg) && !($nextarg =~ /^-/)) {
	    $arg = $nextarg;    # set arg from next arguement
	    shift(@ARGV);
	}
	if ($arg) {$$opt = $arg;}
	elsif ($dflt) {$$opt = $dflt;}
	return ($$opt);
}

sub System {
	my $cmd = shift;
	my $verbs = shift;
	my $debug = shift;
	print STDERR "$cmd\n" if ($verbs > 0);
	return system($cmd) unless ($debug);
}

1;

