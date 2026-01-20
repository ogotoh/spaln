/*****************************************************************************
*
*	Generate xxx.ild[.dat] from xxx.eij or
*	Convert between text and binary forms
*
*	Osamu Gotoh, ph.D.      (-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.      (2001-2023)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#include "ild_data.h"

void usage()
{
	fputs("Description:\n", stdout);
	fputs("\tGenerate intron-length distribution (ild) from ", stdout);
	fputs("exon-intron junction data (eij)\n", stdout);
	fputs("\t  or\n\tConvert between text and binary forms of ild\n", stdout);
	fputs("Usage:\n", stdout);
	fputs("\teij2ild [options] xxx.eij[.gz] or\n", stdout);
	fputs("\teij2ild [options] xxx.ild[.gz|.dgz|.dat[.gz]]\n", stdout);
	fputs("Options:\n", stdout);
	fputs("\t-b S:\tS: binary output file name, [.dat|.dgz] may be added\n", stdout);
	fputs("\t-g:\tgzipped output, .gz may be added to the out-file\n", stdout);
	fputs("\t-o S:\tS: output file name, default: stdout\n", stdout);
	fputs("Warning:\n", stdout);
	fputs("\tExtention of input file must be .eij or .ild[.gz|.dat|.dgz]\n", stdout);
	fputs("\tGzipped eij input is accepted but can slow down calculation\n", stdout);
	exit(1);
}

int main(int argc, const char** argv)
{
const	char*	out_file;
	int	text = 2;
	bool	gzip = false;
	char	str[MAXL];
	while (--argc && **++argv == '-') {
	  switch (argv[0][1]) {
	    case 'b':
		out_file = getarg(argc, argv);
		text = 0;
		break; 
	    case 'g': gzip = true; break;
	    case 'h': usage();
	    case 'o': 
		out_file = getarg(argc, argv);
		break;
	  }
	}
	Ild_data	ild_dt;
	int	read = 0;
	while (argc--) {
const	    char*	fn = *argv++;
	    strcpy(str, fn);
	    if (get_ext(str, "eij")) read += ild_dt.from_eij(str);
	    else	read += ild_dt.from_file(str);
	}
	if (!read)	usage();
	ild_dt.to_file(out_file, text, gzip);
	return (0);
}
