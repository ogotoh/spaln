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

INT	min_samples = 0;

int Ild_data::from_eij(const char* fname)
{
	ReadFile	fp(fname);
	if (fp.gzfd)	read_ild_for_map(fp.gzfd);
	else if (fp.fd)	read_ild_for_map(fp.fd);
	else		fatal(not_found, fname);
	return (header.ndsize);
}

int Ild_data::from_file(const char* fname)
{
	ReadFile	fp(fname);
	int	read = 0;
	if (fp.gzfd) {
	    if (fp.dtype == 1)	read = read_text(fp.gzfd);
	    if (fp.dtype == 2)	read = read_binary(fp.gzfd, fname);
	} else if (fp.fd) {
	    if (fp.dtype == 1)	read = read_text(fp.fd);
	    if (fp.dtype == 2)	read = read_binary(fp.fd, fname);
	}
	return (read);
}

void Ild_data::to_file(const char* oname, const int& text, const bool& gzip)
{
	WriteFile	fp(oname, text, gzip);
	if (text) {
	    if (fp.gzfd) {		// gzipped text output
		write_text(fp.gzfd);
	    } else if (fp.fd) {	// plain text output
		write_text(fp.fd);
	    } else {
		fatal(no_file, oname);
	    }
	} else {
	    if (fp.gzfd) {		// gzipped binary output
		write_binary(fp.gzfd, oname);
	    } else if (fp.fd) {	// plain binary output
		write_binary(fp.fd, oname);
	    } else {
		fatal(no_file, oname);
	    }
	}
}
