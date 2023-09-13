/*****************************************************************************
*
*	Execution of a program with a series of data
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-2023)
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

#include "calcserv.h"

InFiles::InFiles(const char* fname, int ac, const char** arv, 
		int imode, int cmode)
	: name(0), grp2(0), idx_g2(0)
{
	int	ctgnmbr = 0;	// members in catalog
	int	argsize = 0;	// total lengths of argx
	int	ctgsize = 0;	// size of catalog file
	const	char**	av = arv;
	for (int i = 0; i++ < ac; ++av) argsize += strlen(*av) + 1;
	FILE*	fd = fname? fopen(fname, "r"): 0;
	if (fd) {
	    char	str[MAXL];
	    while (fgets(str, MAXL, fd)) {
		if (*str == _LCOMM) flush_line(fd);
		else	++ctgnmbr;
	    }
	    ctgsize = ftell(fd);
	}
	membersize = ac + ctgnmbr;
	if (!membersize) return;	// no input seq
	char**	ptr = name = new char*[membersize + 2];
	char*	wk = *name = new char[argsize + ctgsize + 2];
	if (fd) {			// read catalog
	    rewind(fd);
	    if (fread(*name, ctgsize, 1, fd) != 1)
		fatal("Can't read %s !\n", fname);
	    if (wk[ctgsize - 1] != '\n') wk[ctgsize] = '\n';
	    for (int i = 0; i++ < ctgnmbr; ++ptr) {
		while (*wk == _LCOMM)
		    while (*wk++ != '\n') ; // skip comment
		for (*ptr = wk; *wk != '\n'; ++wk) ;
		*wk++ = '\0';
       		if (!**ptr) grp2 = ptr + 1;	// blank line
		else if (**ptr == HARD_DELIM) {
		    **ptr = '\0';
		    grp2 = ptr + 1;
		    break;		// hard delimiter
		}
	    }
	    fclose(fd);
	}
	av = arv;
	for (int i = 0; i++ < ac; ) {	// argments
	    if (!grp2 && ptr > name && imode != IM_SNGL) {
		*ptr = wk - 1;		// delimiter
		grp2 = ++ptr;		// 2nd group members
	    }
	    strcpy(wk, *av);
	    *ptr++ = wk;
	    wk += strlen(*av++) + 1;
	}
	*ptr = 0;
	if (cmode == IM_GRUP) {
	    if (!(grp2 && *grp2)) fatal("2nd group is missing !\n");
	    idx_g2 = grp2 - name - 1;
	} else	grp2 = 0;
}

InFiles::~InFiles()
{
	if (name) {
	    delete[] *name;
	    delete[] name;
	};
}

