/*****************************************************************************
*
*	header for pairwise calculation of distance
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

#ifndef  _DIST2_
#define  _DIST2_

#include "cmn.h"

enum DistMethod {EC, JS, JA, KL, MH, UA, E2, CS, XS};
enum OutMode {DMX, PAIR};

const	float	def_scale = 100.;

template <class var_t>
class Dist2 {
template <typename file_t>
	void	write_text(file_t fd);
template <typename file_t>
	void	write_binary(file_t fd, const char* fname);
template <typename file_t>
	void	read_text(file_t fd);
template <typename file_t>
	void	read_binary(file_t fd, const char* fname);
	void	(*optjob)(var_t* pvar) = 0;
public:
	INT	mem_num = 0;
	INT	dist_size = 0;
	INT	out_mode = 0;
	int	calc_mode = 0;
	int	grp2 = 0;
	var_t**	vars = 0;
	Strlist*	sname = 0;
	int*	nmidx = 0;
	int	inodr = 0;
	FTYPE*	dist = 0;
	FTYPE*	getdist() {FTYPE* rv = 0; swap(rv, dist); return (rv);}
	void	from_file(const char* fname);
	void	to_file(const char* fname, 
	    const int& text = 1, const bool& gzip = false);
	Dist2<var_t>(int argc, const char** argv, 
	    const int& calc_mode = IM_EVRY, void (*oj)(var_t* var) = 0);
	Dist2<var_t>(const char* fname, const int& calc_mode);
	Dist2<var_t>(int n = 0) : mem_num(n) {
	    sname = new Strlist;
	}
	Dist2<var_t>(const char* fname) {
	    from_file(fname);
	}
	~Dist2<var_t>() {;
	    if (vars)
		for (INT n = 0; n < mem_num; ++n) delete vars[n];
	    delete[] vars; delete sname;
	    delete[] nmidx; delete[] dist;
	}
};

template <class var_t>
template <typename file_t>
void Dist2<var_t>::write_binary(file_t fd, const char* fname)
{
	sname->write_binary(fd);
	if (fwrite(&dist_size, sizeof(INT), 1, fd) != 1)
	    fatal(write_error, fname);
	if (fwrite(dist, sizeof(float), dist_size, fd) != dist_size)
	    fatal(write_error, fname);
}

template <class var_t>
template <typename file_t>
void Dist2<var_t>::write_text(file_t fd)
{
	char	str[MAXL];
	for (INT i = 0; i < mem_num; ++i) {
	    if (sname && (*sname)[i] && *(*sname)[i])
		fputs((*sname)[i], fd);
	    else {
		snprintf(str, MAXL, "OTU%d", i);
		fputs(str, fd);
	    }
	    fputc('\n', fd);
	}
	fputc('\n', fd);
	for (INT j = 1, k = 0; j < mem_num; ++j) {
	    int	n = 0;
	    for (INT i = 0; i < j; ++i, ++k) {
		snprintf(str, MAXL, " %15.7e", dist[k]);
		fputs(str, fd);
		if (++n == 5) {
		    fputc('\n', fd);
		    n = 0;
		}
	    }
	    if (n) fputc('\n', fd);
	}
}

template <class var_t>
template <typename file_t>
void Dist2<var_t>::read_binary(file_t fd, const char* fname)
{
	sname = new Strlist(fd);
	if (fread(&dist_size, sizeof(INT), 1, fd) != 1)
	    fatal(read_error, fname);
	dist = new FTYPE[dist_size];
	if (fread(dist, sizeof(FTYPE), dist_size, fd) != dist_size)
	    fatal(read_error, fname);
	mem_num = sname->size();
}

template <class var_t>
template <typename file_t>
void Dist2<var_t>::read_text(file_t fd)
{
	char	str[MAXL];
	sname = new Strlist;
	while (fgets(str, MAXL, fd)) {
	    if (*str == '\n') break;
	    str[strlen(str) - 1] = '\0';
	    sname->push(str);
	}
	mem_num = sname->size();
	dist_size = elem(mem_num, 0);
	dist = new FTYPE[dist_size];
	FTYPE*	wk = dist;
	FTYPE*	dt = dist + dist_size;
	while (wk < dt && fgetw(str, MAXL, fd))
	    *wk++ = atof(str);
}

template <class var_t>
void Dist2<var_t>::to_file(const char* fname, const int& text, const bool& gzip)
{
	WriteFile	fp(fname, text, gzip);
	if (fp.fd) {
	    if (text)	write_text(fp.fd);
	    else	write_binary(fp.fd, fname);
	} else if (fp.gzfd) {
	    if (text)	write_text(fp.gzfd);
	    else	write_binary(fp.gzfd, fname);
	} else {
	    fatal(no_file, fname);
	}
}

template <class var_t>
void Dist2<var_t>::from_file(const char* fname)
{
	ReadFile	fp(fname);
	if (fp.fd) {
	    if (fp.dtype == 1)	read_text(fp.fd); else
	    if (fp.dtype == 2)	read_binary(fp.fd, fname);
	    else	fatal(read_error, fname);
	} else if (fp.gzfd) {
	    if (fp.dtype == 1)	read_text(fp.gzfd); else
	    if (fp.dtype == 2)	read_binary(fp.gzfd, fname);
	    else	fatal(read_error, fname);
	} else {
	    fatal(not_found, fname);
	}
}

template <class var_t>
Dist2<var_t>::Dist2(int argc, const char** argv, const int& calc_mode, 
	void (*oj)(var_t* var)) : 
	optjob(oj), mem_num(0)
{
	vars = new var_t*[argc];
	sname = new Strlist;
	char	str[MAXL];
	for (int i = 0; i < argc; ++i, ++argv) {
	    vars[mem_num] = new var_t(*argv, mem_num);
	    if (vars[mem_num]->good()) ++mem_num;
	    else {
		delete vars[mem_num];
		continue;
	    }
	    strcpy(str, *argv);
	    char*	sls = strrchr(str, '/');
	    sls = sls? sls + 1: str;
	    char*	dot = strchr(sls, '.');
	    if (dot) *dot = '\0';
	    sname->push(sls);
	    if (optjob) optjob(vars[mem_num]);
	}
	switch (calc_mode) {
	    case IM_ALTR: case IM_PARA:
		dist_size = mem_num / 2;
		break;
	    case IM_EVRY:
		dist_size = elem(mem_num, 0);
		break;
	    default:
		dist_size = mem_num;
		break;
	}
	dist = new FTYPE[dist_size];
}

template <class var_t>
Dist2<var_t>::Dist2(const char* fname, const int& calc_mode)
{
	sname = new Strlist;
	Mfile	mfd(sizeof(var_t*));
	FILE*	fd = fopen(fname, "r");
	while (!feof(fd)) {
	    var_t*	var = new var_t;
	    var->fget(fd);
	    sname->push(var->entry());
	    mfd.write(&var);
	}
	fclose(fd);
	mem_num = mfd.size();
	vars = (var_t**) mfd.flush();
	switch (calc_mode) {
	    case IM_ALTR: case IM_PARA:
		dist_size = mem_num / 2;
		break;
	    case IM_EVRY:
		dist_size = elem(mem_num, 0);
		break;
	    default:
		dist_size = mem_num;
		break;
	}
	dist = new FTYPE[dist_size];
}

#endif	// _DIST2_

