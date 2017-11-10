#	Makefile for seqdb
#	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
#	last modified 2014/12/12

.SUFFIXES: .idx .faa .fn .cf .gf .mfa .grp .bka .bkn .bkp .spn .srd

# for this example, the following four commands will
# generate cDNA-genome mapping in the simplified 
# Nagasaki format
#	% make dictdisc.cf
# 	% makeidx.pl -in dictdisc_g.gf
# 	% make dictdisc.srd

clear:
	rm -f *.ent *.idx *.grp *.odr *.seq *.bka *.bkn *.bkp 
	rm -f *.srd *.err *.grd *.erd *.qrd
	rm -f *.gf *.cf *.faa *.spn

.cf.idx:
	mv $< $*
	makdbs -KD $*
	mv $* $<

.fn.idx:
	mv $< $*
	makdbs -KD $*
	mv $* $<

.gf.idx:
	mv $< $*
	makdbs -KD $*
	mv $* $<

.gf.bkp:
	./makblk.pl -W$@ $(OPT) $<

.gf.bkn:
	./makblk.pl -W$@ $(OPT) $<

.faa.idx:
	mv $< $*
	makdbs -KP $*
	mv $* $<

.faa.bka:
	./makblk.pl -W$@ $(OPT) $<

.mfa.idx:
	mv $< $*
	makdbs -KD $*
	mv $* $<

.mfa.bkn:
	./makblk.pl -W$@ $(OPT) $<

.mfa.bkp:
	./makblk.pl -W$@ $(OPT) $<

.grp.bka:
	./makblk.pl -W$@ $(OPT) $<

.grp.bkn:
	./makblk.pl -W$@ -G $(OPT) $<

.grp.bkp:
	./makblk.pl -W$@ -G $(OPT) $<

.faa.spn:
	(spaln -Q7 -O7 -d$*_g -T$* $< > $@) 2> $*spp.err

.cf.srd:
	(time spaln -Q7 -O12 -LS -M -pq -d$*_g $<) 2> $*.err
	sortgrcd $* > $@

dictdisc_g.gf:	dictdisc_g.gf.gz
	gzip -cd dictdisc_g.gf.gz > $@

dictdisc.faa:	 dictdisc.faa.gz
	gzip -cd  dictdisc.faa.gz > $@

dictdisc.cf:	 dictdisc.cf.gz
	gzip -cd  dictdisc.cf.gz > $@
