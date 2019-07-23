# SPALN information

### Map and align a set of cDNA/EST or protein sequences onto a genome
#### Present Version 2.3.3e
#### Last updated: 2019-07-23

- [Overview](#Ov)
- [Install](#Inst)
- [Format](#Format)
- [Execution](#Exec)
- [Example](#Exam)
- [References](#Ref)
- [Changes](#Changes)

* * *

## <a name="Ov">Overview</a>

**Spaln** (space-efficient spliced alignment) is a
stand-alone program that maps and aligns a set of cDNA or protein sequences
onto a whole genomic sequence in a single job. **Spaln** also performs spliced
or ordinary alignment after rapid similarity search against a protein sequence
 database, if a genomic segment or an amino acid sequence is given as a query.
From Version 1.4, spaln supports a combination of protein sequence database and
a given genomic segment. From Version 2.2, spaln also performs rapid similarity
search and (semi-)global alignment of a set of protein sequence queries against
a protein sequence database. **Spaln** adopts multi-phase
heuristics that makes it possible to perform the job on a conventional personal
computer running under Unix/Linux and MacOS with limited memory. The program is
written in C++ and distributed as source codes and also as executables for a few
platforms. Unless binaries are not provided, users must compile the program on
their own system. Although the program has been tested only on a Linux operating
system, it is likely to be portable to most Unix systems with little or no
modifications. The accessory program **sortgrcd** sorts the gene loci found
by **spaln** in the order of chromosomal position and orientation. From version
2.3.2, **spaln** and **sortgrcd** can handle gzipped genome/database files and
'block' files without prior expansion if USE_ZLIB mode is activated upon
compilation. From version 2.3.2a, compressed query sequence file(s) may also be
accepted.

## <a name="Inst">Install</a>

To compile the source codes in the default settings, follow the instructions below. When you download the source file in the directory _download_, five directories will be generated under _download/spalnXX/_ after installation, where XX is a version code. We assume _work_ is your workspace, which may or may not be identical to _download_

 * bin : binaries
 * doc : documents
 * seqdb : sample sequences. In this directory you should format genomic or database files
 * src : source codes
 * table : parameter files used by **spaln**

To modify the location of executables and/or other settings, run 'configure --help' at step 6 below. (**Warning**: Full path name rather than relative path name must be given for executables or other directories as the arguments of the **configure** command.) These locations are hard coded in **spaln**. The locations of the 'seqdb' and 'table' directories will be respectively denoted by _seqdb_ and _table_ below. Hence, _seqdb_=_download/spalnXX/seqdb_, and _table_=_download/spalnXX/table_ in the default settings.

1. `% mkdir download`
2. `% cd download`
3. Download spalnXX.tar.gz
4. `% tar xfz spalnXX.tar.gz`
5. `% cd ./spalnXX/src`
6. <a name="compile">For compilation</a>  
    `% ./configure [--help]`  
    Please manually edit Makefile if $(CXX) does not indicate a C++ compiler or  
    `% CXX=g++ ./configure [other options]`
    `% ./configure [other options] --use_zlib=1` to use zlib facilities.  
7. `% make`
8. `% make install`
 Executables are copied to ../bin  
 makmdm program makes mutation data matrices of various PAM levels in the ../table directory
9. `% make clearall`
10. Add _download/spalnXX/bin_ to your PATH
    % setenv PATH $PATH:download/spalnXX/bin (csh/tsh)
    $ export PATH=$PATH:download/spalnXX/bin (sh/bsh)
 Preferably, you may add the above line in your start up rc file (e.g. ~/.bashrc)  
   Alternatively, move or copy _download/spalnXX/bin/\*_ to a directory on your PATH, if you have not specified the location of executables at step 6 above.
11. If you have changed the location of _table_ and/or _seqdb_ directory after installation, set the env variables ALN_TAB and/or ALN_DBS as follws:
   * % setenv ALN_TAB New_Aln_Tab (csh/tsh)
   * $ export ALN_TAB=New_Aln_Tab (sh/bsh)
   * % setenv ALN_DBS New_Aln_Dbs (csh/tsh)
   * $ export ALN_DBS=New_Aln_Dbs (sh/bsh)  
Add the above lines to your rc file, so that you don't have to repeat the commands at every login time.
12. Proceed to [Format](#Format).

## <a name="Format">Format</a>

If you do not need genome mapping or database search, you may skip this section. All sequence files should be in (multi-)fasta format.

To perform genome mapping, the genomic sequence must be formatted before use. Formatting is optional for amino acid sequence database search.

1. `% cd seqdb`
2. Download or copy genomic sequences or protein database sequence in multi-fasta format. If **spaln** is [compiled](#compile) accordingly, gzipped file need not be uncompressed (the file name must be _X_.gz).
3. Chromosomal sequences should be concatenated into a single file. Alternatively, you can use multiple chromosomal files without concatenation. This procedure will be described at the end of this section. To render the 'make' command effective, the extension of the genomic sequence file should be '.mfa' or '.gf', and protein database sequence should be '.faa'. Hereafter, the file name is assumed to be xxxgnm.mfa or prosdb.faa.
4. `% ./makeidx.pl -i[n|p|np] [-g] [spaln options] xxxgnm.mfa(.gz)` or  
   `% ./makeidx.pl -i[a] [-g] [spaln options] prosdb.faa(.gz)`  
 These commands are shortcuts that replace the following series of operations 
5-7, if the input is a single sequence file. In that case, you can skip following 
instructions. The block size and *k*-mer size are estimated from the 
genome size. The -ix option specifies the "block file(s)" .bkx to be constructed, 
where *x* is 'a', 'n' or 'p'. The -inp option will construct both .bkn and .bkp 
files together with the .idx and associated files. If -ix is omitted or *x* is 
empty, no block file is constructed. The -g option specifies gzipped output.
5. `% make xxxgnm.idx` (for genomic sequence) or  
   `% make prosdb.idx` (for protein database sequence)  
  This command converts the sequence into a binary format. Four or five files, xxxgnm.seq, xxxgnm.idx, xxxgnm.ent, xxxgnm.grp, and optionally xxxgnm.odr are constructed (prosdb instead of xxxgnm in case of make prosdb.idx). It may take several tens of minutes to construct the files for mammalian genome.
6. `% make xxxgnm.bkn` (for cDNA queries) or  
   `% make xxxgnm.bkp` (for protein queries) or  
   `% make prosdb.bka` (for protein database)  
 * This command makes the block index table. This process may take another several tens of minutes.  
 * Internally, the make command invokes  
    `spaln -Wxxxgnm.bkn -KD [Options] xxxgnm.mfa` or  
    `spaln -Wxxxgnm.bkp -KP [Options] xxxgnm.mfa`  or  
    `spaln -Wprosdb.bka -KA [Options] prosdb.faa`  
 * If xxxgnm.grp or prosdb.grp were successfully constructed at step 5 above, the option values below would be automatically calculated by script makblk.pl. **WARNING:** The estimated maximal gene size can be inadequately small if only a part of the genome (*e.g* a single chromosome) is formatted. At that time, explicitly specify the maximal gene size by the -XG*N* option of makblk.pl or at the runtime of **spaln**. *N* can have suffix 'k' and 'M' to indicate that the number is measured in kbp and Mbp, respectively.  
 * Options: (default value)
  * -XA *N*: alphabet size of the reduced amino acids: 6 < *N* <= 20 (20)
  * -XB *S*: bit patterns of the spaced seeds. The pattern should be asymmetric  when the number of patterns > 2.
  * -XC *N*: number of seed patterns: 0 <= *N* <= 5 (0: contiguous seed)
  * -XG *N*: maximum gene length (262144)
  * -Xa *N*: a parameter used to filter excessively abundant words (10)
  * -Xb *N*: block size (4096) An estimate of *N* is sqrt(genome size). For mammals, *N* is nearly equal 54000.
  * -Xg *N*: maximal distance in block number between 5' terminal and 3' terminal blocks (16)
  * -Xk *N*: word size (11 for DNA, 5 for protein)
  * -Xs *N*: distance between neighboring seeds (= *k*)
7. It is possible to generate xxxgnm.idx and other three files directly from the input files without concatenation:  
    `% makdbs -nxxxgnm -KD file1 ... fileN` and  
    `% make xxxgnm.bkn` (for cDNA queries) or  
    `% make xxxgnm.bkp` (for protein queries)  
 This method is particularly useful when the concatenation might yield a file too large to be dealt with by the OS.

## <a name="Exec">Execution</a>

1. Prepare protein, cDNA, or genomic segment sequence(s) in (multi-)fasta format
(denoted by *query* below). From 2.3.2a, zgipped fasta file(s) may be used as 
the query without prior expansion. <u>Note, however, that compressed query can considerably slow down the execution rate.</u>
2. Store *query* to _work_.
3. `% cd work`
4. Run **spaln** in one of the following four modes. **Spaln**
    does not support comparison between two genomic segments.  
```
    % spaln -Q[0|1|2|3] [-ON] [other options] genome_segment query
    % spaln -Q[4|5|6|7] [-ON] [other options] -d xxxgnm query
    % spaln -Q[4|5|6|7] [-ON] [other options] -a prosdb query
    % spaln -Q[4|5|6|7] [-ON] [other options] prosdb.faa query
```
   * In the last case, *prosdb.faa* will be internally formatted, and the formatted results will be discarded after the end of execution.
   * Only a subset of queries may be examined if *query* is replaced with '*query* (from to)', where 'from' and 'to' are the first and last entry numbers in *query* to be examined.  
   * Options: (default value)
     * -C *N*:	Use the genetic code specified by the "transl_table number" defined in [NCBI transl_table](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) (1).
     * -H *N*: Output is suppressed if the alignment score is less than *N*. See also -pw. (35)
     * -K[D|P|A]: Format either genomic DNA for sequence search with DNA (D) or Protein (P) queries or protein database sequences for search with a genomic segment (A). Use in combination with -W option below.
     * -LS: Smith-Waterman-type local alignment. This option may prune out weakly matched terminal regions.
     * -M[*N*]: Single or multiple output for each query
       * No option (default): single locus
       * No argument: Multiple loci up to the maximum number specified by the program (4 in the present implementation)
       * *N*=1: Re-search the *query* region not aligned in the first trial. May be useful to detect chimera or fragmented genomic region
       * *N*\>1: Output multiple loci maximally up to *N*
     * -O *N*: Select output format for genome vs cDNA or aa (4)
        * *N*=0: Gff3 gene format
        * *N*=1: Alignment
        * *N*=2: Gff3 match format
        * *N*=3: Bed format
        * *N*=4: Exon-oriented format similar to the output of megablast -D 3
        * *N*=5: Intron-oriented output
        * *N*=6: Concatenated exon sequence
        * *N*=7: Translated amino-acid sequence. Presently not very useful for cDNA queries because the entire exon rather than an ORF is translated
        * *N*=8: Mapping (block) information only. Use with -Q4
        * *N*=10: SAM format
        * *N*=12: Output the same information as -O4 in the
	    binary format. If -oOutput is set, three files named
	    Output.grd, Output.erd, and Output.qrd will be created. Otherwise,
	    query.grd, query.erd, and query.qrd will be created. If the -g 
            option is set, gzip-compressed outputs will be generated.
     * -O *N*: Select output format for aa vs aa (4)
        * *N*=0: statistics (%divergence alignment_score #match, #mismatch #gap #unpaired)
        * *N*=1: Alignment
        * *N*=2: Sugar format
        * *N*=3: Psl format
        * *N*=4: XYL = Coordinate + match length
        * *N*=5: statistics + XYL
        * *N*=8: Cigar format
        * *N*=9: Vulgar format
        * *N*=10: SAM format
     * -Q *N*: Select algorithm (3)
        * 0<=*N*<=3: Genomic segment in the fasta format given by the first argument vs. *query* given by the second argument. One may skip the formatting step described above if only this mode of operation is used.
        * 4<=*N*<=7: Genome mapping and alignment. The genomic sequence must be formatted beforehand.
        * *N*=0,4: DP procedure without HSP search. Considerably slow
        * *N*=1,2,3,5,6,7: Recursive HSP search up to the level of (*N* % 4)
     * -R *S*:	Read block index table from file *S*.
      If omitted, the xxxgnm.bkn, xxxgnm.bkp,
      or prosdb.bka file will be read depending on the type of query. The
      appropriate file is searched for in the current directory, the directory
      specified by the _env_ variable 'ALN_DBS', and the
      'seqdb' directory specified at the compile time in this order.
     * -S *N*:	Orientation of the DNA query sequence (0)
        * *N*=0: The
       orientation is inferred from the phrases (e.g. 5' end) in the header
       line of each entry within a fasta file. If no information is available,
       both orientations are examined, and the result with the better score is
       reported.
        * *N*=1: Forward orientation only
        * *N*=2: Reverse-complement orientation only
        * *N*=3: Examine both orientations
     * -T *xxx*:	Specify the species. For genome vs. DNA comparison, -yS flag should also be set in combination with this option. *xxx* corresponds to the subdirectory in the _table_ directory.
     * -V *N*:	Minimum space to induce Hirschberg's algorithm (16M)
     * -W *S*:	Write block index table to file *S*.
     * -g: gzipped output used in combination with -W or -O12 option.
     * -i[a|p]: Input mode with -Q[0<=N<=3].
        * -ia: Alternative mode; a genomic segment of an odd numbered entry in the input file is aligned with the query of the following entry.
        * -ip: Parallel mode; the i-th entry in the file specified by the first argument is aligned with the i-th entry in the file specified by the second argument.
        *  	default: The genomic segment specified by the first argument is aligned with each entry in the file specified by the second argument.
     * -o *S*:	Destination of output file name (stdout)
     * -pa:	Terminal polyA or polyT sequence is not trimmed.
     * -pi:	Mark exon-intron junctions by color in the alignment (-O1).
     * -pq:	Suppress warning messages sent to *stderr*.
     * -pw:	Report result even if alignment score is below threshold value.
     * -px:	Suppress self-comparisons in the execution mode (C) or (D).
     * -xB *S*: 	Bit pattern of seeds used for HSP search at level 1
     * -xb *S*:	Bit pattern of seeds used for HSP search at level 3
     * -u *N*:	Gap-extension penalty (3, 2, 2)
     * -v *N*:	Gap-opening penalty (8, 6, 9)
     * -ya *N*:	Dinucleotide pairs at the ends of an intron (0)
        * *N*=0: Accept only the canonical pairs (GT..AG,GC..AG,AT..AC)
        * *N*=1: accept also AT..AN
        * *N*=2: allow up to one mismatch from GT..AG
        * *N*=3: accept any pairs. An omission of *N* implies *N* = 3
     * -yi *N*:	Intron penalty (11, 8, 11)
     * -yj *N*:	Incline of long gap penalty (0.6)
     * -yk *N*:	Flex point where the incline of gap penalty changes (7)
     * -yl *N*:	Double affine gap penalty if *N*=3; otherwise affine gap penalty
     * -ym *N*:	Score for a nucleotide match (2, 2)
     * -yn *N*:	Penalty for nucleotide mismatch (6, 2)
     * -yo *N*:	Penalty for an in-frame termination codon (100)
     * -yp *N*:	PAM level used in the alignment (third) phase (150)
     * -yq *N*:	PAM level used in the second phase (50)
     * -yx *N*:	Penalty for a frame shift (100)
     * -yy *N*:	Relative contribution of splicing signal (8)
     * -yz *N*:	Relative contribution of coding potential (2)
     * -yA *N*:	Relative contribution of the translational initiation or termination signal (8)
     * -yB *N*:	Relative contribution of branch point signal (0)
     * -yE *N*:	Minimum exon length (2)
     * -yI *S*:	Intron distribution parameters
     * -yJ *N*:	Relative contribution of the bonus given
      to a conserved intron position
     * -yL *N*:	Minimum intron length (30, 30)
     * -yS:	For
      cDNA queries, use species-specific exon-intron boundary signals. For
      protein queries, invoke the 'salvage' procedure to examine all blocks
      with positive scores.
     * -yS _N_:	*N* specifies the percentile
      contribution of the species-specific splice signal. The other part is
      derived from the universal signal given to the dinucleotides at the ends
      of an intron. An omission of *N* implies *N* = 100.
     * -yX:	For a DNA query, this option sets parameter
      values for cross-species comparison. The actual values are given as the
      second number in the parentheses of the above lines. -yS100 is also
      automatically invoked. Conversely, this option specifies an intra-species
      mode for a protein query, whereby -yS30 is invoked.
     * -yY *N*:	Relative contribution of length-dependent part of intron penalty (8)
     * -yZ *N*:	Relative contribution of oligomer composition within an intron (0)

5. **Sortgrcd**
  * **Sortgrcd** is used to recover the output of **spaln** with -O12 option, to apply some filtering, and also to rearrange the output of multiple **spaln** runs.
  * Run **sortgrcd** as follows:  
      `% sortgrcd [options] xxx.grd(.gz)`
  * Options:
    * -C _N_: Minimum cover rate = % nucleotides in predicted exons / length of *query* (x 3 if query is protein) (0-100)
    * -E _N_: Report only the best (*N*=1) or all (*N*=2) results per gene locus (1)
    * -F _N_: Filter level (*N*=0: no; *N*=1: mild; *N*=2: medium; *N*=3: stringent)
    * -I _N_: Minimum sequence identity (0-100)
    * -H _N_: Minimum alignment score (35)
    * -O _N_: Output mode. Same as that of **spaln** except that *N*=1, 2, and >=8 are not supported. -O15 reports -O5 format information for only unique introns.
    * -S _C_: Sort chromosomes/contigs in the order of *C*=a: alphabetical, b: abundance, c: appearance in genome database, r: reverse order for minus strand
    * -V _N_: Internal memory size used for core sort. If the
	     data size is greater than *N*, the sorting procedure will be
	     done in pieces.
    * -m _N_:	Maximum number of mismatches within 20 bp from the nearest exon-intron boundary
    * -n _N_:	Maximum number of non-canonical (other than GT..AG, GC..AG, AT..AC) intron ends
    * -u _N_:	Maximum number of unpaired (gap) sites within 20 bp from the nearest exon-intron boundary
  * By default, no filter listed above is applied.
  * When the output of **spaln** is separated into several files, the combined
results are subjected to the sorting. Although xxx.grd (or xxx.grd.gz) files are assigned as the
argument, there must be corresponding xxx.erd and xxx.qrd (or xxx.efd.gz and
xxx.qrd.gz) files in the same directory.
  * In the default output format, the gene structure corresponding to each
transcript is delimited by a line starting with '@', whereas each gene locus is
delimited by a line starting with '!' [4]. Two transcripts belong to the same
locus if their corresponding genomic regions overlap by at least one nucleotide
on the same strand.
  * With -O0 option, the outputs follow the instruction of [Gff3](http://www.sequenceontology.org/gff3.shtml) where a gene locus is defined as described above.


## <a name="Exam">Example</a>
  * To experience the flow of procedures with the samples in _seqdb_, type in the following series of commands after moving to _seqdb_.
```
    % make dictdisc.cf
    % make dictdisc.faa
    % make dictdisc_g.gf
    % perl makeidx.pl -inp dictdisc_g.gf
    % make dictdisc.srd
    % make dictdisc.spn
```
  * Alternatively, you may try below if USE_ZLIB is activated..
```
    % perl makeidx.pl -inp [-g] dictdisc_g.gf.gz
    % spaln -Q7 -d dictdisc_g -T dictdisc [-t10] dictdisc.faa.gz
    % spaln -Q7 -d dictdisc_g -yS -T dictdisc -O12 -g [-t10] dictdisc.cf.gz
    % sortgrcd -O15 -F2 dictdisc.grd.gz
```
    
## <a name="Changes">Changes from previous version</a>
1. The heuristic alignment engine has been updated, resulting in marginal but significant improvement in speed and accuracy, 
especially with -Q3/7 option.
2. chachr.pl has been extended to accept Ganbank/DDBJ and EMBL-formatted files in 
addition to FASTA files. Maybe used as a format convertor.
3. Update help and error messages of **Spaln** and **Sortgrcd**.

## Changes in version 2.3.3</a>
1. The maxmal path size has been extented from 255 to 2047 characters.
2. The 'NEVSEL' constant value has been changed to avoid underflow of 2 * NEVSEL.
3. In utilseq.c and .h, a member variable in class PatMat has been moved to a local variable
to recover thread safety.
4. When **Spaln** is run: `% spaln protein genome`, the order of the 1st and 2nd arguemnets is exchanged with a warning message.

## Changes in version 2.3.2a</a>
1. From this version, query fasta file(s) may be compressed.
2. The new option of <b>spaln</b> '-g' directly generates compressed output(s) when used in combination with -W or -O12 option.
3. makeidx.pl and makblk.pl have been modified to accord with gzipped genome/database fasta files.
4. A small bug in makdbs.c has been fixed.

## Changes in version 2.3.2
1. From this version, input genome/database fasta files (X.mfa, X.gf, or X.faa), formatted data files (X.seq,X.bka, X.bkn, and X.bkp) for **spaln**, and X.grd, X.erd, and X.qrd for **sortgrcd** may be gzipped if USE_ZLIB mode is activated upon [compilation](#compile). **Note:** other data files (X.ent, X.grp, X.idx, and X.odr) must not be compressed.
2. A serious bug concerning with multiple queries has been fixed. This has considerably improved mapping sensitivity particularly when -M option is set under single thread operation mode.
3. Fixation of several small bugs and fine tuning of codes further enhanced mapping sensitivity and specificity particularly for short protein queries.
4. -O *N* option of **sortgrcd** has been extended. -O0: Gff3: -O3: BED; -O4: exon-oriented; -O5: intron-oriented; -O6: concatenated exons; -O7: translated amino acid sequence; -O15: unique introns.

## <a name="Ref">References</a>

<a name="Ref1">[[1]](http://nar.oxfordjournals.org/cgi/content/abstract/gkn105?ijkey=N2yLVza41RuShAg&keytype=ref) Gotoh, O.
A space-efficient and accurate method for mapping and aligning cDNA sequences onto genomic sequence. *Nucleic Acids Research* **36** (8) 2630-2638 (2008).

<a name="Ref2">[[2]](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btn460?ijkey=XajuzvyHlcQZoQd&keytype=ref) Gotoh, O.
Direct mapping and alignment of protein sequences onto genomic sequence. *Bioinformatics* **24** (21) 2438-2444 (2008).

<a name="Ref3">[[3]](http://nar.oxfordjournals.org/content/40/20/e161) Iwata, H. and Gotoh, O.
Benchmarking spliced alignment programs including  Spaln2, an extended version of Spaln that incorporates additional species-specific features. *Nucleic Acids Research* **40** (20) e161 (2012)

<a name="Ref4">[[4]](https://academic.oup.com/bioinformatics/article/22/10/1211/236993) Nagasaki, H., Arita, M., Nishizawa, T., Suwa, M., Gotoh, O.
Automated classification of alternative splicing and transcriptional initiation and construction of a visual database of the classified patterns. *Bioinformatics* **22** (10) 1211-1216 (2006).

* * *

Copyright (c) 1997-2019 Osamu Gotoh (o.gotoh@aist.go.jp) All Rights Reserved.
