# SPALN information

### Map and align a set of cDNA/EST or protein sequences onto a genome
#### Present Version 2.4.9
#### Last updated: 2022-05-10

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
accepted. From version 2.4.0, multiple files corresponding to different output 
forms can be generated at a single run. 

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
    To use zlib facilities, confirm that libz.* are installed in the load 
    library path.  Then,  
    `% ./configure [other options] --use_zlib=1`  
7. `% make`
8. `% make install`
 Executables are copied to ../bin  
 makmdm program makes mutation data matrices of various PAM levels in the ../table directory
9. `% make clearall`
10. Add _download/spalnXX/bin_ to your PATH  
    `% setenv PATH $PATH:download/spalnXX/bin (csh/tsh)`  
    `$ export PATH=$PATH:download/spalnXX/bin (sh/bsh)`  
 Preferably, you may add the above line in your start up rc file (e.g. ~/.bashrc)  
   Alternatively, move or copy _download/spalnXX/bin/\*_ to a directory on your PATH, if you have not specified the location of executables at step 6 above.
11. If you have changed the location of _table_ and/or _seqdb_ directory after installation, set the env variables ALN_TAB and/or ALN_DBS as follws:
   * `% setenv ALN_TAB New_Aln_Tab (csh/tsh)`
   * `$ export ALN_TAB=New_Aln_Tab (sh/bsh)`
   * `% setenv ALN_DBS New_Aln_Dbs (csh/tsh)`
   * `$ export ALN_DBS=New_Aln_Dbs (sh/bsh)`  
    Add the above lines to your rc file, so that you don't have to repeat the commands at every login time.
12. Proceed to [Format](#Format).

## <a name="Format">Format</a>

If you do not need genome mapping or database search, you may skip this section.
 All sequence files should be in (multi-)fasta format.

To perform genome mapping, the genomic sequence must be formatted before use. 
Formatting is optional for amino acid sequence database search.

1. `% cd seqdb`
2. Download or copy genomic sequences or protein database sequence in multi-fasta format. If **spaln** is 
[compiled](#compile) accordingly, gzipped file need not be uncompressed (the file name must be _X_.gz).
3. To use 'makeidx.pl' command, chromosomal sequences must be concatenated into a single file. The extension of the genomic sequence file must be '.mfa' or '.gf', and protein database sequence must be '.faa', to render 'make' command effective. With 'spaln -W' command, these restrictions are not obligatory. Hereafter, the file name is assumed to be xxxgnm.gf or prosdb.faa. 
4. To format xxxgnm.gf(.gz), run either of the following two commands, which  are equivalent to each other except that the former is faster, accepts multiple input files, and does not need Makefile.  
   `% spaln -W -K[D|P] [-XGMAX_GENE] [spaln options] xxxgnm.gf(.gz) ...`  
   `% ./makeidx.pl -i[n|p|np] [-XGMAX_GENE] [spaln options] xxxgnm.gf(.gz)`  
To format prosdb.faa(.gz), run either of the following two commands, which  are equivalent to each other except that the former accepts multiple input files.  
   `% spaln -W -KA [spaln options] prosdb.faa(.gz) ...`  
   `% ./makeidx.pl -ia [spaln options] prosdb.faa(.gz)`  
 * -K*X* (or corresponding -i*x*) option specifies the "block file" xxxgnm.bk*x* to be constructed, where *X* is 'A', 'D' or 'P' and *x* is 'a', 'n' or 'p'. The -inp option will construct both xxxgnm.bkn (for cDNA queries)
and xxxgnm.bkp (for protein queries) files together with the xxxgnm.idx and associated files. 
-K*X* option is mandatory. If -i*x* is omitted or *x* is empty, xxxgnm.idx and associated files are created but no block file is constructed.  
 * The block size and *k*-mer size are estimated from the 
genome size unless explicitly specified (see below).  
 * If *MAX_GENE* (the length of the plausibly longest gene on the genome) is 
not specified, *MAX_GENE* is also estimated from the genome size.
  <u>Don't forget to specify *MAX_GENE* if xxxgnm.gf represents only a part of the genome!!</u>  Otherwise, *MAX_GENE* may be seriously underestimated.  
 * Options : (default value)
   * -g: The outputs except for X.grp are gzipped.
   * -t*N*: Number of threads. (1)
   * ~~-E: Generate local lookup table.~~
   * -yX:    Format for remote queries (more sensitive but less economic than default)
   * -XA*N*: Alphabet size of the reduced amino acids: 6 < *N* <= 20 (20)
   * -XB*S*: Bit patterns of the spaced seeds concatenated with commas. The pattern should be asymmetric  when the number of patterns > 2.
   * -XC*N*: Number of seed patterns: 0 <= *N* <= 5 (0: contiguous seed)
   * -XG*N*: Maximum gene length (inferred from genome size)
   * -Xa*N*: A parameter used to filter excessively abundant words (10)
   * -Xb*N*: Block size (inferred from genome size)
   * -Xk*N*: Word size (inferred from block size)
   * -Xs*N*: Distance between neighboring seeds (= *k*)

## <a name="Exec">Execution</a>

1. Prepare protein, cDNA, or genomic segment sequence(s) in (multi-)fasta format
(denoted by *query* below). From 2.3.2a, gzipped fasta file(s) may be used as 
the query without prior expansion. 
2. Store *query* to _work_.
3. `% cd work`
4. Run **spaln** in one of the following four modes. **Spaln**
    does not support comparison between two genomic segments.  
```
    (A) % spaln -Q[0|1|2|3] [-ON] [other options] genome_segment query
    (B) % spaln -Q[4|5|6|7] [-ON] [other options] -[d|D] xxxgnm query
    (C) % spaln -Q[4|5|6|7] [-ON] [other options] -[a|A] prosdb query
    (D) % spaln -Q[4|5|6|7] [-ON] [other options] prosdb.faa query
```
   * In the last case, *prosdb.faa* will be internally formatted, and the formatted results will be discarded after the end of execution.
   * Only a subset of queries may be examined if *query* is replaced with '*query* (from to)' (quotations are necessary), where 'from' and 'to' are the first and last entry numbers in *query* to be examined.  
   * Options: (default value)
     * -C *N*:	Use the genetic code specified by the "transl_table number" defined in [NCBI transl_table](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) (1).
     * ~~-E: Use local lookup table.~~
     * -H *N*: Output is suppressed if the alignment score is less than *N*. See also -pw. (35)
     * -K[D|P|A]: Format either genomic DNA for sequence search with DNA (D) or Protein (P) queries or protein database sequences for search with a genomic segment or protein queries (A). Use in combination with -W option below.
     * -LS: Smith-Waterman-type local alignment. This option may prune out weakly matched terminal regions.
     * -M[*N*[.*M*]]: Single or multiple output for each query
       * No option (default): single locus
       * No argument: Multiple loci up to the maximum number specified by the program (4 in the present implementation)
       * *N*=1: Re-search the *query* region not aligned in the first trial. May be useful to detect chimera or fragmented genomic region
       * *N*\>1: Output multiple loci maximally up to *N*
       * *M*: Maximal number of candidate loci to be subjected to spliced alignment (4). If *M* < *N*, *M* is reset to *N*.
     * -O *N*[,N<sub>2</sub>,N<sub>3</sub>...]: Select output format for genome vs cDNA or aa (4)  
      It is possible to output multiple files with extensions of .O*N* at a run by multiply applying this option. 
      Or by concatenating the format numbers with commas or colons, ex. -O0,1,4. See also -o option.

       * *N*=0: Gff3 gene format
       * *N*=1: Alignment
       * *N*=2: Gff3 match format
       * *N*=3: Bed format
       * *N*=4: Exon-oriented format similar to the output of megablast -D 3
       * *N*=5: Intron-oriented output
       * *N*=6: Concatenated exon sequence in extended (multi-)fasta format, in which the exon-intron structure of the parental gene is supplied by one
or more comment lines starting with ';C', such as  
       `;C complement(join(1232555..1232760,1233786..1233849,1233947..1234119,`  
       `;C 1234206..1234392))`
       * *N*=7: Translated amino-acid sequence. Presently not very useful for cDNA queries because the entire exon rather than an ORF is translated
       * *N*=8: Mapping (block) information only. Use with -Q4
       * *N*=10: SAM format
       * *N*=12: Output the same information as -O4 in the
	    binary format. If -oOutput is set, three files named
	    Output.grd, Output.erd, and Output.qrd will be created. Otherwise,
	    query.grd, query.erd, and query.qrd will be created. If the -g 
            option is set, gzip-compressed outputs will be generated.
       * *N*=15: Copy of the query sequence supplemented with the inferred gene structural information.  
  
     * -O *N*[,N<sub>2</sub>,N<sub>3</sub>...]: Select output format for aa vs aa (4)  
      It is possible to output multiple files with extensions of .O*N* at a run by multiply applying this option. 
      Or by concatenating the format numbers with commas or colons, ex. -O0,1,4. See also -o option.
       * *N*=0: statistics (%divergence alignment_score #match #mismatch #gap #unpaired)
       * *N*=1: Alignment
       * *N*=2: Sugar format
       * *N*=3: Psl format
       * *N*=4: XYL = Coordinate + match length
       * *N*=5: statistics + XYL
       * *N*=8: Cigar format
       * *N*=9: Vulgar format
       * *N*=10: SAM format
     * -Q *N*: Select algorithm (3)
        * 0<=*N*<=3: Genomic segment in the fasta format given by the first argument vs. *query* given by the following arguments. See also -i option below.  One may skip the formatting step described above if only this mode of operation is used.
        * 4<=*N*<=7: Genome mapping and alignment. The genomic sequence must be formatted beforehand.
        * *N*=0,4: DP procedure without HSP search. Considerably slow
        * *N*=1,2,3,5,6,7: Recursive HSP search up to the level of (*N* % 4)
     * -R *S*:	Read block index table from file *S*.
      If omitted, the xxxgnm.bkn, xxxgnm.bkp,
      or prosdb.bka file will be read depending on the type of query. The
      appropriate file is searched for in the current directory, the directory
      specified by the _env_ variable 'ALN_DBS', and the
      'seqdb' directory specified at the compile time in this order.
     * -S *N*:	Orientation of the DNA query sequence (3)
        * *N*=0: The
       orientation is inferred from the phrases (e.g. 5' end) in the header
       line of each entry within a fasta file. If no information is available,
       both orientations are examined, and the result with the better score is
       reported. Terminal polyA or polyT sequence is not trimmed.
        * *N*=1: Forward orientation only. PolyA tail may be trimmed off.
        * *N*=2: Reverse-complement orientation only. Leading polyT sequence may be trimmed off.
        * *N*=3: Examine both orientations. Terminal polyA or polyT sequence may be trimmed off.
     * -T *S*: Specify the species-specific parameter set. *S* corresponds to the subdirectory in the _table_ directory. Alternatively, *S* may be the 1st or the 3rd term in _table/gnm2tab_ file, where the 2nd term on the line indicates the subdirectory.
     * -V *N*:	Minimum space to induce Hirschberg's algorithm (16M)
     * -W *S*:	Write block index table to files *S*.bk*x*. if *S* is omitted, the file name (without directory and extension) of the first argument is used as *S*.
     * -g: gzipped output used in combination with -O12 option.
     * -i[a|p]: Input mode with -Q[0<=N<=3].
        * -ia: Alternative mode; a genomic segment of an odd numbered entry in the input file is aligned with the query of the following entry.
        * -ip: Parallel mode; the i-th entry in the file specified by the first argument is aligned with the i-th entry in the file specified by the second argument.
        *  	default: The genomic segment specified by the first argument is aligned with each entry in the file specified by the second argument.
     * -o *S*:	Destination of output file name (stdout). If multiple output formats are specified by -O option(s), *S* specifies the directory or prefix to which the file names with .O*N* extensions are concatenated.
      
     * -pa*N*:	Terminal polyA or polyT sequence longer than *N* (12) is trimmed off and the orientation is fixed accordingly. If *N* = 0 or empty, these functionalities are disabled.
     * -pi:	Mark exon-intron junctions by color in the alignment (-O1).
     * -pn:	Prohibit ovewrite of existing output files (ask).
     * -po:	Allow ovewrite of existing output files (ask).
     * -pq:	Suppress warning messages sent to *stderr*.
     * -pw:	Report result even if alignment score is below threshold value.
     * -px:	Suppress self-comparisons in the execution mode (C) or (D).
     * -pT:	Exclude termination codon from output (include).
     * -u *N*:	Gap-extension penalty (3, 2, 2)
     * -v *N*:	Gap-opening penalty (8, 6, 9)
     * -xB *S*: 	Bit pattern of seeds used for HSP search at level 1
     * -xb *S*:	Bit pattern of seeds used for HSP search at level 3
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
      to a conserved intron position (10)
     * -yL *N*:	Minimum intron length (30, 30)
     * -yS *N*:	*N* specifies the percentile
      contribution of the species-specific splice signal. The other part is
      derived from the universal signal given to the dinucleotides at the ends
      of an intron. An omission of *N* implies *N* = 100.
     * -yX *N*:	*N* = 0: set parameter values for intra-species comparison.
      *N* = 1: set parameter values for cross-species comparison. The default value for *N* is 0 or 1 for DNA or protein query, respectively.
     * -yY *N*:	Relative contribution of length-dependent part of intron penalty (8)
     * -yZ *N*:	Relative contribution of oligomer composition within an intron (0)
     * -XS: Activete salvage mode. Considerably slow.

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
## Changes in version 2.4.9
1. The -g option at the format time compresses all output files except for X.grp. According to this change, the format of X.grp is slightly modified, which does not necessitate reformat of existing files.
2. To guard existing files from accidental loss, the propriety of overwrite of each existing file is asked by default. The -pn (always skip) or -po (always overwrite) option evades this inquiry.
3. The default memory size used for core sort in *sortgrcd* has been changed from 1MB to 16MB.

## Changes in version 2.4.8
1. By default, termination codon is included in the output of any form with a protein query. Use -pT option to get the previous forms (without termination codon) of outputs.
2. Enable to read long (>= 2**31B) gzipped files by *spaln* and *sortgrcd*.
3. Fix a bug at the format time with -g option.

## Changes in version 2.4.7
1. The lower limit of query length has been halved, which improved mappability of short queries.
2. A serious bug for DNA queries with -LS option has been fixed.
3. A bug at the format mode (-W option) with gzipped input has been fixed.
4. Several other bugs in **spaln** and **sortgrcd** that cause occasional segmentation faults and memory leaks have been fixed.

## Changes in version 2.4.6
1. Fix several bugs that caused occasional segmentation faults and memory leaks.
2. The -E option has been abolished from both format and search modes.
3. Avoid DP calculation when realistic alignment is desperately obtainable, particularly in the local alignment mode with -LS option.
4. "Bad molecular type ..." error message at the run time is suppressed when the query molecular type is specified by -K option.

## Changes in version 2.4.5
1. Fix a serious bug for DNA query with -S3 option.
2. Back to the original default settings in the formatting mode. Use -yX option to format genomic sequence for more sensitive yet expensive block search than the standard setting.

## Changes in version 2.4.4
1. Fix a bug in preventing orphan exons in DP for DNA query.
2. Revival of the salvage mode with -XS option (examine all positively scored blocks).

## Changes in version 2.4.3
1. Small improvement in memory usage.
2. A few bugs related to the unidirectional Hirschberg algorithm have been fixed.
3. Illegal memory access in generating splice-junction sequence has been fixed.

## Changes in version 2.4.2
1. Minor improvements in the DP-based spliced alignment engines.
2. A few bugs related to the unidirectional Hirschberg algorithm have been fixed.
3. A bug in the Boyer-Moore algorithm has been fixed. This is relevant to protein queries with -yX0 option, and when the search goes beyond a block boundary.

## Changes in version 2.4.1
1. The algorithm for delimiting a genic region has been modified to find remote terminal coding exon(s) separated by long (up to 99.6% quantile) intron(s) from the main body of the gene.
2. The -yx0 option now tries to search for missing internal micro exons and terminal very short coding exons.
3. Selenocysteine (denoted by U) is now regarded as the 21th amino acid which favorably matches an in-frame TGA termination codon (U in the Tron code) upon DNA vs amino acid sequence alignment.
4. Gene candidates are now sorted according to the final alignment score rather than the intermediate chained HSP score. This modification has improved the chance of true orthologous hits rather than paralog hits at an expense of a slight increase in computational load.
5. Compared with the previous versions, a larger number of species-specific parameter sets (247 <- 102) are provided to support more species (1495 <- 688). **Note** that some parameter-set identifiers are changed. Please use eight-digit species identifies (e.g. zea_mays) rather than former parameter-set identifiers (e.g. Magnolio) as the argument of -T option.

## Changes in version 2.4.0
1. **Spaln** can now directly format genomic sequences without relying on 'make' command. See [Format](#Format).
2.  The internal format of index files is slightly modified. Although previously-formatted files can be used by the new version, the opposite is not true. Note that use of older files with the new version can lead to a slight loss in sensitivity.
3. The above change has been done to facilitate multi-thread operation at the format time, although the acceleration rate by multi-threading is only marginal.
4. Multiple output forms can be produced at a single run. See -O and -o options.
5. The traditional bidirectional Hirschberg algorithm is changed to the unidirectional variant.
6. Also, the bidirectional 'sandwich' or 'attack by both sides' spliced alignment algorithm has been changed to unidirectional 'skipped' spliced alignment algorithm. This and the above changes have considerably reduced code complexity.
7. Local lookup table (xxxgnm.lun or xxxgnm.lup) is generated and used with -E option. Be cautious to use this option, as a large disk space is required to store the generated file, and a large memory is required at the runtime.
8. Paired-ends mode has been removed.
9. Many small bugs have been fixed.

## Changes in version 2.3.3
1. The heuristic alignment engine has been updated, resulting in marginal but significant improvement in speed and accuracy, 
especially with -Q3/7 option.
2. chachr.pl has been extended to accept Ganbank/DDBJ and EMBL-formatted files in 
addition to FASTA files. Maybe used as a format convertor.
3. Update help and error messages of **Spaln** and **Sortgrcd**.
4. Prevent segmentation fault invoked with -ia or -ip option.
5. The maximal path size has been extended from 255 to 2047 characters.
6. The 'NEVSEL' constant value has been changed to avoid underflow of 2 * NEVSEL.
7. In utilseq.c and utilseq.h, a member variable in class PatMat has been moved to a local variable
to recover thread safety.
8. When **Spaln** is run: `% spaln protein genome`, the order of the 1st and 2nd arguments is exchanged with a warning message.

## Changes in version 2.3.2
1. From this version, query fasta file(s) may be compressed.
2. The new option of <b>spaln</b> '-g' directly generates compressed output(s) when used in combination with -W or -O12 option.
3. makeidx.pl and makblk.pl have been modified to accord with gzipped genome/database fasta files.
4. A small bug in makdbs.c has been fixed.
5. From this version, input genome/database fasta files (X.mfa, X.gf, or X.faa), formatted data files (X.seq,X.bka, X.bkn, and X.bkp) for **spaln**, and X.grd, X.erd, and X.qrd for **sortgrcd** may be gzipped if USE_ZLIB mode is activated upon [compilation](#compile). **Note:** other data files (X.ent, X.grp, X.idx, and X.odr) must not be compressed.
6. A serious bug concerning with multiple queries has been fixed. This has considerably improved mapping sensitivity particularly when -M option is set under single thread operation mode.
7. Fixation of several small bugs and fine tuning of codes further enhanced mapping sensitivity and specificity particularly for short protein queries.
8. -O *N* option of **sortgrcd** has been extended. -O0: Gff3: -O3: BED; -O4: exon-oriented; -O5: intron-oriented; -O6: concatenated exons; -O7: translated amino acid sequence; -O15: unique introns.

## <a name="Ref">References</a>

<a name="Ref1">[[1]](http://nar.oxfordjournals.org/cgi/content/abstract/gkn105?ijkey=N2yLVza41RuShAg&keytype=ref) Gotoh, O.
A space-efficient and accurate method for mapping and aligning cDNA sequences onto genomic sequence. *Nucleic Acids Research* **36** (8) 2630-2638 (2008).

<a name="Ref2">[[2]](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btn460?ijkey=XajuzvyHlcQZoQd&keytype=ref) Gotoh, O.
Direct mapping and alignment of protein sequences onto genomic sequence. *Bioinformatics* **24** (21) 2438-2444 (2008).

<a name="Ref3">[[3]](http://nar.oxfordjournals.org/content/40/20/e161) Iwata, H. and Gotoh, O.
Benchmarking spliced alignment programs including  Spaln2, an extended version of Spaln that incorporates additional species-specific features. *Nucleic Acids Research* **40** (20) e161 (2012)

<a name="Ref4">[[4]](https://doi.org/10.1093/bioinformatics/16.3.190) Gotoh, O.
Homology-based gene structure prediction: simplified matching algorithm using a translated codon (tron) and improved accuracy by allowing for long gaps. *Bioinformatics* **16** (3) 190-202 (2000)

<a name="Ref5">[[5]](https://academic.oup.com/bioinformatics/article/22/10/1211/236993) Nagasaki, H., Arita, M., Nishizawa, T., Suwa, M., Gotoh, O.
Automated classification of alternative splicing and transcriptional initiation and construction of a visual database of the classified patterns. *Bioinformatics* **22** (10) 1211-1216 (2006).

<a name="Ref6">[[6]](https://doi.org/10.1007/978-1-0716-1036-7_5) Gotoh, O.
Cooperation of Spaln and Prrn5 for construction of gene-structure-aware multiple sequence alignment. In: Katoh K. (eds) Multiple Sequence Alignment. *Methods in Molecular Biology* **2231**, Humana, New York, NY. (2021).

* * *

Copyright (c) 1997-2022 Osamu Gotoh (o.gotoh@aist.go.jp) All Rights Reserved.
