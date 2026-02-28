# Species-specific parameters used by aln/spaln

### Generate species-specific parameter sets
#### Present Version 1.0.2
#### Last updated: 2026-02-28

- [Summary](#Sm)
- [Preparation](#Pp)
- [Execution](#Exec)
- [Changes](#Chng)


* * *

## <a name="Sm">Summary</a>

Generate species-specific parameter files used by **Aln** and **Spaln** when 
the genomic DNA sequence (G_g.fna) and a bunch of transcripts (cDNA, CDS, 
EST, TSA or long NGS read) nucleotide sequences (C_c.cf) are given both 
in FASTA format. Short reads in FASTQ format are not accepted. Two directories, 
`Table` and `SeqDb`, play special roles as explained in README.md. Note that 
it is possible to do everything in a single working directory `Work`, i.e., 
`Work` = `Table` = `SeqDb`. 

## <a name="Pp">Preparation</a>

1. `make all` instead of `make` must be conducted at the [installation](README.md). 
Note that GNU gsl library is needed to make fitild.
2. Make sure that `BIN` and `PERL`, where executables and perl scripts are 
installed, are properly set within the PATH.
3. To visualize intron-length distribution with 'fitild', 'gnuplot' and X Window
   system must also be installed.
4. `% cd SeqDb` and store G_g.fna(.gz) there.
5. Format G_g.fna for DNA queries. Typicall,  
   `% spaln -W -KD -g -t4 G_g.fna(.gz)`.  
   Follow the instruction in [README.md](README.md), when G_g.fna represents 
only a part of the whole genome.

## <a name="Exec">Execution</a>

1. `% cd Table`.
2. `% mkdir yourspec`, where `yourspec` is the identifier of your species.
3. `% cd yourspec` and store C_c.cf(.gz) there. 
4. `% gunzip C_c.cf.gz` if your transcript sequence is compressed. Otherwise,
   the following calculation can be much prolonged.
5. `% make_eij.pl -d G_g [additional spaln options] C_c.cf`.
   The default spaln options are "-Q7 -KD -O12 -LC -pq".  
   It is a good idea to use -T option to specify an existing parameter set 
   taxonomically closest to your species. Otherwise, the predefined 
   'generic' parameter set will be used. Other common additional options are 
   "-t*N*" and "-A*N*".  
   Confirm that C_c.eij file has been generated. If the file consists of too 
   small a number of lines (number of unique introns), something went wrong. 
   Check the above-mentioned procedures carefully. It is desired that C_c.eij 
   consists of more than 5000 lines.  
6. `% make_ssp.pl -d G_g -S -g C_c.eij`  
   This command will generate Splice[3|5].dat and IntronPotTab.dat 
together with a few associated files. The optional -g option generates gzipped 
files (Splice[3|5].dgz and IntronPotTab.dgz) instead of the binary '.dat' files. 
'G.cano' shows the numbers of 
canonical and non-canonical splice junction tetramers. An excessively large 
percentage of non-canonical tetramers warns that something unusual has happened.  
7. If your species uses non-standard genetic code, or the nucleotide 
composition of your genome is highly biased, species-specific CodePotTab.[dat|dgz] 
might have to be obtained. For this purpose, run  
   `% make_ssp.pl -d G_g -S -g [-CN] -c CDS.fna`,  
   where 'CDS.fna' is a FASTA file containing CDS sequences, and *N* is
   the "transl_table number" defined in 
   [NCBI transl_table](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) (0 by default).
8. Conventions:
   * To visualize the observed and fitted intron-length distributions, run  
     `% fitild -g -a -d G.ildp G.ild`  
   * To convert Splice[3|5].[dat|dbz] to readable text form, run  
     `% npssm -f Splice[3|5].[dat|dbz] > Splice[3|5]`  or  
     `% npssm -f Splice[3|5].[dat|dbz] -o Splice[3|5]`  
   * Likewise, to convert CodePotTab.[dat|dbz] or IntronPotTab.[dat|dbz] to text form,
     run,  
     `% exinpot -f [CodePotTab|IntronPotTab].[dat|dbz] > [CodePotTab|IntronPotTab]` or  
     `% exinpot -f [CodePotTab|IntronPotTab].[dat|dbz] -o [CodePotTab|IntronPotTab]`  
9. **Spaln** Ver.3.0.8 accepts all forms of parameter files. However, **aln** and
earlier version of **Spaln** may fail to handle the binary or gzipped binary files. 
In that case, please take advantage of the above-mentioned conversion procedure.  
10. Although usually unnecessary, you may repeat 4-7 once again using the newly 
generated parameter files. This is feasible as the parameter files in the 
current directory are preferentially used.
11. Of the files generated, **aln/spaln** use AlnParam, CodePotTab[|.dat|.dgz], 
   IntronPotTab[|.dat|.dgz], and Splice[3|5][|.dat|.dgz]. 
   Not all of them must exist; 'generic' parameter set will cover the missing part.
   If desired, you may delete all or any other intermediate files.
12. `% cd Work`. Now you can use the newly generated parameter sets by adding 
`-T yourspec` as an option of **Aln** or **Spaln**.

## <a name="Chng">Changes from previous version</a>
## Changes in version 1.0.2
1. In addition to binary .dat form, gzipped binary .dgz form is usable. 

## Changes in version 1.0.1
1. In addition to text form, more compact binary form (.dat) of parameter files
are now generated and accessible.

## <a name="Ref">References</a>

<a name="Ref1">[[1]](https://doi.org/10.1093/bioinformatics/bty353) Gotoh, O. 
Modeling one thousand intron length distributions with Fitild, *Bioinformatics*
**34** (19) 3258-3264 (2018) 

<a name="Ref1">[[2]](https://doi.org/10.1093/bioinformatics/btae517) Gotoh, O. 
Spaln3: improvement in speed and accuracy of genome mapping and spliced 
alignment of protein query sequences, *Bioinformatics* **40** (8) btae517 (2024)

\* * *

Copyright (c) 2023-2026 Osamu Gotoh (gotoh.osamu.67a@st.kyoto-u.ac.jp) All Rights Reserved.
