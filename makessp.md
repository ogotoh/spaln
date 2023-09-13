# Species-specific parameters used by aln/spaln

### Generate species-specific parameter sets
#### Present Version 1.0.0
#### Last updated: 2023-09-12

- [Summary](#Sm)
- [Preparation](#Pp)
- [Execution](#Exec)

* * *

## <a name="Sm">Summary</a>

Generate species-specific parameter files used by **Aln** and **Spaln** when 
the genomic DNA sequence (G_g.fna) and a bunch of transcript (cDNA, CDS, 
EST, or TSA) nucleotide sequences (C_c.cf) are given both in FASTA format. 
Short reads in FASTQ format are not accepted. Two directories, `Table` and 
`SeqDb`, play special roles as explained in README.md. Note that it 
is possible to do everything in a single working directory `Work`, i.e., 
`Work` = `Table` = `SeqDb`.

## <a name="Pp">Preparation</a>

1. `make all` instead of `make` must be conducted at the [installation](README.md). 
Note that GNU gsl library is needed to make fitild .
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
4. `% make_eij.pl -d G_g [additional spaln options] C_c.cf(.gz)`.  
   The default spaln options are "-Q7 -O12 -yX0 -LS -pq".  
   It is a good idea to use -T option to specify an existing parameter set 
   taxonomically closest to your species. Otherwise, the predefined 
   'generic' parameter set will be used. Other common additional options are 
   "-t _N_" and "-A _N_".  
   Confirm that C_c.eij file has been generated. If the file consists of too 
   small a number of lines (number of unique introns), something went wrong. 
   Check the above-mentioned procedures carefully. It is desired that C_c.eij 
   consists of more than 5000 lines.  
5. `% make_ssp.pl -d G_g -S -e9,16 C_c.eij`  
   This command will generate 'Splice5.dat', 'Splice3.dat', and 'AlnParam', 
together with a few associated files. 'G.cano' shows the numbers of canonical 
and non-canonical splice junction tetramers. An excessively large percentage of 
non-canonical tetramers warns that something unusual has happened.  
   To visualize the observed and fitted intron-length distributions, run  
   `% fitild -g -a -d G.ildp G.ild`
6. To obtain intron potential, '13,' should be added to the -e option,
i.e. -e9,13,16.
7. If your species uses non-standard genetic code, or the nucleotide 
composition of your genome is highly biased, species-specific CodePotTab.dat 
might have to be obtained. For this purpose, run  
   `% make_ssp.pl -d G_g -S [-C transl_table_number] -c CDS.fna`,  
   where 'CDS.fna' is a FASTA file containing CDS sequences.
8. Although usually unnecessary, you may repeat 4-7 once again using the newly 
generated parameter files. This is feasible as the parameter files in the 
current directory are preferentially used.
9. Of the files generated, **aln/spaln** use AlnParam, CodePotTab.dat, 
   IntronPotTab.dat, Splice3.dat, and Splice5.dat. Not all of them must 
   exist; 'generic' parameter set will cover the lacking part. You may 
   delete here other intermediate files.
10. `cd Work`. Now you can use the newly generated parameter sets by adding 
`-T yourspec` as an option of **Aln** or **Spaln**.

## <a name="Ref">References</a>

<a name="Ref1">[[1]](https://doi.org/10.1093/bioinformatics/bty353) Gotoh, O. 
Modeling one thousand intron length distributions with Fitild, *Bioinformatics*
**34** (19) 3258-3264 (2018) 

\* * *

Copyright (c) 2023 Osamu Gotoh (gotoh.osamu.67a@st.kyoto-u.ac.jp) All Rights Reserved.
