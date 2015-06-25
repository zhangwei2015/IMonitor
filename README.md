#IMonitor
Introduction

IMonitor - analyze the sequence data of immune repertoire sequenced by NGS. If paired-end reads(FASTQ) as input, it will be merged to single sequence according to overlap region. FASTA sequence as input is acceptable. It provides re-alignment algorithm to identify accurately V,D,J alleles.Sequencing error will be corrected. CDR3 region, deletion/insertion of VDJ will be identied and translate nucleotide into amino acid. Multiple statistics and graphs will be provied.

System Requirement

It runs on 64-bit Linux systems. for 1Gb raw data as input, about maximum 2Gb memory would be required.
Perl and R need to be installed for you system

Installation

   1. Before use it, perl(https://www.perl.org/get.html)
   2. Before use it,R(http://www.r-project.org/) need to be installed. and provide the installation paths for parameter -Rs
   3. Download the IMonitor.tar.gz to your directory, uncompress it.

Usage

   1. Create shell
   perl IMonitor.pl
        Compulsory: for FASTQ format(paired-end read), -a -b -A1 -A2 -o -n -t -Rs; for FASTA format(single-end read), -i -o -n -t -Rs
        Optionally: others
        all the parameters have the detail introduction if you run "perl IMonitor.pl"
   this step will create multiple directory and shells in Bin/.

