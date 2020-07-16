Manual of IMonitor
=====
# Introduction

IMonitor - analyze the sequence data of immune repertoire sequenced by NGS. If paired-end reads(FASTQ) as input, it will be merged to single sequence according to overlap region. FASTA sequence as input is acceptable. It provides re-alignment algorithm to identify accurately V,D,J alleles.Sequencing error will be corrected.CDR3 region can be identified by both VJ gene assignment and conserative region. deletion/insertion of VDJ will be identied and translate nucleotide into amino acid. Multiple statistics and graphs will be provied.

# System Requirement

It runs on 64-bit Linux systems. for 1Gb raw data as input, about maximum 2Gb memory would be required.
Perl and R need to be installed for you system

# Installation

   1. Before use it, perl(https://www.perl.org/get.html)
   2. Before use it,R(http://www.r-project.org/) need to be installed. and provide the installation paths for parameter -Rs
   3. Download the IMonitor-1.4.1.tar.gz to your directory, uncompress it.
      tar -zxvf IMonitor-1.4.1.tar.gz

# Version 1.4.1
   

# Usage
   !!!Note:  for the usage of 1.4.1, please run "perl IMonitor.pl" to see the details!!!

###1. Create shell
   
   perl IMonitor.pl
        Compulsory: for FASTQ format(paired-end read), -a -b -o -n -t -Rs; for FASTA format(single-end read), -i -o -n -t -Rs
        Optionally: others
        all the parameters have the detail introduction if you run "perl IMonitor.pl"
   this step will create multiple directory and shells in Bin/
       
         [parameters]
            -a      <S> full path of input fq file 1
            -b      <S> full path of input fq file 2
            -i      <S> single reads with FASTA format file
            -iq     <S> single reads with FASTQ format file
            -o      <S> output directory path
            -n      <S> sample name, used for prefix of ouput file
            -t      <S> gene type. e.g. TRB, IGH 
            -k      <I> read length [100]

            -d      <F> add the paremeter means consider D genes for analysis. For IGH,TRB is necessary
            -c          logical value, for without V or J sequence, find the CDR3 by conservative region. V(YXC),J([WF]GXG)
            -jif    <I> J alignment identity for filtering [80]
            -vif    <I> V alignment identity for filtering [80]
            -r      <S> The reference directory [Bin/Ref/gene_type]
            -Q      <I> sequencing quality for filtering [15]
            -RQ     <Fl> qulity filter rate, used with -Q [0.05]
            -f1     <I> CDR3(nucleotide) abundance filter [0]
            -f2     <I> Clonotype(nucleotide, ful-length) abundance filter [0]
            -m          logical value, used to analyze hyper-mutation
            -ew         logical value, sequencing error correction for whole sequence. but this need a long time to run,only for FASTQ files as input
            -ec         logical value, sequencing error correction for only CDR3.only for FASTQ files as input
            -Qe     <I> the quality(as a cutoff) is used for sequencing error correction. [20](-ew or -ec is required) 
            -adseq  <S>     3' adapter sequence
            -adcut  <I>     cut sequence from adaptor index,unless performed -f/-r also in use
                                    discard the read when the adaptor index of the read is less than INT

            the next two options only for sequencing type
            -v      <I> used to calculate the base quality [64]
            -seqType <I> Sequence fq type, 0->old fastq(Hiseq2000), 1->new fastq(Hiseq4000)[default: 0]
            old fastq id: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2
            new fastq id: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC
        
            -mul    <I> split the fa into multiple parts for alignment [3]
            -Rs     the R script directory.for a new system, this parameter is need to change[/ifs1/ST_MED/USER/zhangwei/software/R-3.0.2/bin/Rscript]
            -h      print help information

   Note:
   #1. If  Pair-end(PE) sequencing FASTQ format as input, then:
           perl IMonitor.pl
        Compulsory: -a -b -o -n -t -Rs
       Optionally: others
   #2. If Single-end(SE) sequencing FASTA format as input, then:
      perl IMonitor.pl
      Compulsory: -i -o -n -t -Rs
      Optionally: others, but -ew,-ec are invalid here
   #3. If Single-end(SE) sequencing FASTQ format as input, then:
      perl IMonitor.pl
      Compulsory: -iq -o -n -t -Rs
      Optionally: others
   
   For Zebra sequencing(Single-end), Compulsory: -iq -o -n -t -Rs -v 33 [-Qe 25 (-ew or -ec is required)]


   Note:
   The rate of IMonitor output is multipled 100%


###2. Run shell
   
   1. it can easy to run the general sh in Bin/ directory 'Execute_all.sh'
    
            sh Execute_all.sh

   2. run multiple shells in Bin/ directory seperately,so the *.blast.sh could be run in parallel.
   
            sh *.merge_fq_fq2fa.sh
            sh *.blast.*.sh
            sh *.structure.sh
            sh *.statistics.graph.sh
            sh *.rm.intermediate.file.sh

   3. Warning:
      some warning is normal when run blast alignment because of blast parameter setting. such as "[blastall] WARNING: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options".`
   4. Reference(/Ref/*)
   In the directory /Ref/*, it provides Human(IGH,TRB,IGK/L) germline sequences as the reference for analysis. These germline were dowoloaded from IMGT(www.imgt.org/).
   However, you can provide your own reference, but it need to use a script for initialization. see example Ref/TRB/run.sh
            
            sh Ref/bin/run.sh <V_sort.fa> <J_sort.fa> <cdr3region> <primer.txt> <gene_type> <out_dir> 1 [<D_sort.fa>]
            If create the normal reference, then <flag> = 1 && <primer.txt> is invail and use any file input is OK
            <V_sort.fa>/<J_sort.fa>/<D_sort.fa>: V/D/J germline FASTA sequences, nucleotide sequences with gaps according to the IMGT unique numbering.
                    gaps criterion: the conservative region( V(YXC),J([WF]GXG)) must be at the same positioins among all sequences. We need to determine the start and end positions for CDR3 region(<cdr3region>).
                    <V_sort.fa>: you can download from "IMGT reference directory"(http://www.imgt.org/vquest/refseqh.html)
                    <J_sort.fa>: you can add the gaps by yourself according to the IMGT database(http://www.imgt.org/IMGTrepertoire/Proteins/index.php#B or http://www.imgt.org/IMGTrepertoire/Proteins/)
                    
            <cdr3region>: the position of CDR3 start in V and CDR3 end in J. e.g. V310J25, means CDR3 start from the 310the position of V, end at the 25th position of J
            <gene_type>: gene name, e.g. IGH,TRB,TRA,IGKL

# Output

###1. Directory
   
             --Bin
             --Align
             --Result
              --Figures

    the output details as follow:
    
            |-- Align
            |-- Bin
            |   |-- Execute_all.sh
            |   |-- *.blast.1.sh
            |   |-- *.blast.2.sh
            |   |-- *.blast.3.sh
            |   |-- *.merge_fq_fq2fa.sh
            |   |-- *.rm.intermediate.file.sh
            |   |-- *.statistics.graph.sh
            |   |-- *.structure.sh
            |   `-- error.log
            |-- Figures
            |   |-- *_CDR3_freq.pdf
            |   |-- *_CDR3_nt_length_dis.pdf
            |   |-- *_Insert_size_dis.pdf
            |   |-- *_J_NT_composition.pdf
            |   |-- *_J_usage.pdf
            |   |-- *_VJ_pairing_3d.pdf
            |   |-- *_V_NT_composition.pdf
            |   |-- *_V_usage.pdf
            |   |-- *_mutation_dis.pdf
            |   |-- *_overall_plot.pdf
            |   |-- *_rarefract_curve.pdf
            |   |-- *_vdj_del_len.pdf
            |   |-- *_vdj_ins_len.pdf
            |   `-- *_vdj_len_inCDR3.pdf
            |-- Result
            |   |-- Statistics.pl
            |   |-- *.change_id.backup.gz
            |   |-- *.discarded.gz
            |   |-- *.lowabundance.filter.gz
            |   |-- *.merged.fa.gz
            |   |-- *.structure.gz
            |   |-- *_CDR3_AA.frequency.gz
            |   |-- *_CDR3_AA_section.stat
            |   |-- *_CDR3_NT.frequency.gz
            |   |-- *_CDR3_NT.length
            |   |-- *_Clonotype_AA.frequency.gz
            |   |-- *_Clonotype_NT.frequency.gz
            |   |-- *_J_NT_composition.txt
            |   |-- *_J_gene.usage
            |   |-- *_VDJ_deletion_nt_len.dis
            |   |-- *_VDJ_insertion_nt_len.dis
            |   |-- *_VDJ_length_inCDR3.dis
            |   |-- *_VJ_pairing.usage
            |   |-- *_V_NT_composition.txt
            |   |-- *_V_gene.usage
            |   |-- *_bascial_filter_stat.txt
            |   |-- *_further.stat.txt
            |   |-- *_hypermutation.stat
            |   |-- *_insert_size_len
            |   |-- *_rarefraction_curve.txt
            |   `-- *_seq_error_discard.gz

###2. file format
1. Result/*.structure.gz format
            
            ID     fuction V_ref   D_ref   J_ref   CDR3_start      CDR3_end        CDR3(dna)       CDR3(aa)        3'V_del 5'D_del 3'D_del 5'J_del  VD_ins  DJ_ins  VJ_ins  strand  sequence        amino_acid      alignment_record

2. Figures/*_overall_plot.pdf: a general display for the statistics.

# Testing

directory Test/ has a data for testing

            FQ_run.sh:
            perl ../IMonitor.pl -a data/XHS_1.fq.gz -b data/XHS_2.fq.gz -o . -n XHS -T TRB -k 100 -r ../Ref/TRB -d -m -Rs /opt/blc/genome/biosoft/R/bin/Rscript
            FA_run.sh:
            perl ../IMonitor.pl -i data/XHS.merged.fa.gz  -o . -n XHS -T TRB -k 100 -r ../Ref/TRB -d -m -Rs /opt/blc/genome/biosoft/R/bin/Rscript
            SE_run.sh:
            perl ../IMonitor.pl -iq data/Zebra_test.fq.gz -o . -n Zebra -t TRB --seqType 1 -v 33 -d -ec -mul 1 -Rs /data/Public_tools/R-4.0.2/bin/Rscript
            
#Reference
Zhang W, Du Y, Su Z, Wang C, Zeng X, Zhang R, Hong X, Nie C, Wu J, Cao H, et al: IMonitor: A Robust Pipeline for TCR and BCR Repertoire Analysis. Genetics 2015, 201:459-472.


