Manual of IMonitor
=====
# Introduction

IMonitor - analyze the sequence data of immune repertoire sequenced by NGS. If paired-end reads(FASTQ) as input, it will be merged to single sequence according to overlap region. FASTA sequence as input is acceptable. It provides re-alignment algorithm to identify accurately V,D,J alleles.Sequencing error will be corrected. CDR3 region, deletion/insertion of VDJ will be identied and translate nucleotide into amino acid. Multiple statistics and graphs will be provied.

# System Requirement

It runs on 64-bit Linux systems. for 1Gb raw data as input, about maximum 2Gb memory would be required.
Perl and R need to be installed for you system

# Installation

   1. Before use it, perl(https://www.perl.org/get.html)
   2. Before use it,R(http://www.r-project.org/) need to be installed. and provide the installation paths for parameter -Rs
   3. Download the IMonitor.tar.gz to your directory, uncompress it.

# Usage

   1. Create shell
   perl IMonitor.pl
        Compulsory: for FASTQ format(paired-end read), -a -b -A1 -A2 -o -n -t -Rs; for FASTA format(single-end read), -i -o -n -t -Rs
        Optionally: others
        all the parameters have the detail introduction if you run "perl IMonitor.pl"
   this step will create multiple directory and shells in Bin/
        [parameters]
            \-a      <S> full path of input fq file 1
            \-b      <S> full path of input fq file 2
            -A1     <S> adaptor list file 1
            -A2     <S> adaptor list file 2
            -i      <S> single reads with FA format file
            -o      <S> output directory path
            -n      <S> sample name, used for prefix of ouput file
            -t      <S> gene type. e.g. TRB, IGH
            -k      <I> read length [100]

            -d      <F> add the paremeter means consider D genes for analysis. For IGH,TRB is necessary
            -jif    <I> J alignment identity for filtering,for C region sequenced,10 is recommended, or 80 is recommended[10]
            -vif    <I> V alignment identity for filtering [80]
            -r      <S> The reference directory [Bin/Ref/gene_type]
            -Q      <I> sequencing quality for filtering [15]
            -RQ     <Fl> qulity filter rate, used with -Q [0.05]
            -f1     <I> CDR3(nucleotide) abundance filter [0]
            -f2     <I> Clonotype(nucleotide, ful-length) abundance filter [0]
            -m          logical value, used to analyze hyper-mutation
            -ew         logical value, sequencing error correction for whole sequence. but this need a long time to run,only for FASTQ files as input
            -ec         logical value, sequencing error correction for only CDR3.only for FAST
            Q files as input
            -v      <I> used to calculate the base quality [64]

            -mul    <I> split the fa into multiple parts for alignment [3]
            -Rs     the R script directory.for a new system, this parameter is need to change[/opt/blc/genome/biosoft/R/bin/Rscript]
            -h      print help information

Note
            1. If  Pair-end(PE) sequencing FASTQ format as input, then:
                    perl IMonitor.pl
                    Compulsory: -a -b -A1 -A2 -o -n -t
                    Optionally: others
            2. If Single-end(SE) sequencing FASTA format as input, then:
                    perl IMonitor.pl
                    Compulsory: -i -o -n -t
                    Optionally: others, but -ew,-ec are invalid here

            Note: if sequence C region, -jif 10 is recommended, or -jif 80 is recommended
              The rate of IMonitor output is multipled 100%


   2. Run shell
   2.1 it can easy to run the general sh in Bin/ directory 'Execute_all.sh': sh Execute_all.sh
   2.2 run multiple shells in Bin/ directory seperately,so the *.blast.sh could be run in parallel.
        sh *.merge_fq_fq2fa.sh
        sh *.blast.*.sh
        sh *.structure.sh
        sh *.statistics.graph.sh
        sh *.rm.intermediate.file.sh

   2.3 some warning is normal when run blast alignment because of blast parameter setting. such as "[blastall] WARNING: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options".

# Output

   1. Directory
        --Bin
        --Align
        --Result
        --Figures

    use /Test data as a example, the output as follow:
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


  2. file format:
        2.1. Result/*.structure.gz format
        ID     fuction V_ref   D_ref   J_ref   CDR3_start      CDR3_end        CDR3(dna)       CDR3(aa)        3'V_del 5'D_del 3'D_del 5'J_del  VD_ins  DJ_ins  VJ_ins  strand  sequence        amino_acid      alignment_record
        2.2 Figures/*_overall_plot.pdf
        a general display for the statistics.

# Testing

      directory Test/ has a data for testing
      FQ_run.sh:
      perl ../IMonitor.pl -a data/XHS_1.fq.gz -b data/XHS_2.fq.gz -A1 data/1.adapter.list.gz -A2 data/2.adapter.list.gz -o . -n XHS -T TRB -k 100 -r /ifs1/ST_MED/USER/zhangwei/Immunity/Bioinf/pipeline/IMonitor_for_submit/Ref/TRB -d -m
      FA_run.sh:
      perl ../IMonitor.pl -i data/XHS.merged.fa.gz  -o . -n XHS -T TRB -k 100 -r /ifs1/ST_MED/USER/zhangwei/Immunity/Bioinf/pipeline/IMonitor_for_submit/Ref/TRB -d -m



