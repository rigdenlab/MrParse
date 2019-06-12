TWOUVO_SEQ = 'ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG'

PHMMER_LOG_TXT = """# phmmer :: search a protein sequence against a protein database
# HMMER 3.1b1 (May 2013); http://hmmer.org/
# Copyright (C) 2013 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             ../data/2uvoA.fasta
# target sequence database:        /opt/ccp4/ccp4-7.0/share/mrbump/data/pdb95.txt
# MSA of hits saved to file:       phmmerAlignment.log
# per-seq hits tabular output:     phmmerTblout.log
# per-dom hits tabular output:     phmmerDomTblout.log
# max ASCII text line length:      unlimited
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       2UVO:A|PDBID|CHAIN|SEQUENCE  [L=171]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------
   6.8e-105  351.3  93.4   7.6e-105  351.2  93.4    1.0  1  2x3t_C    
    2.1e-26   94.8  53.8    5.2e-26   93.5  53.8    1.5  1  1ulk_B    
    1.2e-18   69.4  33.2    1.3e-18   69.3  33.2    1.0  1  1ulm_B    
    7.1e-12   47.3  10.8    7.1e-12   47.3  10.8    1.0  1  4wp4_A    
    1.3e-11   46.4  37.0    1.5e-11   46.2  37.0    1.1  1  1eis_A    
    5.6e-11   44.4  36.1    6.3e-11   44.2  36.1    1.1  1  1iqb_B    
      7e-09   37.5  16.2    7.4e-09   37.4  16.2    1.0  1  4mpi_A    
    4.9e-08   34.7  13.5      5e-08   34.7  13.5    1.0  1  5wuz_A    
    1.1e-07   33.6  17.2    1.1e-07   33.6  17.2    1.0  1  5xdi_A    
    1.8e-06   29.7  20.9    1.8e-06   29.6  20.9    1.0  1  2lb7_A    
    6.5e-06   27.8   7.2    6.5e-06   27.8   7.2    1.0  1  1zuv_A    
    2.3e-05   26.0   8.1    2.3e-05   26.0   8.1    1.0  1  1mmc_A    
     0.0022   19.5  12.6      0.003   19.1  12.6    1.3  1  2kus_A    
     0.0039   18.7  11.7     0.0039   18.7  11.7    1.1  1  2n1s_A    
  ------ inclusion threshold ------
      0.023   16.2  18.7      0.023   16.2  18.7    1.1  1  1p9z_A    


Domain annotation for each sequence (and alignments):
>> 2x3t_C  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  351.2  93.4  2.4e-108  7.6e-105       2     170 ..       1     166 []       1     166 [] 1.00

  Alignments for each domain:
  == domain 1  score: 351.2 bits;  conditional E-value: 2.4e-108
  2UVO:A|PDBID|CHAIN|SEQUENCE   2 rcgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggcd 170
                                  rcgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgyc   cqsggcd
                       2x3t_C   1 RCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYC---CQSGGCD 166
                                  8**************************************************************************************************************************************************************...******9 PP

>> 1ulk_B  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   93.5  53.8   1.6e-29   5.2e-26      46     166 ..       4     120 ..       1     124 [. 0.54

  Alignments for each domain:
  == domain 1  score: 93.5 bits;  conditional E-value: 1.6e-29
  2UVO:A|PDBID|CHAIN|SEQUENCE  46 cgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqs 166
                                  cg +a+g  c +  ccsq+gycg   eycg gcq+  c  + +cg + ggk c + lccsq+g+cg +   cg gcqs  cs    cgkd ggr+ct + ccs++g cg+   +c  gcqs
                       1ulk_B   4 CGVRASGRVCPDGYCCSQWGYCGTTEEYCGKGCQS-QCDYN-RCGKEFGGKECHDELCCSQYGWCGNSDGHCGEGCQS-QCSYW-RCGKDFGGRLCTEDMCCSQYGWCGLTDDHCEDGCQS 120
                                  55556666666666666666666666666666654.34432.566666666666666666666666555566666655.35443.366666666666666666666666666666666655 PP

>> 1ulm_B  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   69.3  33.2   4.1e-22   1.3e-18      88     169 ..       3      81 ..       1      82 [] 0.88

  Alignments for each domain:
  == domain 1  score: 69.3 bits;  conditional E-value: 4.1e-22
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggc 169
                                  +cg +a+gk cpn  ccsqwg+cg    +cg gcqs  c     cg+d ggr+c  + ccsk+g cg +  +c  gcqs  c
                       1ulm_B   3 ECGERASGKRCPNGKCCSQWGYCGTTDNYCGQGCQSQ-CDY-WRCGRDFGGRLCEEDMCCSKYGWCGYSDDHCEDGCQSQ-C 81 
                                  6999999999999999999999999999999999985.765.46999999999999999999999999999999999984.5 PP

>> 4wp4_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   47.3  10.8   2.2e-15   7.1e-12      88     124 ..       2      40 ..       1      43 [] 0.87

  Alignments for each domain:
  == domain 1  score: 47.3 bits;  conditional E-value: 2.2e-15
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgg..gcqsg 124
                                  +cg qaggklcpnnlccsqwg+cg   e+c+   +cqs+
                       4wp4_A   2 QCGRQAGGKLCPNNLCCSQWGWCGSTDEYCSPdhNCQSN 40 
                                  6*****************************852268875 PP

>> 1eis_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   46.2  37.0   4.7e-15   1.5e-11       2      79 ..       1      82 [.       1      85 [] 0.76

  Alignments for each domain:
  == domain 1  score: 46.2 bits;  conditional E-value: 4.7e-15
  2UVO:A|PDBID|CHAIN|SEQUENCE  2 rcgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtsk....rcgsqaggatctnnqccsqygycgfgaeycgag.cq 79
                                 rcg qg    cp   ccs +g+cg +  ycg+ c+n  cw+ +    rcg+  g+  c  ++ccs +g+cg g +yc+ g cq
                       1eis_A  1 RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCEN-KCWSGErsdhRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGnCQ 82
                                 788888888888888888888888888888888887.5887542223688888888888888888888888888888554255 PP

>> 1iqb_B  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   44.2  36.1     2e-14   6.3e-11       2      79 ..       1      82 [.       1      88 [] 0.82

  Alignments for each domain:
  == domain 1  score: 44.2 bits;  conditional E-value: 2e-14
  2UVO:A|PDBID|CHAIN|SEQUENCE  2 rcgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtsk....rcgsqaggatctnnqccsqygycgfgaeyc.gagcq 79
                                 rcg qg    cp   ccs +g+cg +  ycg+ c+n  cw+ +    rcg+  g+  c  ++ccs +g+cg g +yc g+ cq
                       1iqb_B  1 RCGSQGGGGTCPALWCCSIWGWCGDSEPYCGRTCEN-KCWSGErsdhRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCsGSKCQ 82
                                 899999999999999999999999999999999998.6997652223799999999999999999999999999999445566 PP

>> 4mpi_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   37.4  16.2   2.3e-12   7.4e-09      88     124 ..       3      39 ..       1      43 [. 0.79

  Alignments for each domain:
  == domain 1  score: 37.4 bits;  conditional E-value: 2.3e-12
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg 124
                                  +cg qagg lcp  lccsq+g+c+   e+cg+gcqs 
                       4mpi_A   3 QCGRQAGGALCPGGLCCSQYGWCANTPEYCGSGCQSQ 39 
                                  5888888888888888888888888888888888874 PP

>> 5wuz_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   34.7  13.5   1.6e-11     5e-08      88     124 ..       2      40 ..       1      43 [] 0.70

  Alignments for each domain:
  == domain 1  score: 34.7 bits;  conditional E-value: 1.6e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcg..ggcqsg 124
                                   cg qag++ c n lccsq+gfcg  se+c+  +gcqs+
                       5wuz_A   2 NCGRQAGNRACANQLCCSQYGFCGSTSEYCSraNGCQSN 40 
                                  478888888888888888888888888888533577775 PP

>> 5xdi_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   33.6  17.2   3.5e-11   1.1e-07      45      79 ..       2      37 ..       1      40 [] 0.88

  Alignments for each domain:
  == domain 1  score: 33.6 bits;  conditional E-value: 3.5e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE 45 rcgsqaggatctnnqccsqygycgfgaeycgag.cq 79
                                 +cg qagga c+n  ccsq+gycg    ycgag cq
                       5xdi_A  2 QCGRQAGGARCSNGLCCSQFGYCGSTPPYCGAGqCQ 37
                                 69*****************************99555 PP

>> 2lb7_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   29.6  20.9   5.6e-10   1.8e-06      44      84 ..       2      42 ..       1      44 [] 0.64

  Alignments for each domain:
  == domain 1  score: 29.6 bits;  conditional E-value: 5.6e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE 44 krcgsqaggatctnnqccsqygycgfgaeycgagcqggpcr 84
                                 +rcg qa ga c n  cc +yg+cg g  ycgag   + cr
                       2lb7_A  2 QRCGDQARGAKCPNCLCCGKYGFCGSGDAYCGAGSCQSQCR 42
                                 56777777777777777777777777777777765445565 PP

>> 1zuv_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   27.8   7.2   2.1e-09   6.5e-06      97     118 ..       8      29 ..       1      30 [] 0.60

  Alignments for each domain:
  == domain 1  score: 27.8 bits;  conditional E-value: 2.1e-09
  2UVO:A|PDBID|CHAIN|SEQUENCE  97 lcpnnlccsqwgfcglgsefcg 118
                                   cp+ +ccsqwg+cg g ++cg
                       1zuv_A   8 RCPSGMCCSQWGYCGKGPKYCG 29 
                                  3666666666666666666665 PP

>> 1mmc_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   26.0   8.1   7.3e-09   2.3e-05      11      33 ..       8      30 .]       2      30 .] 0.85

  Alignments for each domain:
  == domain 1  score: 26.0 bits;  conditional E-value: 7.3e-09
  2UVO:A|PDBID|CHAIN|SEQUENCE 11 ecpnnlccsqygycgmggdycgk 33
                                  cp+ +ccsq+gycg g  ycg+
                       1mmc_A  8 RCPSGMCCSQFGYCGKGPKYCGR 30
                                 69999999999999999999996 PP

>> 2kus_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   19.1  12.6   9.4e-07     0.003      45      75 ..       6      34 ..       2      35 .] 0.70

  Alignments for each domain:
  == domain 1  score: 19.1 bits;  conditional E-value: 9.4e-07
  2UVO:A|PDBID|CHAIN|SEQUENCE 45 rcgsqaggatctnnqccsqygycgfgaeycg 75
                                 +cg   gg  c    ccsqygycg g +yc+
                       2kus_A  6 QCGPGWGG--CRGGLCCSQYGYCGSGPKYCA 34
                                 56654444..778888888888888888884 PP

>> 2n1s_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   18.7  11.7   1.2e-06    0.0039      11      32 ..       9      30 .]       1      30 [] 0.50

  Alignments for each domain:
  == domain 1  score: 18.7 bits;  conditional E-value: 1.2e-06
  2UVO:A|PDBID|CHAIN|SEQUENCE 11 ecpnnlccsqygycgmggdycg 32
                                  c   lccs+ygycg g  ycg
                       2n1s_A  9 RCSGGLCCSKYGYCGSGPAYCG 30
                                 3555555555555555555554 PP

>> 1p9z_A  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   16.2  18.7   7.1e-06     0.023      53      84 ..       8      39 ..       1      40 [] 0.74

  Alignments for each domain:
  == domain 1  score: 16.2 bits;  conditional E-value: 7.1e-06
  2UVO:A|PDBID|CHAIN|SEQUENCE 53 atctnnqccsqygycgfgaeycgagcqggpcr 84
                                   c    ccs ygycg ga ycgag     cr
                       1p9z_A  8 RPCNAGLCCSIYGYCGSGAAYCGAGNCRCQCR 39
                                 46888999999999999999999985555555 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (171 nodes)
Target sequences:                        47841  (11034331 residues searched)
Passed MSV filter:                      1970  (0.0411781); expected 956.8 (0.02)
Passed bias filter:                      641  (0.0133985); expected 956.8 (0.02)
Passed Vit filter:                        51  (0.00106603); expected 47.8 (0.001)
Passed Fwd filter:                        15  (0.000313539); expected 0.5 (1e-05)
Initial search space (Z):              47841  [actual number of targets]
Domain search space  (domZ):              15  [number of targets reported over threshold]
# CPU time: 0.33u 0.02s 00:00:00.35 Elapsed: 00:00:00.14
# Mc/sec: 13477.65
//
# Alignment of 14 hits satisfying inclusion thresholds saved to: phmmerAlignment.log
[ok]
"""
