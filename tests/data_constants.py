TWOUVO_SEQ = 'ERCGEQGSNMECPNNLCCSQYGYCGMGGDYCGKGCQNGACWTSKRCGSQAGGATCTNNQCCSQYGYCGFGAEYCGAGCQGGPCRADIKCGSQAGGKLCPNNLCCSQWGFCGLGSEFCGGGCQSGACSTDKPCGKDAGGRVCTNNYCCSKWGSCGIGPGYCGAGCQSGGCDG'
FIVEHXG_SEQ = 'MRYFFMAEPIRAMEGDLLGVEIITHFASSPARPLHPEFVISSWDNSQKRRFLLDLLRTIAAKHGWFLRHGLFCIVNIDRGMAQLVLQDKDIRALLHAMLFVELQVAEHFSCQDNVLVDPLIHALHKQPNPLWLGDLGVGNATAAPLVCGCFSGVKLDRSFFVSQIEKMTFPLLVKHIRHYCDKIVVGGQENARYLPALKTAGIWATQGTLFPSVALEEIETLLLGSRMNTLRESNMGTMHTSELLKHIYDINLSYLLLAQRLIVQDKASAMFRLGINEEMANTLGALSLPQMVKLAETNQLVCHFRFDDHQTITRLTQDSRVDDLQQIHTGIMLSTRLLNEVDDTARKKRA'

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

PHMMER_AF_LOG_TXT = """# phmmer :: search a protein sequence against a protein database
# HMMER 3.2 (June 2018); http://hmmer.org/
# Copyright (C) 2018 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             ../data/2uvoA.fasta
# target sequence database:        /Users/adamsimpkin/opt/clean/ccp4-7.1/share/mrparse/data/af2_sequences.fasta
# MSA of hits saved to file:       phmmerAlignment_af2.log
# per-seq hits tabular output:     phmmerTblout_af2.log
# per-dom hits tabular output:     phmmerDomTblout_af2.log
# max ASCII text line length:      unlimited
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       2UVO:A|PDBID|CHAIN|SEQUENCE  [L=171]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence              Description
    ------- ------ -----    ------- ------ -----   ---- --  --------              -----------
    5.4e-74  253.2  98.2    7.1e-74  252.8  98.2    1.1  1  AFDB:AF-Q0JF21-F1      Lectin
    2.8e-11   48.2  10.0    2.8e-11   48.2  10.0    2.2  2  AFDB:AF-Q7Y1Z1-F1      Chitinase 7
    3.9e-10   44.5  16.7    3.9e-10   44.5  16.7    2.8  1  AFDB:AF-A0A1D6LMS5-F1  Chitinase
    1.1e-09   43.0  15.2    1.1e-09   43.0  15.2    2.4  2  AFDB:AF-Q7DNA1-F1      Chitinase 2
    3.4e-09   41.4  10.5    3.4e-09   41.4  10.5    3.0  3  AFDB:AF-Q9SDY6-F1      Chitinase
    1.1e-08   39.8   8.0    1.1e-08   39.8   8.0    1.7  2  AFDB:AF-C6T7J9-F1      Chitinase
      2e-08   38.9  10.9      2e-08   38.9  10.9    2.7  2  AFDB:AF-I1MMY2-F1      Chitinase
      3e-08   38.3   8.5      3e-08   38.3   8.5    2.3  2  AFDB:AF-Q6K8R2-F1      Chitinase 6
    1.3e-07   36.3  15.4    1.3e-07   36.3  15.4    1.9  1  AFDB:AF-B6TR38-F1      Chitinase
    1.5e-07   36.0  10.1    1.5e-07   36.0  10.1    2.2  2  AFDB:AF-B6TT00-F1      Chitinase
    1.7e-07   35.9  22.6    1.7e-07   35.9  22.6    2.7  2  AFDB:AF-Q42993-F1      Chitinase 1
      2e-07   35.7  21.9      2e-07   35.7  21.9    2.8  1  AFDB:AF-P25765-F1      Chitinase 12
    2.7e-07   35.2   9.6    2.7e-07   35.2   9.6    2.6  3  AFDB:AF-O24603-F1      Endochitinase CHI
      3e-07   35.1  14.9      3e-07   35.1  14.9    2.0  2  AFDB:AF-P43082-F1      Hevein-like preproprotein
    4.5e-07   34.5  19.1    4.5e-07   34.5  19.1    2.5  1  AFDB:AF-A0A1D6LMS1-F1  Chitinase
    4.8e-07   34.4  17.6    4.8e-07   34.4  17.6    2.4  2  AFDB:AF-A0A1D6GWN3-F1  Chitinase
      5e-07   34.3  12.1      5e-07   34.3  12.1    2.0  2  AFDB:AF-O24598-F1      Endochitinase At2g43580
    5.3e-07   34.2   8.6    5.3e-07   34.2   8.6    1.9  2  AFDB:AF-Q9M2U5-F1      Endochitinase EP3
    6.8e-07   33.9  11.4    6.8e-07   33.9  11.4    1.9  2  AFDB:AF-I1M587-F1      Chitinase
      2e-06   32.4  10.9      2e-06   32.4  10.9    2.0  2  AFDB:AF-O22841-F1      Endochitinase At2g43620
    2.3e-06   32.2  13.0    2.3e-06   32.2  13.0    1.9  1  AFDB:AF-O22842-F1      Endochitinase At2g43610
    3.6e-06   31.5  11.0    3.6e-06   31.5  11.0    2.4  2  AFDB:AF-A0A1D6GWN1-F1  Chitinase
      5e-06   31.1  27.6      5e-06   31.1  27.6    2.7  1  AFDB:AF-P24626-F1      Chitinase 3
    6.2e-06   30.8  14.1    6.2e-06   30.8  14.1    2.9  3  AFDB:AF-O24658-F1      Endochitinase At2g43590
    7.7e-06   30.4  17.6    9.7e-06   30.1  17.6    1.2  1  AFDB:AF-Q0JC38-F1      Os04g0493600 protein
      8e-06   30.4  12.2      8e-06   30.4  12.2    2.6  2  AFDB:AF-P19171-F1      Basic endochitinase B
    1.1e-05   29.9  21.2    2.4e-05   28.8  21.2    1.5  1  AFDB:AF-A0A1X7YIJ7-F1  Chitinase
    1.6e-05   29.4  15.2    1.6e-05   29.4  15.2    2.3  2  AFDB:AF-O04138-F1      Chitinase 4
    5.8e-05   27.6  25.6    5.8e-05   27.6  25.6    2.1  2  AFDB:AF-P29023-F1      Endochitinase B
    0.00012   26.5  24.9    0.00012   26.5  24.9    2.2  3  AFDB:AF-C0P451-F1      Chitinase
    0.00036   25.0  27.8    0.00036   25.0  27.8    1.9  2  AFDB:AF-I1NCA0-F1      Chitinase
     0.0013   23.2   0.3     0.0013   23.2   0.3    1.0  1  AFDB:AF-A0A0P0Y930-F1  Os12g0238550 protein
     0.0014   23.0   7.4     0.0014   23.0   7.4    2.1  2  AFDB:AF-O24654-F1      Inactive endochitinase At2g43600
     0.0015   23.0   6.8     0.0015   23.0   6.8    2.2  2  AFDB:AF-A0A1P8AME8-F1  Chitinase family protein
     0.0022   22.4  33.5     0.0022   22.4  33.5    2.5  2  AFDB:AF-Q688M5-F1      Chitinase 9
     0.0023   22.4  25.4     0.0023   22.4  25.4    1.9  2  AFDB:AF-I1NCA1-F1      Chitinase


Domain annotation for each sequence (and alignments):
>> AFDB:AF-Q0JF21-F1  Lectin
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  252.8  98.2     7e-78   7.1e-74       2     169 ..      30     197 ..      29     201 .. 0.99

  Alignments for each domain:
  == domain 1  score: 252.8 bits;  conditional E-value: 7e-78
  2UVO:A|PDBID|CHAIN|SEQUENCE   2 rcgeqgsnmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggc 169
                                   cg+q+  m cp+nlccsq+gycg+g dycg gcq+gac +s+rcgsq ggatc+nnqccsqygycgfg+eycg+gcq+gpcradikcg  a+g+lcpnn+ccsqwg+cglgsefcg+gcqsgac  +k cgk agg  c nn+ccs  g cg+g +ycg+gcqsggc
            AFDB:AF-Q0JF21-F1  30 TCGKQNDGMICPHNLCCSQFGYCGLGRDYCGTGCQSGACCSSQRCGSQGGGATCSNNQCCSQYGYCGFGSEYCGSGCQNGPCRADIKCGRNANGELCPNNMCCSQWGYCGLGSEFCGNGCQSGACCPEKRCGKQAGGDKCPNNFCCSAGGYCGLGGNYCGSGCQSGGC 197
                                  6*********************************************************************************************************************************************************************** PP

>> AFDB:AF-Q7Y1Z1-F1  Chitinase 7
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   48.2  10.0   2.8e-15   2.8e-11      83     136 ..      30      81 ..      20      97 .. 0.86
   2 ?   -2.2   1.0       7.6   7.7e+04      34      53 ..     179     198 ..     159     230 .. 0.48

  Alignments for each domain:
  == domain 1  score: 48.2 bits;  conditional E-value: 2.8e-15
  2UVO:A|PDBID|CHAIN|SEQUENCE  83 cradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkda 136
                                  +ra+ +cg qagg  cpn lccs+wg+cgl  ++c ggcqs  c   +  g d 
            AFDB:AF-Q7Y1Z1-F1  30 ARAE-QCGRQAGGARCPNRLCCSRWGWCGLTDDYCKGGCQS-QCRVSRDGGDDD 81 
                                  5555.8**********************************8.588888777764 PP

  == domain 2  score: -2.2 bits;  conditional E-value: 7.6
  2UVO:A|PDBID|CHAIN|SEQUENCE  34 gcqngacwtskrcgsqagga 53 
                                  g+ +  c  + r     g a
            AFDB:AF-Q7Y1Z1-F1 179 GATSDFCVPNARWPCAPGKA 198
                                  33333333333333333333 PP

>> AFDB:AF-A0A1D6LMS5-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   44.5  16.7   3.9e-14   3.9e-10      82     127 ..      20      67 ..      12      82 .. 0.80

  Alignments for each domain:
  == domain 1  score: 44.5 bits;  conditional E-value: 3.9e-14
  2UVO:A|PDBID|CHAIN|SEQUENCE  82 pcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqs...gacs 127
                                  p+ra+ +cgsqagg lcpn lccsq+g+cg  s++cg+gcqs   g+c 
        AFDB:AF-A0A1D6LMS5-F1  20 PARAE-QCGSQAGGALCPNCLCCSQFGWCGSTSDYCGSGCQSqcsGSCG 67 
                                  88887.8**********************************72225553 PP

>> AFDB:AF-Q7DNA1-F1  Chitinase 2
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   43.0  15.2   1.1e-13   1.1e-09      84     128 ..      31      73 ..      20      91 .. 0.85
   2 ?   -1.4   1.3       4.1   4.2e+04      52      78 ..     288     314 ..     269     331 .. 0.66

  Alignments for each domain:
  == domain 1  score: 43.0 bits;  conditional E-value: 1.1e-13
  2UVO:A|PDBID|CHAIN|SEQUENCE  84 radikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacst 128
                                  ra+ +cg+qagg  cpn lccs+wg+cg  s+fcg gcqs  cs 
            AFDB:AF-Q7DNA1-F1  31 RAE-QCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ-CSG 73 
                                  555.8**********************************85.543 PP

  == domain 2  score: -1.4 bits;  conditional E-value: 4.1
  2UVO:A|PDBID|CHAIN|SEQUENCE  52 gatctnnqccsqygycgfgaeycgagc 78 
                                  g  c +       +  gf   ycga  
            AFDB:AF-Q7DNA1-F1 288 GLECGHGPDDRVANRIGFYQRYCGAFG 314
                                  444544444444555667777777643 PP

>> AFDB:AF-Q9SDY6-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   41.4  10.5   3.4e-13   3.4e-09      88     129 ..      25      65 ..      19      87 .. 0.75
   2 ?   -2.2   5.4       7.3   7.3e+04      49      86 ..     159     185 ..     121     211 .. 0.61
   3 ?   -2.3   0.2       8.2   8.2e+04      51      58 ..     277     284 ..     251     315 .. 0.61

  Alignments for each domain:
  == domain 1  score: 41.4 bits;  conditional E-value: 3.4e-13
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstd 129
                                  +cg+qagg lcpn lccs++g+cg    +cg gcqs   s  
            AFDB:AF-Q9SDY6-F1  25 QCGTQAGGALCPNRLCCSKFGWCGDTDSYCGEGCQSQCKS-A 65 
                                  7999999999999999999999999999999999985322.2 PP

  == domain 2  score: -2.2 bits;  conditional E-value: 7.3
  2UVO:A|PDBID|CHAIN|SEQUENCE  49 qaggatctnnqccsqygycgfgaeycgagcqgg..pcrad 86 
                                   a             +gyc ++ +  +  c gg  pc a 
            AFDB:AF-Q9SDY6-F1 159 YA-------------WGYCFINEQNQATYCDGGnwPCAAG 185
                                  33.............3444444444444444433334333 PP

  == domain 3  score: -2.3 bits;  conditional E-value: 8.2
  2UVO:A|PDBID|CHAIN|SEQUENCE  51 ggatctnn 58 
                                  gg  c + 
            AFDB:AF-Q9SDY6-F1 277 GGLECGHG 284
                                  33333333 PP

>> AFDB:AF-C6T7J9-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   39.8   8.0   1.1e-12   1.1e-08       7      44 ..      34      71 ..      29      88 .. 0.87
   2 ?   -2.0   0.3       6.3   6.3e+04      29      42 ..     159     172 ..     146     194 .. 0.59

  Alignments for each domain:
  == domain 1  score: 39.8 bits;  conditional E-value: 1.1e-12
  2UVO:A|PDBID|CHAIN|SEQUENCE  7 gsnmecpnnlccsqygycgmggdycgkgcqngacwtsk 44
                                 + n  c   lccs+ygycg g dycgkgc+ g c+ + 
            AFDB:AF-C6T7J9-F1 34 AQNCGCEAELCCSKYGYCGSGDDYCGKGCKEGPCYGTA 71
                                 57899*****************************9764 PP

  == domain 2  score: -2.0 bits;  conditional E-value: 6.3
  2UVO:A|PDBID|CHAIN|SEQUENCE  29 dycgkgcqngacwt 42 
                                  dyc k  ++  c  
            AFDB:AF-C6T7J9-F1 159 DYCDKTNRHYPCAH 172
                                  44444333333322 PP

>> AFDB:AF-I1MMY2-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   38.9  10.9     2e-12     2e-08      87     138 ..      24      75 ..      18      89 .. 0.76
   2 ?    0.7   1.0      0.94   9.4e+03     107     127 ..     158     178 ..     117     204 .. 0.67

  Alignments for each domain:
  == domain 1  score: 38.9 bits;  conditional E-value: 2e-12
  2UVO:A|PDBID|CHAIN|SEQUENCE  87 ikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdagg 138
                                    cg+q gg +cpn lccsq+g+cg     cg gcqs       p     +g
            AFDB:AF-I1MMY2-F1  24 QNCGTQVGGVICPNGLCCSQYGWCGNTEAHCGRGCQSQCTPGSTPTPTTPSG 75 
                                  47************************999*******9866555555443333 PP

  == domain 2  score: 0.7 bits;  conditional E-value: 0.94
  2UVO:A|PDBID|CHAIN|SEQUENCE 107 wgfcglgsefcgggcqsgacs 127
                                  wg+c ++    +  c sg   
            AFDB:AF-I1MMY2-F1 158 WGYCFINERNQADYCTSGTRW 178
                                  555555555444555544322 PP

>> AFDB:AF-Q6K8R2-F1  Chitinase 6
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   38.3   8.5     3e-12     3e-08      51      85 ..      26      60 ..      15      78 .. 0.73
   2 ?   -0.0   0.5       1.6   1.6e+04     144     163 ..     148     169 ..     140     194 .. 0.73

  Alignments for each domain:
  == domain 1  score: 38.3 bits;  conditional E-value: 3e-12
  2UVO:A|PDBID|CHAIN|SEQUENCE 51 ggatctnnqccsqygycgfgaeycgagcqggpcra 85
                                  +  c ++qccs++g+cg g++ycg gcq+gpc  
            AFDB:AF-Q6K8R2-F1 26 QSCGCASDQCCSKWGFCGTGSDYCGTGCQAGPCDV 60
                                 34557888888888888888888888888888854 PP

  == domain 2  score: -0.0 bits;  conditional E-value: 1.6
  2UVO:A|PDBID|CHAIN|SEQUENCE 144 nyc..cskwgscgigpgycgag 163
                                  nyc   s    c  g gy g g
            AFDB:AF-Q6K8R2-F1 148 NYCdeTSTQWPCMAGKGYYGRG 169
                                  6663322333577788887776 PP

>> AFDB:AF-B6TR38-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   36.3  15.4   1.2e-11   1.3e-07      43      84 ..      31      71 ..      23      89 .. 0.80

  Alignments for each domain:
  == domain 1  score: 36.3 bits;  conditional E-value: 1.2e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE 43 skrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcr 84
                                  ++cg qaggatc +  ccs++g+cg  +eycgagcq+  c 
            AFDB:AF-B6TR38-F1 31 GQQCGQQAGGATCRDCLCCSRFGFCGDTSEYCGAGCQS-QCT 71
                                 57899999999999999999999999999999999996.343 PP

>> AFDB:AF-B6TT00-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   36.0  10.1   1.5e-11   1.5e-07      51      84 ..      26      59 ..      19      74 .. 0.64
   2 ?    0.3   0.3       1.3   1.3e+04     143     163 ..     147     169 ..      99     192 .. 0.81

  Alignments for each domain:
  == domain 1  score: 36.0 bits;  conditional E-value: 1.5e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE 51 ggatctnnqccsqygycgfgaeycgagcqggpcr 84
                                  +  c +  ccs++gycg g +ycgagcq+gpc 
            AFDB:AF-B6TT00-F1 26 QNCGCASGLCCSRFGYCGTGEDYCGAGCQSGPCD 59
                                 3445666666666666666666666666666664 PP

  == domain 2  score: 0.3 bits;  conditional E-value: 1.3
  2UVO:A|PDBID|CHAIN|SEQUENCE 143 nnycc...skwgscgigpgycgag 163
                                   nyc    ++w  c  g gy g g
            AFDB:AF-B6TT00-F1 147 KNYCDrnnTQW-PCQAGKGYYGRG 169
                                  45554211344.577777777766 PP

>> AFDB:AF-Q42993-F1  Chitinase 1
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   35.9  22.6   1.7e-11   1.7e-07      88     138 ..      22      75 ..      16      83 .. 0.82
   2 ?    0.8   1.3      0.88   8.9e+03     143     158 ..     173     189 ..     142     221 .. 0.57

  Alignments for each domain:
  == domain 1  score: 35.9 bits;  conditional E-value: 1.7e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqs...gacstdkpcgkdagg 138
                                  +cgsqagg lcpn lccsq+g+cg  s +cg+gcqs   g+c    p     gg
            AFDB:AF-Q42993-F1  22 QCGSQAGGALCPNCLCCSQYGWCGSTSAYCGSGCQSqcsGSCGGGGPTPPSGGG 75 
                                  7*********************************96333777777666655555 PP

  == domain 2  score: 0.8 bits;  conditional E-value: 0.88
  2UVO:A|PDBID|CHAIN|SEQUENCE 143 nnycc..skwgscgigpg 158
                                  ++yc   s+w  c+ g  
            AFDB:AF-Q42993-F1 173 SDYCVqsSQW-PCAAGKK 189
                                  2222111112.2223333 PP

>> AFDB:AF-P25765-F1  Chitinase 12
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   35.7  21.9     2e-11     2e-07      88     133 ..      23      66 ..      15      82 .. 0.59

  Alignments for each domain:
  == domain 1  score: 35.7 bits;  conditional E-value: 2e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcg 133
                                  +cgsqagg +cpn lccsq+g+cg  s++cg+gcqs  cs    cg
            AFDB:AF-P25765-F1  23 QCGSQAGGAVCPNCLCCSQFGWCGSTSDYCGAGCQSQ-CSA-AGCG 66 
                                  4666666666666666666666666666666666653.333.1233 PP

>> AFDB:AF-O24603-F1  Endochitinase CHI
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   35.2   9.6   2.7e-11   2.7e-07      47      91 ..      30      72 ..      23     112 .. 0.78
   2 ?   -0.1   1.0       1.7   1.7e+04     143     164 ..     155     178 ..     137     196 .. 0.60
   3 ?   -2.1   0.1       6.9   6.9e+04     145     145 ..     259     259 ..     227     273 .. 0.50

  Alignments for each domain:
  == domain 1  score: 35.2 bits;  conditional E-value: 2.7e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE 47 gsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgs 91
                                 +sq     c ++ ccs+ygycg   e+cg gcq+gpcr+    g 
            AFDB:AF-O24603-F1 30 ASQN--CGCASDFCCSKYGYCGTTDEFCGEGCQAGPCRSSGGGGD 72
                                 5554..458999999999999999999999999999998765554 PP

  == domain 2  score: -0.1 bits;  conditional E-value: 1.7
  2UVO:A|PDBID|CHAIN|SEQUENCE 143 nnyccskw..gscgigpgycgagc 164
                                    yc ++     c+ g gy g g+
            AFDB:AF-O24603-F1 155 GEYCDTEKpeFPCAQGKGYYGRGA 178
                                  444433320123555555555543 PP

  == domain 3  score: -2.1 bits;  conditional E-value: 6.9
  2UVO:A|PDBID|CHAIN|SEQUENCE 145 y 145
                                  y
            AFDB:AF-O24603-F1 259 Y 259
                                  1 PP

>> AFDB:AF-P43082-F1  Hevein-like preproprotein
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   35.1  14.9   2.9e-11     3e-07       2      48 ..      23      70 ..      22      88 .. 0.80
   2 ?   -2.5   2.7       9.2   9.3e+04      23      51 ..     120     146 ..     105     207 .. 0.61

  Alignments for each domain:
  == domain 1  score: 35.1 bits;  conditional E-value: 2.9e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  2 rcgeqgsnmecpnnlccsqygycgmggdycg..kgcqngacwtskrcgs 48
                                 +cg qg    cp n+ccsqygycg  +dyc+   +cq+  cw s   g 
            AFDB:AF-P43082-F1 23 QCGRQGGGRTCPGNICCSQYGYCGTTADYCSptNNCQS-NCWGSGPSGP 70
                                 8*****************************83357987.59*9876654 PP

  == domain 2  score: -2.5 bits;  conditional E-value: 9.2
  2UVO:A|PDBID|CHAIN|SEQUENCE  23 ycgmggdycgkgcqngacwtskrcgsqag 51 
                                  +cg +g     +c  g c   k+  + a+
            AFDB:AF-P43082-F1 120 FCGPAGPRGQASC--GKCLRVKNTRTNAA 146
                                  2333333222222..45555555444444 PP

>> AFDB:AF-A0A1D6LMS1-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   34.5  19.1   4.5e-11   4.5e-07      88     135 ..      27      72 ..      20      81 .. 0.83

  Alignments for each domain:
  == domain 1  score: 34.5 bits;  conditional E-value: 4.5e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkd 135
                                  +cg+qagg lcp+ lccsqwg+cg   ++c  gcqs        cg  
        AFDB:AF-A0A1D6LMS1-F1  27 QCGTQAGGALCPDCLCCSQWGYCGSTPDYCTDGCQSQCFG--SGCGGG 72 
                                  7**********************************97544..457754 PP

>> AFDB:AF-A0A1D6GWN3-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   34.4  17.6   4.8e-11   4.8e-07      87     124 ..      22      59 ..      12      69 .. 0.88
   2 ?   -0.8   1.5       2.8   2.8e+04      22      37 ..     171     186 ..     136     199 .. 0.57

  Alignments for each domain:
  == domain 1  score: 34.4 bits;  conditional E-value: 4.8e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  87 ikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg 124
                                   +cg  a gklcp+ lccs+wg+cg  s++cg gcqs 
        AFDB:AF-A0A1D6GWN3-F1  22 AQCGDGADGKLCPDCLCCSKWGYCGSTSDYCGDGCQSQ 59 
                                  58**********************************95 PP

  == domain 2  score: -0.8 bits;  conditional E-value: 2.8
  2UVO:A|PDBID|CHAIN|SEQUENCE  22 gycgmggdycgkgcqn 37 
                                   yc m g+y+   c  
        AFDB:AF-A0A1D6GWN3-F1 171 DYCDMTGEYAQWPCVA 186
                                  3444444444333333 PP

>> AFDB:AF-O24598-F1  Endochitinase At2g43580
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   34.3  12.1     5e-11     5e-07       9      90 ..      26      61 ..      20      98 .. 0.53
   2 ?   -1.8   0.2       5.6   5.7e+04      25      33 ..     246     254 ..     225     264 .. 0.51

  Alignments for each domain:
  == domain 1  score: 34.3 bits;  conditional E-value: 5e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  9 nmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcg 90
                                 n +                                           c  n ccsq+gycg  a+ycg+ cq+gpcr     g
            AFDB:AF-O24598-F1 26 NCD-------------------------------------------CAPNLCCSQFGYCGTTADYCGSTCQSGPCRVG---G 61
                                 444...........................................55555555555555555555555555555542...2 PP

  == domain 2  score: -1.8 bits;  conditional E-value: 5.6
  2UVO:A|PDBID|CHAIN|SEQUENCE  25 gmggdycgk 33 
                                  g   dycg+
            AFDB:AF-O24598-F1 246 GYYRDYCGQ 254
                                  333344443 PP

>> AFDB:AF-Q9M2U5-F1  Endochitinase EP3
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   34.2   8.6   5.3e-11   5.3e-07       8      43 ..      29      64 ..      24     103 .. 0.70
   2 ?   -1.7   0.4       5.4   5.4e+04      34      83 ..     157     163 ..     137     193 .. 0.59

  Alignments for each domain:
  == domain 1  score: 34.2 bits;  conditional E-value: 5.3e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  8 snmecpnnlccsqygycgmggdycgkgcqngacwts 43
                                  n  c + lccsq+g+cg  +dycg gcq g c+  
            AFDB:AF-Q9M2U5-F1 29 QNCGCSSELCCSQFGFCGNTSDYCGVGCQQGPCFAP 64
                                 455666666666666666666666666666666654 PP

  == domain 2  score: -1.7 bits;  conditional E-value: 5.4
  2UVO:A|PDBID|CHAIN|SEQUENCE  34 gcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpc 83 
                                  ++                                              pc
            AFDB:AF-Q9M2U5-F1 157 NA-------------------------------------------TQYPC 163
                                  22...........................................22222 PP

>> AFDB:AF-I1M587-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   33.9  11.4   6.8e-11   6.8e-07       8      46 ..      28      66 ..      23      94 .. 0.77
   2 ?   -1.4   0.3       4.3   4.4e+04      57      70 ..     161     174 ..     137     195 .. 0.67

  Alignments for each domain:
  == domain 1  score: 33.9 bits;  conditional E-value: 6.8e-11
  2UVO:A|PDBID|CHAIN|SEQUENCE  8 snmecpnnlccsqygycgmggdycgkgcqngacwtskrc 46
                                  n  c   lccsq+gycg g +ycg gc+ g c++s   
            AFDB:AF-I1M587-F1 28 QNCGCAEGLCCSQHGYCGNGEEYCGTGCKQGPCYSSTPS 66
                                 577788888888888888888888888888888887655 PP

  == domain 2  score: -1.4 bits;  conditional E-value: 4.3
  2UVO:A|PDBID|CHAIN|SEQUENCE  57 nnqccsqygycgfg 70 
                                     c s  gy g g
            AFDB:AF-I1M587-F1 161 QYPCLSNRGYYGRG 174
                                  33444445554444 PP

>> AFDB:AF-O22841-F1  Endochitinase At2g43620
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   32.4  10.9     2e-10     2e-06      48      86 ..      29      67 ..      21      86 .. 0.76
   2 ?   -2.0   0.8       6.4   6.5e+04      20      34 ..     169     184 ..     152     198 .. 0.66

  Alignments for each domain:
  == domain 1  score: 32.4 bits;  conditional E-value: 2e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE 48 sqaggatctnnqccsqygycgfgaeycgagcqggpcrad 86
                                  q g   c  n ccs+ygycg    ycg gc++gpc + 
            AFDB:AF-O22841-F1 29 QQCGTTGCAANLCCSRYGYCGTTDAYCGTGCRSGPCSSS 67
                                 345556688888888888888888888888888888754 PP

  == domain 2  score: -2.0 bits;  conditional E-value: 6.4
  2UVO:A|PDBID|CHAIN|SEQUENCE  20 qygy.cgmggdycgkg 34 
                                    +y c  g dy g+g
            AFDB:AF-O22841-F1 169 STAYpCTPGKDYYGRG 184
                                  3333366666777666 PP

>> AFDB:AF-O22842-F1  Endochitinase At2g43610
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   32.2  13.0   2.3e-10   2.3e-06      93     128 ..      31      66 ..      21      90 .. 0.59

  Alignments for each domain:
  == domain 1  score: 32.2 bits;  conditional E-value: 2.3e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE  93 aggklcpnnlccsqwgfcglgsefcgggcqsgacst 128
                                   g   c  n+ccs+wg+cg    +cg gcqsg c++
            AFDB:AF-O22842-F1  31 CGTNGCKGNMCCSRWGYCGTTKAYCGTGCQSGPCNS 66 
                                  333445555555555555555555555555555543 PP

>> AFDB:AF-A0A1D6GWN1-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   31.5  11.0   3.6e-10   3.6e-06      87     124 ..      27      64 ..      15      74 .. 0.87
   2 ?   -1.1   0.3       3.4   3.4e+04     129     168 ..     225     236 ..     202     271 .. 0.69

  Alignments for each domain:
  == domain 1  score: 31.5 bits;  conditional E-value: 3.6e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE  87 ikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsg 124
                                   +cg+ +   lcp  lccs+wgfcg    +cg+gcqs 
        AFDB:AF-A0A1D6GWN1-F1  27 PQCGANSTTALCPYCLCCSKWGFCGSTEAYCGNGCQSQ 64 
                                  47**********************************95 PP

  == domain 2  score: -1.1 bits;  conditional E-value: 3.4
  2UVO:A|PDBID|CHAIN|SEQUENCE 129 dkpcgkdaggrvctnnyccskwgscgigpgycgagcqsgg 168
                                  d  c                            g g  +gg
        AFDB:AF-A0A1D6GWN1-F1 225 DAEC----------------------------GRGPDAGG 236
                                  3333............................22222222 PP

>> AFDB:AF-P24626-F1  Chitinase 3
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   31.1  27.6     5e-10     5e-06      88     127 ..      20      62 ..      14      81 .. 0.74

  Alignments for each domain:
  == domain 1  score: 31.1 bits;  conditional E-value: 5e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqs...gacs 127
                                  +cgsqagg lcpn lccsq+g+cg  s++cg+gcqs   g c 
            AFDB:AF-P24626-F1  20 QCGSQAGGALCPNCLCCSQYGWCGSTSDYCGAGCQSqcsGGCG 62 
                                  6999999999999999999999999999999999862224443 PP

>> AFDB:AF-O24658-F1  Endochitinase At2g43590
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   30.8  14.1   6.2e-10   6.2e-06      53      97 ..      27      68 ..      17      97 .. 0.67
   2 ?    0.5   1.6       1.1   1.1e+04     141     163 ..     141     165 ..     135     183 .. 0.73
   3 ?   -1.0   0.6       3.1   3.1e+04      54      76 ..     231     253 ..     223     263 .. 0.52

  Alignments for each domain:
  == domain 1  score: 30.8 bits;  conditional E-value: 6.2e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE 53 atctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggkl 97
                                   c  n ccsq+gycg    ycg gc++gpcr     g+  gg +
            AFDB:AF-O24658-F1 27 CGCAPNLCCSQFGYCGTDDAYCGVGCRSGPCRGS---GTPTGGSV 68
                                 4577888888888888888888888888888864...44444443 PP

  == domain 2  score: 0.5 bits;  conditional E-value: 1.1
  2UVO:A|PDBID|CHAIN|SEQUENCE 141 ctnnyccsk..wgscgigpgycgag 163
                                  +t nyc s      c+ g gy g g
            AFDB:AF-O24658-F1 141 ATRNYCQSSntQYPCAPGKGYFGRG 165
                                  5778887651134688888888876 PP

  == domain 3  score: -1.0 bits;  conditional E-value: 3.1
  2UVO:A|PDBID|CHAIN|SEQUENCE  54 tctnnqccsqygycgfgaeycga 76 
                                   c      +  +  g+  +ycg 
            AFDB:AF-O24658-F1 231 ECNGGNSGAVNARIGYYRDYCGQ 253
                                  44444444444445555555553 PP

>> AFDB:AF-Q0JC38-F1  Os04g0493600 protein
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   30.1  17.6   9.6e-10   9.7e-06       8      47 ..      28      65 ..      22      82 .. 0.80

  Alignments for each domain:
  == domain 1  score: 30.1 bits;  conditional E-value: 9.6e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE  8 snmecpnnlccsqygycgmggdycgkgcqngacwtskrcg 47
                                  n  c +  ccsq+gycg    ycg+gcq+g cw s   g
            AFDB:AF-Q0JC38-F1 28 QNCGCQDGYCCSQWGYCGTTEAYCGQGCQSGPCWGSG--G 65
                                 6888999999999999999999999999999999874..2 PP

>> AFDB:AF-P19171-F1  Basic endochitinase B
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   30.4  12.2     8e-10     8e-06      88     133 ..      35      81 ..      30      96 .. 0.80
   2 ?   -2.6   0.7       9.6   9.6e+04     130     142 ..     286     298 ..     268     321 .. 0.58

  Alignments for each domain:
  == domain 1  score: 30.4 bits;  conditional E-value: 8e-10
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcg.ggcqsgacstdkpcg 133
                                  +cg qagg lcpn lccs++g+cg    +c   gcqs       p g
            AFDB:AF-P19171-F1  35 QCGRQAGGALCPNGLCCSEFGWCGNTEPYCKqPGCQSQCTPGGTPPG 81 
                                  7*****************************6369*997666666665 PP

  == domain 2  score: -2.6 bits;  conditional E-value: 9.6
  2UVO:A|PDBID|CHAIN|SEQUENCE 130 kpcgkdaggrvct 142
                                    cg+   grv+ 
            AFDB:AF-P19171-F1 286 LECGRGQDGRVAD 298
                                  3444444444443 PP

>> AFDB:AF-A0A1X7YIJ7-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   28.8  21.2   2.4e-09   2.4e-05      10      80 ..      34      99 ..      28     104 .. 0.76

  Alignments for each domain:
  == domain 1  score: 28.8 bits;  conditional E-value: 2.4e-09
  2UVO:A|PDBID|CHAIN|SEQUENCE 10 mecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqg 80
                                   c   +ccs+ygycg  + ycg+gc++g cw s  cg   gga+       ++  + g+   ++g+ c+g
        AFDB:AF-A0A1X7YIJ7-F1 34 CGCQPGFCCSKYGYCGKTSAYCGEGCKSGPCWGSAGCGG--GGASVARV--VTKSFFNGIK-SHAGSWCEG 99
                                 568899*******************************95..77776543..3333344443.456666665 PP

>> AFDB:AF-O04138-F1  Chitinase 4
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   29.4  15.2   1.6e-09   1.6e-05       8      44 ..      28      64 ..      22      88 .. 0.82
   2 ?   -0.6   0.4       2.4   2.5e+04      25      45 ..     160     180 ..     146     200 .. 0.70

  Alignments for each domain:
  == domain 1  score: 29.4 bits;  conditional E-value: 1.6e-09
  2UVO:A|PDBID|CHAIN|SEQUENCE  8 snmecpnnlccsqygycgmggdycgkgcqngacwtsk 44
                                  n  c +  ccsq+gycg    ycg+gcq+g cw s 
            AFDB:AF-O04138-F1 28 QNCGCQDGYCCSQWGYCGTTEAYCGQGCQSGPCWGSG 64
                                 6888999999999999999999999999999999874 PP

  == domain 2  score: -0.6 bits;  conditional E-value: 2.4
  2UVO:A|PDBID|CHAIN|SEQUENCE  25 gmggdycgkgcqngacwtskr 45 
                                  g + dyc k+ +   c   k+
            AFDB:AF-O04138-F1 160 GANMDYCDKSNKQWPCQPGKK 180
                                  555567776666666665554 PP

>> AFDB:AF-P29023-F1  Endochitinase B
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   27.6  25.6   5.7e-09   5.8e-05       9     123 ..      22      91 ..      12      99 .. 0.50
   2 ?   -2.1   0.3       6.8   6.8e+04     153     163 ..     159     169 ..     143     193 .. 0.59

  Alignments for each domain:
  == domain 1  score: 27.6 bits;  conditional E-value: 5.7e-09
  2UVO:A|PDBID|CHAIN|SEQUENCE   9 nmecpnnlccsqygycgmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqs 123
                                  n  c  n+ccs++gycg   +ycg gcq+g c                                           r+    g   gg     ++  s + f g+ ++ +g+gc+ 
            AFDB:AF-P29023-F1  22 NCGCQPNVCCSKFGYCGTTDEYCGDGCQSGPC-------------------------------------------RSGRGGGGSGGGGANVASVVTSSF-FNGIKNQ-AGSGCEG 91 
                                  44455555555555555555555555555555...........................................555444444444444334433332.4444433.3444443 PP

  == domain 2  score: -2.1 bits;  conditional E-value: 6.8
  2UVO:A|PDBID|CHAIN|SEQUENCE 153 cgigpgycgag 163
                                  c+ g  y g g
            AFDB:AF-P29023-F1 159 CAAGQKYYGRG 169
                                  33333333322 PP

>> AFDB:AF-C0P451-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   26.5  24.9   1.2e-08   0.00012      51     123 ..      34     103 ..      23     111 .. 0.59
   2 ?   -2.3   0.4       8.2   8.2e+04     153     163 ..     171     181 ..     156     204 .. 0.59
   3 ?   -2.4   0.3       8.8   8.8e+04      51      75 ..     245     269 ..     238     278 .. 0.64

  Alignments for each domain:
  == domain 1  score: 26.5 bits;  conditional E-value: 1.2e-08
  2UVO:A|PDBID|CHAIN|SEQUENCE  51 ggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqs 123
                                   +  c  n ccs++gycg   eycg gcq+gpcr+    gs  gg     ++      f g+ s+ +g+gc+ 
            AFDB:AF-C0P451-F1  34 QNCGCQPNVCCSKFGYCGTTDEYCGDGCQSGPCRSGG-GGSSGGGGANVASVVTGS-FFNGIKSQ-AGSGCEG 103
                                  3445677777777777777777777777777777654.344444443333333332.35566555.5666654 PP

  == domain 2  score: -2.3 bits;  conditional E-value: 8.2
  2UVO:A|PDBID|CHAIN|SEQUENCE 153 cgigpgycgag 163
                                  c+ g  y g g
            AFDB:AF-C0P451-F1 171 CAAGQKYYGRG 181
                                  33333333322 PP

  == domain 3  score: -2.4 bits;  conditional E-value: 8.8
  2UVO:A|PDBID|CHAIN|SEQUENCE  51 ggatctnnqccsqygycgfgaeycg 75 
                                  g+  c  n  +   +  g+  +yc 
            AFDB:AF-C0P451-F1 245 GALECGGNNPAQMNARVGYYRQYCR 269
                                  4456777777777777777777774 PP

>> AFDB:AF-I1NCA0-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   25.0  27.8   3.6e-08   0.00036       1      62 [.      21      83 ..      21     103 .. 0.81
   2 ?   -1.8   0.7       5.6   5.6e+04      27      40 ..     131     138 ..     113     163 .. 0.48

  Alignments for each domain:
  == domain 1  score: 25.0 bits;  conditional E-value: 3.6e-08
  2UVO:A|PDBID|CHAIN|SEQUENCE  1 ercgeqgsnmecpnnlccsqygycgmggdycg..kgcqngacwtskrcgsqaggatctnnqccs 62
                                 e+cg q+  + cpnnlccsqyg+cg   +yc+  k+cq++ cw     g   gg    +n  ++
            AFDB:AF-I1NCA0-F1 21 EQCGRQAGGQTCPNNLCCSQYGWCGNTEEYCSpsKNCQSN-CWGGGGGGGGGGGGESASNVRAT 83
                                 78*****************************544899975.99998888877777666664443 PP

  == domain 2  score: -1.8 bits;  conditional E-value: 5.6
  2UVO:A|PDBID|CHAIN|SEQUENCE  27 ggdycgkgcqngac 40 
                                  g d c      g c
            AFDB:AF-I1NCA0-F1 131 GRDSC------GKC 138
                                  22222......333 PP

>> AFDB:AF-A0A0P0Y930-F1  Os12g0238550 protein
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   23.2   0.3   1.3e-07    0.0013      65      85 ..       8      28 ..       2      45 .. 0.74

  Alignments for each domain:
  == domain 1  score: 23.2 bits;  conditional E-value: 1.3e-07
  2UVO:A|PDBID|CHAIN|SEQUENCE 65 gycgfgaeycgagcqggpcra 85
                                 ++cg g++y g gcq+gpc  
        AFDB:AF-A0A0P0Y930-F1  8 AFCGTGSDYYGTGCQAGPCDV 28
                                 678888888888888888854 PP

>> AFDB:AF-O24654-F1  Inactive endochitinase At2g43600
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   23.0   7.4   1.4e-07    0.0014      98     128 ..      30      61 ..      20      77 .. 0.76
   2 ?    1.0   0.4      0.77   7.7e+03      14      42 ..     152     182 ..     145     196 .. 0.67

  Alignments for each domain:
  == domain 1  score: 23.0 bits;  conditional E-value: 1.4e-07
  2UVO:A|PDBID|CHAIN|SEQUENCE  98 cpn.nlccsqwgfcglgsefcgggcqsgacst 128
                                  cp    ccs+wgfcg   e+cg  c sg c+ 
            AFDB:AF-O24654-F1  30 CPGlKECCSRWGFCGTKDEYCGFFCFSGPCNI 61 
                                  66534699999999999999999999999975 PP

  == domain 2  score: 1.0 bits;  conditional E-value: 0.77
  2UVO:A|PDBID|CHAIN|SEQUENCE  14 nnlccsq.ygy.cgmggdycgkgcqngacwt 42 
                                  n   cs+   y c  g +y g+g   +  w 
            AFDB:AF-O24654-F1 152 NERYCSKsKKYpCEPGKNYYGRGLLQSITWN 182
                                  4444443123338888888888888888886 PP

>> AFDB:AF-A0A1P8AME8-F1  Chitinase family protein
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   23.0   6.8   1.5e-07    0.0015      57      86 ..      62      91 ..      47     101 .. 0.79
   2 ?   -0.1   0.7       1.7   1.7e+04      16      37 ..     182     204 ..     168     217 .. 0.57

  Alignments for each domain:
  == domain 1  score: 23.0 bits;  conditional E-value: 1.5e-07
  2UVO:A|PDBID|CHAIN|SEQUENCE 57 nnqccsqygycgfgaeycgagcqggpcrad 86
                                  n+ccs  gycg + e+cg  c +gpc+  
        AFDB:AF-A0A1P8AME8-F1 62 INECCSHTGYCGTNVEHCGFWCLSGPCQLS 91
                                 489999999999999999999999999865 PP

  == domain 2  score: -0.1 bits;  conditional E-value: 1.7
  2UVO:A|PDBID|CHAIN|SEQUENCE  16 lccsqygy.cgmggdycgkgcqn 37 
                                   c s   y c  g  y g+g   
        AFDB:AF-A0A1P8AME8-F1 182 YCSSSKTYpCQSGKKYYGRGLLQ 204
                                  23333333244444555555444 PP

>> AFDB:AF-Q688M5-F1  Chitinase 9
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   22.4  33.5   2.2e-07    0.0022      88     138 ..      25      71 ..      12      85 .. 0.79
   2 ?   -0.9   0.4         3     3e+04      64      86 ..     166     192 ..     142     208 .. 0.56

  Alignments for each domain:
  == domain 1  score: 22.4 bits;  conditional E-value: 2.2e-07
  2UVO:A|PDBID|CHAIN|SEQUENCE  88 kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdagg 138
                                  +cgsqagg lcpn lccs +g+cg  s++cg gcqs  c     cg   gg
            AFDB:AF-Q688M5-F1  25 QCGSQAGGALCPNCLCCSSYGWCGSTSDYCGDGCQSQ-CD---GCGGGGGG 71 
                                  7999999999999999999999999999999999985.43...24443333 PP

  == domain 2  score: -0.9 bits;  conditional E-value: 3
  2UVO:A|PDBID|CHAIN|SEQUENCE  64 ygyc.....gfgaeycgagcqggpcrad 86 
                                  +gyc     g  a yc  +++  pc  d
            AFDB:AF-Q688M5-F1 166 WGYCfkeeiGATASYCVPSAE-WPCAPD 192
                                  333322222333444443333.344444 PP

>> AFDB:AF-I1NCA1-F1  Chitinase
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   22.4  25.4   2.3e-07    0.0023       2      63 ..      23      82 ..      22     109 .. 0.76
   2 ?   -1.6   1.1       4.8   4.8e+04     124     139 ..     120     135 ..      85     145 .. 0.59

  Alignments for each domain:
  == domain 1  score: 22.4 bits;  conditional E-value: 2.3e-07
  2UVO:A|PDBID|CHAIN|SEQUENCE  2 rcgeqgsnmecpnnlccsqygycgmggdycg..kgcqngacwtskrcgsqaggatctnnqccsq 63
                                 +cg q+  + c nnlccsqyg+cg + d+c+  k+cq+  cw s       gg +++n   ++ 
            AFDB:AF-I1NCA1-F1 23 NCGRQAGGQTCGNNLCCSQYGWCGNSEDHCSpsKNCQS-TCWGSGG--GGGGGESASN-VRATY 82
                                 7*****************************64489996.8**9853..3345555444.33334 PP

  == domain 2  score: -1.6 bits;  conditional E-value: 4.8
  2UVO:A|PDBID|CHAIN|SEQUENCE 124 gacstdkpcgkdaggr 139
                                  + c    p g+da g+
            AFDB:AF-I1NCA1-F1 120 AFCGPVGPRGRDACGK 135
                                  3455555555555554 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (171 nodes)
Target sequences:                     362094  (155907410 residues searched)
Passed MSV filter:                     20812  (0.0574768); expected 7241.9 (0.02)
Passed bias filter:                     5905  (0.0163079); expected 7241.9 (0.02)
Passed Vit filter:                       550  (0.00151894); expected 362.1 (0.001)
Passed Fwd filter:                        38  (0.000104945); expected 3.6 (1e-05)
Initial search space (Z):             362094  [actual number of targets]
Domain search space  (domZ):              36  [number of targets reported over threshold]
# CPU time: 2.44u 0.10s 00:00:02.54 Elapsed: 00:00:00.88
# Mc/sec: 30188.23
//
# Alignment of 36 hits satisfying inclusion thresholds saved to: phmmerAlignment_af2.log
[ok]
"""