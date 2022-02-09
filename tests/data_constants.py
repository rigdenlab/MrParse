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

PHMMER_AF_JSON = '{"results":{"hits":[{"archScore":"8","ph":"Streptophyta","arch":"PF00187.22 PF00187.22 PF00187.22 P' \
                 'F00187.22","kg":"Eukaryota","ndom":1,"extlink":"https://alphafold.ebi.ac.uk/entry/Q0JF21","acc2":"Q' \
                 '0JF21","taxid":"39947","acc":"Q0JF21","taxlink":"http://www.uniprot.org/taxonomy/","desc":"Lectin",' \
                 '"pvalue":-181.501522091517,"flags":"3","nregions":1,"niseqs":16,"name":"Q0JF21","species":"Oryza sa' \
                 'tiva subsp. japonica","score":"253.2","bias":"98.2","sindex":"236295369","nincluded":"1","domains":' \
                 '[{"aliniseqs":1,"alisqacc":"F1","aliIdCount":125,"alirfline":"","is_included":1,"alihmmname":">Seq"' \
                 ',"bitscore":252.817687988281,"ievalue":"7.1e-74","alisqto":197,"aliSim":0.875,"jali":197,"bias":"98' \
                 '.19","ienv":29,"cevalue":"7.0e-78","alimline":" cg+q+  m cp+nlccsq+gycg+g dycg gcq+gac +s+rcgsq gga' \
                 'tc+nnqccsqygycgfg+eycg+gcq+gpcradikcg  a+g+lcpnn+ccsqwg+cglgsefcg+gcqsgac  +k cgk agg  c nn+ccs  g ' \
                 'cg+g +ycg+gcqsggc","alihmmfrom":2,"aliL":227,"is_reported":1,"alintseq":"","alisindex":"236295369",' \
                 '"jenv":201,"alimmline":"","alihmmacc":"","oasc":"0.99","aliaseq":"TCGKQNDGMICPHNLCCSQFGYCGLGRDYCGTG' \
                 'CQSGACCSSQRCGSQGGGATCSNNQCCSQYGYCGFGSEYCGSGCQNGPCRADIKCGRNANGELCPNNMCCSQWGYCGLGSEFCGNGCQSGACCPEKRCG' \
                 'KQAGGDKCPNNFCCSAGGYCGLGGNYCGSGCQSGGC","alihmmto":169,"aliId":0.744047619047619,"alippline":"6******' \
                 '***************************************************************************************************' \
                 '**************************************************************","alimodel":"rcgeqgsnmecpnnlccsqygyc' \
                 'gmggdycgkgcqngacwtskrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqs' \
                 'gacstdkpcgkdaggrvctnnyccskwgscgigpgycgagcqsggc","aliM":171,"iali":30,"alicsline":"","aliSimCount":1' \
                 '47,"alihmmdesc":"","alisqdesc":"Lectin","alisqname":"Q0JF21","alisqfrom":30,"aliN":168}],"evalue":"' \
                 '5.4e-74","nreported":1,"archindex":"88390276629085 39947"},{"archScore":"4","ph":"Streptophyta","ar' \
                 'ch":"PF00187.22 PF00182.22","kg":"Eukaryota","ndom":2,"extlink":"https://alphafold.ebi.ac.uk/entry/' \
                 'Q7Y1Z1","acc2":"Q7Y1Z1","taxid":"39947","acc":"Q7Y1Z1","taxlink":"http://www.uniprot.org/taxonomy/"' \
                 ',"desc":"Chitinase 7","pvalue":-37.0974407605701,"flags":"3","nregions":2,"niseqs":16,"name":"Q7Y1Z' \
                 '1","species":"Oryza sativa subsp. japonica","score":"48.2","bias":"10.0","sindex":"236343183","ninc' \
                 'luded":"1","domains":[{"aliniseqs":1,"alisqacc":"F1","aliIdCount":29,"alirfline":"","is_included":1' \
                 ',"alihmmname":">Seq","bitscore":48.2214241027832,"ievalue":"2.8e-11","alisqto":81,"aliSim":0.711538' \
                 '461538462,"jali":81,"bias":"9.96","ienv":20,"cevalue":"2.8e-15","alimline":"+ra+ +cg qagg  cpn lccs' \
                 '+wg+cgl  ++c ggcqs  c   +  g d ","alihmmfrom":83,"aliL":340,"is_reported":1,"alintseq":"","alisinde' \
                 'x":"236343183","jenv":97,"alimmline":"","alihmmacc":"","oasc":"0.86","aliaseq":"ARAE-QCGRQAGGARCPNR' \
                 'LCCSRWGWCGLTDDYCKGGCQS-QCRVSRDGGDDD","alihmmto":136,"aliId":0.557692307692308,"alippline":"5555.8**' \
                 '********************************8.588888777764","alimodel":"cradikcgsqaggklcpnnlccsqwgfcglgsefcgggc' \
                 'qsgacstdkpcgkda","aliM":171,"iali":30,"alicsline":"","aliSimCount":37,"alihmmdesc":"","alisqdesc":"' \
                 'Chitinase 7","alisqname":"Q7Y1Z1","alisqfrom":30,"aliN":54},{"aliniseqs":1,"alisqacc":"F1","aliIdCo' \
                 'unt":5,"alirfline":"","is_included":0,"alihmmname":">Seq","bitscore":-2.23642659187317,"ievalue":"7' \
                 '7000","alisqto":198,"aliSim":0.4,"jali":198,"bias":"1.04","ienv":159,"cevalue":"7.6","alimline":"g+' \
                 ' +  c  + r     g a","alihmmfrom":34,"aliL":340,"is_reported":0,"alintseq":"","alisindex":"236343183' \
                 '","jenv":230,"alimmline":"","alihmmacc":"","oasc":"0.48","aliaseq":"GATSDFCVPNARWPCAPGKA","alihmmto' \
                 '":53,"aliId":0.25,"alippline":"33333333333333333333","alimodel":"gcqngacwtskrcgsqagga","aliM":171,"' \
                 'iali":179,"alicsline":"","aliSimCount":8,"alihmmdesc":"","alisqdesc":"Chitinase 7","alisqname":"Q7Y' \
                 '1Z1","alisqfrom":179,"aliN":20}],"evalue":"2.8e-11","nreported":1,"archindex":"55525311637573 39947' \
                 '"},{"archScore":"4","ph":"Streptophyta","arch":"PF00187.22 PF00182.22","kg":"Eukaryota","ndom":3,"e' \
                 'xtlink":"https://alphafold.ebi.ac.uk/entry/A0A1D6LMS5","acc2":"A0A1D6LMS5","taxid":"4577","acc":"A0' \
                 'A1D6LMS5","taxlink":"http://www.uniprot.org/taxonomy/","desc":"Chitinase","pvalue":-34.464617586289' \
                 '1,"flags":"3","nregions":3,"niseqs":16,"name":"A0A1D6LMS5","species":"Zea mays","score":"44.5","bia' \
                 's":"16.7","sindex":"237720498","nincluded":"1","domains":[{"aliniseqs":1,"alisqacc":"F1","aliIdCoun' \
                 't":31,"alirfline":"","is_included":1,"alihmmname":">Seq","bitscore":44.4842414855957,"ievalue":"3.9' \
                 'e-10","alisqto":67,"aliSim":0.869565217391304,"jali":67,"bias":"16.70","ienv":12,"cevalue":"3.9e-14' \
                 '","alimline":"p+ra+ +cgsqagg lcpn lccsq+g+cg  s++cg+gcqs   g+c ","alihmmfrom":82,"aliL":310,"is_rep' \
                 'orted":1,"alintseq":"","alisindex":"237720498","jenv":82,"alimmline":"","alihmmacc":"","oasc":"0.80' \
                 '","aliaseq":"PARAE-QCGSQAGGALCPNCLCCSQFGWCGSTSDYCGSGCQSqcsGSCG","alihmmto":127,"aliId":0.6739130434' \
                 '78261,"alippline":"88887.8**********************************72225553","alimodel":"pcradikcgsqaggklc' \
                 'pnnlccsqwgfcglgsefcgggcqs...gacs","aliM":171,"iali":20,"alicsline":"","aliSimCount":40,"alihmmdesc"' \
                 ':"","alisqdesc":"Chitinase","alisqname":"A0A1D6LMS5","alisqfrom":20,"aliN":49},{"aliniseqs":1,"alis' \
                 'qacc":"F1","aliIdCount":17,"alirfline":"","is_included":0,"alihmmname":">Seq","bitscore":-4.5749888' \
                 '420105,"ievalue":"360000","alisqto":167,"aliSim":0.449275362318841,"jali":167,"bias":"10.49","ienv"' \
                 ':92,"cevalue":"36","alimline":"c  n   +  g+ + +  + g g  g+p  +    g   g   +p+      wg+c            ' \
                 '                          +    + gp yc  ++","alihmmfrom":55,"aliL":310,"is_reported":0,"alintseq":"' \
                 '","alisindex":"237720498","jenv":205,"alimmline":"","alihmmacc":"","oasc":"0.47","aliaseq":"CPANGFY' \
                 'TYAGFIAAANAFPGFGTTGAPDTSHETTG---GWATAPDGP--YAWGYCF------------------------------------KEEQGGASGPDYC' \
                 'EPSA","alihmmto":164,"aliId":0.246376811594203,"alippline":"333333333333333333444444444433332222...' \
                 '223334332..2344443....................................33333333333333322","alimodel":"ctnnqccsqygycg' \
                 'fgaeycgagcqggpcradikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdaggrvctnnyccskwgscgigpgycgagc","' \
                 'aliM":171,"iali":99,"alicsline":"","aliSimCount":31,"alihmmdesc":"","alisqdesc":"Chitinase","alisqn' \
                 'ame":"A0A1D6LMS5","alisqfrom":99,"aliN":110},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":0,"alirfli' \
                 'ne":"","is_included":0,"alihmmname":">Seq","bitscore":-2.89088797569275,"ievalue":"120000","alisqto' \
                 '":275,"aliSim":0,"jali":275,"bias":"0.18","ienv":251,"cevalue":"12","alimline":" ","alihmmfrom":137' \
                 ',"aliL":310,"is_reported":0,"alintseq":"","alisindex":"237720498","jenv":304,"alimmline":"","alihmm' \
                 'acc":"","oasc":"0.62","aliaseq":"D","alihmmto":137,"aliId":0,"alippline":"1","alimodel":"g","aliM":' \
                 '171,"iali":275,"alicsline":"","aliSimCount":0,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"' \
                 'A0A1D6LMS5","alisqfrom":275,"aliN":1}],"evalue":"3.9e-10","nreported":1,"archindex":"55525311637573' \
                 ' 4577"},{"archScore":"4","ph":"Streptophyta","arch":"PF00187.22 PF00182.22","kg":"Eukaryota","ndom"' \
                 ':3,"extlink":"https://alphafold.ebi.ac.uk/entry/Q7DNA1","acc2":"Q7DNA1","taxid":"39947","acc":"Q7DN' \
                 'A1","taxlink":"http://www.uniprot.org/taxonomy/","desc":"Chitinase 2","pvalue":-33.4454090806076,"f' \
                 'lags":"3","nregions":3,"niseqs":16,"name":"Q7DNA1","species":"Oryza sativa subsp. japonica","score"' \
                 ':"43.0","bias":"15.2","sindex":"236237656","nincluded":"1","domains":[{"aliniseqs":1,"alisqacc":"F1' \
                 '","aliIdCount":29,"alirfline":"","is_included":1,"alihmmname":">Seq","bitscore":43.0375175476074,"i' \
                 'evalue":"1.1e-09","alisqto":73,"aliSim":0.813953488372093,"jali":73,"bias":"15.22","ienv":20,"ceval' \
                 'ue":"1.1e-13","alimline":"ra+ +cg+qagg  cpn lccs+wg+cg  s+fcg gcqs  cs ","alihmmfrom":84,"aliL":340' \
                 ',"is_reported":1,"alintseq":"","alisindex":"236237656","jenv":91,"alimmline":"","alihmmacc":"","oas' \
                 'c":"0.85","aliaseq":"RAE-QCGAQAGGARCPNCLCCSRWGWCGTTSDFCGDGCQSQ-CSG","alihmmto":128,"aliId":0.674418' \
                 '604651163,"alippline":"555.8**********************************85.543","alimodel":"radikcgsqaggklcpn' \
                 'nlccsqwgfcglgsefcgggcqsgacst","aliM":171,"iali":31,"alicsline":"","aliSimCount":35,"alihmmdesc":"",' \
                 '"alisqdesc":"Chitinase 2","alisqname":"Q7DNA1","alisqfrom":31,"aliN":45},{"aliniseqs":1,"alisqacc":' \
                 '"F1","aliIdCount":1,"alirfline":"","is_included":0,"alihmmname":">Seq","bitscore":-3.32551097869873' \
                 ',"ievalue":"170000","alisqto":184,"aliSim":1,"jali":184,"bias":"0.28","ienv":168,"cevalue":"16","al' \
                 'imline":"c","alihmmfrom":17,"aliL":340,"is_reported":0,"alintseq":"","alisindex":"236237656","jenv"' \
                 ':205,"alimmline":"","alihmmacc":"","oasc":"0.46","aliaseq":"C","alihmmto":17,"aliId":1,"alippline":' \
                 '"2","alimodel":"c","aliM":171,"iali":184,"alicsline":"","aliSimCount":1,"alihmmdesc":"","alisqdesc"' \
                 ':"Chitinase 2","alisqname":"Q7DNA1","alisqfrom":184,"aliN":1},{"aliniseqs":1,"alisqacc":"F1","aliId' \
                 'Count":8,"alirfline":"","is_included":0,"alihmmname":">Seq","bitscore":-1.36343097686768,"ievalue":' \
                 '"42000","alisqto":314,"aliSim":0.37037037037037,"jali":314,"bias":"1.32","ienv":269,"cevalue":"4.1"' \
                 ',"alimline":"g  c +       +  gf   ycga  ","alihmmfrom":52,"aliL":340,"is_reported":0,"alintseq":"",' \
                 '"alisindex":"236237656","jenv":331,"alimmline":"","alihmmacc":"","oasc":"0.66","aliaseq":"GLECGHGPD' \
                 'DRVANRIGFYQRYCGAFG","alihmmto":78,"aliId":0.296296296296296,"alippline":"44454444444455566777777764' \
                 '3","alimodel":"gatctnnqccsqygycgfgaeycgagc","aliM":171,"iali":288,"alicsline":"","aliSimCount":10,"' \
                 'alihmmdesc":"","alisqdesc":"Chitinase 2","alisqname":"Q7DNA1","alisqfrom":288,"aliN":27}],"pdbs":["' \
                 '3iwr_A","3iwr_B","2dkv_A"],"evalue":"1.1e-09","nreported":1,"archindex":"55525311637573 39947"},{"a' \
                 'rchScore":"4","ph":"Streptophyta","arch":"PF00187.22 PF00182.22","kg":"Eukaryota","ndom":3,"extlink' \
                 '":"https://alphafold.ebi.ac.uk/entry/Q9SDY6","acc2":"Q9SDY6","taxid":"3847","acc":"Q9SDY6","taxlink' \
                 '":"http://www.uniprot.org/taxonomy/","desc":"Chitinase","pvalue":-32.2919606696729,"flags":"3","nre' \
                 'gions":3,"niseqs":16,"name":"Q9SDY6","species":"Glycine max","score":"41.4","bias":"10.5","sindex":' \
                 '"243314487","nincluded":"1","domains":[{"aliniseqs":1,"alisqacc":"F1","aliIdCount":24,"alirfline":"' \
                 '","is_included":1,"alihmmname":">Seq","bitscore":41.4002456665039,"ievalue":"3.4e-09","alisqto":65,' \
                 '"aliSim":0.731707317073171,"jali":65,"bias":"10.51","ienv":19,"cevalue":"3.4e-13","alimline":"+cg+q' \
                 'agg lcpn lccs++g+cg    +cg gcqs   s  ","alihmmfrom":88,"aliL":320,"is_reported":1,"alintseq":"","al' \
                 'isindex":"243314487","jenv":87,"alimmline":"","alihmmacc":"","oasc":"0.75","aliaseq":"QCGTQAGGALCPN' \
                 'RLCCSKFGWCGDTDSYCGEGCQSQCKS-A","alihmmto":129,"aliId":0.585365853658537,"alippline":"79999999999999' \
                 '99999999999999999999985322.2","alimodel":"kcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstd","aliM":171,"i' \
                 'ali":25,"alicsline":"","aliSimCount":30,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"Q9SDY6' \
                 '","alisqfrom":25,"aliN":42},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":10,"alirfline":"","is_inclu' \
                 'ded":0,"alihmmname":">Seq","bitscore":-2.17198753356934,"ievalue":"73000","alisqto":185,"aliSim":0.' \
                 '555555555555556,"jali":185,"bias":"5.35","ienv":121,"cevalue":"7.3","alimline":" a             +gyc' \
                 ' ++ +  +  c gg  pc a ","alihmmfrom":49,"aliL":320,"is_reported":0,"alintseq":"","alisindex":"243314' \
                 '487","jenv":211,"alimmline":"","alihmmacc":"","oasc":"0.61","aliaseq":"YA-------------WGYCFINEQNQAT' \
                 'YCDGGnwPCAAG","alihmmto":86,"aliId":0.37037037037037,"alippline":"33.............344444444444444443' \
                 '3334333","alimodel":"qaggatctnnqccsqygycgfgaeycgagcqgg..pcrad","aliM":171,"iali":159,"alicsline":""' \
                 ',"aliSimCount":15,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"Q9SDY6","alisqfrom":159,"ali' \
                 'N":40},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":3,"alirfline":"","is_included":0,"alihmmname":">' \
                 'Seq","bitscore":-2.33587241172791,"ievalue":"82000","alisqto":284,"aliSim":0.5,"jali":284,"bias":"0' \
                 '.21","ienv":251,"cevalue":"8.2","alimline":"gg  c + ","alihmmfrom":51,"aliL":320,"is_reported":0,"a' \
                 'lintseq":"","alisindex":"243314487","jenv":315,"alimmline":"","alihmmacc":"","oasc":"0.61","aliaseq' \
                 '":"GGLECGHG","alihmmto":58,"aliId":0.375,"alippline":"33333333","alimodel":"ggatctnn","aliM":171,"i' \
                 'ali":277,"alicsline":"","aliSimCount":4,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"Q9SDY6' \
                 '","alisqfrom":277,"aliN":8}],"evalue":"3.4e-09","nreported":1,"archindex":"55525311637573 3847"},{"' \
                 'archScore":"3","ph":"Streptophyta","arch":"PF00187.22 PF00182.22 PF00182.22","kg":"Eukaryota","ndom' \
                 '":2,"extlink":"https://alphafold.ebi.ac.uk/entry/C6T7J9","acc2":"C6T7J9","taxid":"3847","acc":"C6T7' \
                 'J9","taxlink":"http://www.uniprot.org/taxonomy/","desc":"Chitinase","pvalue":-31.1536505620012,"fla' \
                 'gs":"3","nregions":2,"niseqs":16,"name":"C6T7J9","species":"Glycine max","score":"39.8","bias":"8.0' \
                 '","sindex":"243339183","nincluded":"1","domains":[{"aliniseqs":1,"alisqacc":"F1","aliIdCount":21,"a' \
                 'lirfline":"","is_included":1,"alihmmname":">Seq","bitscore":39.7844619750977,"ievalue":"1.1e-08","a' \
                 'lisqto":71,"aliSim":0.684210526315789,"jali":71,"bias":"7.98","ienv":29,"cevalue":"1.1e-12","alimli' \
                 'ne":"+ n  c   lccs+ygycg g dycgkgc+ g c+ + ","alihmmfrom":7,"aliL":280,"is_reported":1,"alintseq":"' \
                 '","alisindex":"243339183","jenv":88,"alimmline":"","alihmmacc":"","oasc":"0.87","aliaseq":"AQNCGCEA' \
                 'ELCCSKYGYCGSGDDYCGKGCKEGPCYGTA","alihmmto":44,"aliId":0.552631578947368,"alippline":"57899*********' \
                 '********************9764","alimodel":"gsnmecpnnlccsqygycgmggdycgkgcqngacwtsk","aliM":171,"iali":34,' \
                 '"alicsline":"","aliSimCount":26,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"C6T7J9","alisq' \
                 'from":34,"aliN":38},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":5,"alirfline":"","is_included":0,"a' \
                 'lihmmname":">Seq","bitscore":-1.96306431293488,"ievalue":"63000","alisqto":172,"aliSim":0.5,"jali":' \
                 '172,"bias":"0.27","ienv":146,"cevalue":"6.3","alimline":"dyc k  ++  c  ","alihmmfrom":29,"aliL":280' \
                 ',"is_reported":0,"alintseq":"","alisindex":"243339183","jenv":194,"alimmline":"","alihmmacc":"","oa' \
                 'sc":"0.59","aliaseq":"DYCDKTNRHYPCAH","alihmmto":42,"aliId":0.357142857142857,"alippline":"44444333' \
                 '333322","alimodel":"dycgkgcqngacwt","aliM":171,"iali":159,"alicsline":"","aliSimCount":7,"alihmmdes' \
                 'c":"","alisqdesc":"Chitinase","alisqname":"C6T7J9","alisqfrom":159,"aliN":14}],"evalue":"1.1e-08","' \
                 'nreported":1,"archindex":"113074525985478 3847"},{"archScore":"3","ph":"Streptophyta","arch":"PF001' \
                 '87.22 PF00182.22","kg":"Eukaryota","ndom":2,"extlink":"https://alphafold.ebi.ac.uk/entry/I1MMY2","a' \
                 'cc2":"I1MMY2","taxid":"3847","acc":"I1MMY2","taxlink":"http://www.uniprot.org/taxonomy/","desc":"Ch' \
                 'itinase","pvalue":-30.5115182979563,"flags":"3","nregions":2,"niseqs":16,"name":"I1MMY2","species":' \
                 '"Glycine max","score":"38.9","bias":"10.9","sindex":"243358091","nincluded":"1","domains":[{"alinis' \
                 'eqs":1,"alisqacc":"F1","aliIdCount":24,"alirfline":"","is_included":1,"alihmmname":">Seq","bitscore' \
                 '":38.8729820251465,"ievalue":"2.0e-08","alisqto":75,"aliSim":0.557692307692308,"jali":75,"bias":"10' \
                 '.90","ienv":18,"cevalue":"2.0e-12","alimline":"  cg+q gg +cpn lccsq+g+cg     cg gcqs       p     +g' \
                 '","alihmmfrom":87,"aliL":317,"is_reported":1,"alintseq":"","alisindex":"243358091","jenv":89,"alimm' \
                 'line":"","alihmmacc":"","oasc":"0.76","aliaseq":"QNCGTQVGGVICPNGLCCSQYGWCGNTEAHCGRGCQSQCTPGSTPTPTTP' \
                 'SG","alihmmto":138,"aliId":0.461538461538462,"alippline":"47************************999*******98665' \
                 '55555443333","alimodel":"ikcgsqaggklcpnnlccsqwgfcglgsefcgggcqsgacstdkpcgkdagg","aliM":171,"iali":24' \
                 ',"alicsline":"","aliSimCount":29,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"I1MMY2","alis' \
                 'qfrom":24,"aliN":52},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":6,"alirfline":"","is_included":0,"' \
                 'alihmmname":">Seq","bitscore":0.743120551109314,"ievalue":"9400","alisqto":178,"aliSim":0.476190476' \
                 '190476,"jali":178,"bias":"0.96","ienv":117,"cevalue":"0.94","alimline":"wg+c ++    +  c sg   ","ali' \
                 'hmmfrom":107,"aliL":317,"is_reported":1,"alintseq":"","alisindex":"243358091","jenv":204,"alimmline' \
                 '":"","alihmmacc":"","oasc":"0.67","aliaseq":"WGYCFINERNQADYCTSGTRW","alihmmto":127,"aliId":0.285714' \
                 '285714286,"alippline":"555555555444555544322","alimodel":"wgfcglgsefcgggcqsgacs","aliM":171,"iali":' \
                 '158,"alicsline":"","aliSimCount":10,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"I1MMY2","a' \
                 'lisqfrom":158,"aliN":21}],"evalue":"2.0e-08","nreported":2,"archindex":"55525311637573 3847"},{"arc' \
                 'hScore":"2","ph":"Streptophyta","arch":"PF00187.22 PF00182.22 PF00182.22","kg":"Eukaryota","ndom":2' \
                 ',"extlink":"https://alphafold.ebi.ac.uk/entry/Q6K8R2","acc2":"Q6K8R2","taxid":"39947","acc":"Q6K8R2' \
                 '","taxlink":"http://www.uniprot.org/taxonomy/","desc":"Chitinase 6","pvalue":-30.1115021385755,"fla' \
                 'gs":"3","nregions":2,"niseqs":16,"name":"Q6K8R2","species":"Oryza sativa subsp. japonica","score":"' \
                 '38.3","bias":"8.5","sindex":"236253149","nincluded":"1","domains":[{"aliniseqs":1,"alisqacc":"F1","' \
                 'aliIdCount":18,"alirfline":"","is_included":1,"alihmmname":">Seq","bitscore":38.30517578125,"ievalu' \
                 'e":"3.0e-08","alisqto":60,"aliSim":0.771428571428571,"jali":60,"bias":"8.51","ienv":15,"cevalue":"3' \
                 '.0e-12","alimline":" +  c ++qccs++g+cg g++ycg gcq+gpc  ","alihmmfrom":51,"aliL":271,"is_reported":1' \
                 ',"alintseq":"","alisindex":"236253149","jenv":78,"alimmline":"","alihmmacc":"","oasc":"0.73","alias' \
                 'eq":"QSCGCASDQCCSKWGFCGTGSDYCGTGCQAGPCDV","alihmmto":85,"aliId":0.514285714285714,"alippline":"3455' \
                 '7888888888888888888888888888854","alimodel":"ggatctnnqccsqygycgfgaeycgagcqggpcra","aliM":171,"iali"' \
                 ':26,"alicsline":"","aliSimCount":27,"alihmmdesc":"","alisqdesc":"Chitinase 6","alisqname":"Q6K8R2",' \
                 '"alisqfrom":26,"aliN":35},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":10,"alirfline":"","is_include' \
                 'd":0,"alihmmname":">Seq","bitscore":-0.0346345528960228,"ievalue":"16000","alisqto":169,"aliSim":0.' \
                 '5,"jali":169,"bias":"0.48","ienv":140,"cevalue":"1.6","alimline":"nyc   s    c  g gy g g","alihmmfr' \
                 'om":144,"aliL":271,"is_reported":0,"alintseq":"","alisindex":"236253149","jenv":194,"alimmline":"",' \
                 '"alihmmacc":"","oasc":"0.73","aliaseq":"NYCdeTSTQWPCMAGKGYYGRG","alihmmto":163,"aliId":0.5,"alippli' \
                 'ne":"6663322333577788887776","alimodel":"nyc..cskwgscgigpgycgag","aliM":171,"iali":148,"alicsline":' \
                 '"","aliSimCount":10,"alihmmdesc":"","alisqdesc":"Chitinase 6","alisqname":"Q6K8R2","alisqfrom":148,' \
                 '"aliN":22}],"evalue":"3.0e-08","nreported":1,"archindex":"113074525985478 39947"},{"archScore":"4",' \
                 '"ph":"Streptophyta","arch":"PF00187.22 PF00182.22","kg":"Eukaryota","ndom":2,"extlink":"https://alp' \
                 'hafold.ebi.ac.uk/entry/B6TR38","acc2":"B6TR38","taxid":"4577","acc":"B6TR38","taxlink":"http://www.' \
                 'uniprot.org/taxonomy/","desc":"Chitinase","pvalue":-28.6942419861504,"flags":"3","nregions":2,"nise' \
                 'qs":16,"name":"B6TR38","species":"Zea mays","score":"36.3","bias":"15.4","sindex":"237740814","ninc' \
                 'luded":"1","domains":[{"aliniseqs":1,"alisqacc":"F1","aliIdCount":24,"alirfline":"","is_included":1' \
                 ',"alihmmname":">Seq","bitscore":36.2934341430664,"ievalue":"1.3e-07","alisqto":71,"aliSim":0.780487' \
                 '804878049,"jali":71,"bias":"15.39","ienv":23,"cevalue":"1.2e-11","alimline":" ++cg qaggatc +  ccs++' \
                 'g+cg  +eycgagcq+  c ","alihmmfrom":43,"aliL":334,"is_reported":1,"alintseq":"","alisindex":"2377408' \
                 '14","jenv":89,"alimmline":"","alihmmacc":"","oasc":"0.80","aliaseq":"GQQCGQQAGGATCRDCLCCSRFGFCGDTSE' \
                 'YCGAGCQS-QCT","alihmmto":84,"aliId":0.585365853658537,"alippline":"57899999999999999999999999999999' \
                 '999996.343","alimodel":"skrcgsqaggatctnnqccsqygycgfgaeycgagcqggpcr","aliM":171,"iali":31,"alicsline' \
                 '":"","aliSimCount":32,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname":"B6TR38","alisqfrom":31,"' \
                 'aliN":42},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":4,"alirfline":"","is_included":0,"alihmmname"' \
                 ':">Seq","bitscore":-3.81587338447571,"ievalue":"230000","alisqto":193,"aliSim":0.555555555555556,"j' \
                 'ali":193,"bias":"0.04","ienv":170,"cevalue":"23","alimline":"c+ g  y g","alihmmfrom":24,"aliL":334,' \
                 '"is_reported":0,"alintseq":"","alisindex":"237740814","jenv":201,"alimmline":"","alihmmacc":"","oas' \
                 'c":"0.49","aliaseq":"CAPGKKYFG","alihmmto":32,"aliId":0.444444444444444,"alippline":"222333333","al' \
                 'imodel":"cgmggdycg","aliM":171,"iali":185,"alicsline":"","aliSimCount":5,"alihmmdesc":"","alisqdesc' \
                 '":"Chitinase","alisqname":"B6TR38","alisqfrom":185,"aliN":9}],"evalue":"1.3e-07","nreported":1,"arc' \
                 'hindex":"55525311637573 4577"},{"archScore":"3","ph":"Streptophyta","arch":"PF00187.22 PF00182.22 P' \
                 'F00182.22","kg":"Eukaryota","ndom":2,"extlink":"https://alphafold.ebi.ac.uk/entry/B6TT00","acc2":"B' \
                 '6TT00","taxid":"4577","acc":"B6TT00","taxlink":"http://www.uniprot.org/taxonomy/","desc":"Chitinase' \
                 '","pvalue":-28.4890673188471,"flags":"3","nregions":2,"niseqs":16,"name":"B6TT00","species":"Zea ma' \
                 'ys","score":"36.0","bias":"10.1","sindex":"237844714","nincluded":"1","domains":[{"aliniseqs":1,"al' \
                 'isqacc":"F1","aliIdCount":19,"alirfline":"","is_included":1,"alihmmname":">Seq","bitscore":36.00219' \
                 '7265625,"ievalue":"1.5e-07","alisqto":59,"aliSim":0.735294117647059,"jali":59,"bias":"10.11","ienv"' \
                 ':19,"cevalue":"1.5e-11","alimline":" +  c +  ccs++gycg g +ycgagcq+gpc ","alihmmfrom":51,"aliL":271,' \
                 '"is_reported":1,"alintseq":"","alisindex":"237844714","jenv":74,"alimmline":"","alihmmacc":"","oasc' \
                 '":"0.64","aliaseq":"QNCGCASGLCCSRFGYCGTGEDYCGAGCQSGPCD","alihmmto":84,"aliId":0.558823529411765,"al' \
                 'ippline":"3445666666666666666666666666666664","alimodel":"ggatctnnqccsqygycgfgaeycgagcqggpcr","aliM' \
                 '":171,"iali":26,"alicsline":"","aliSimCount":25,"alihmmdesc":"","alisqdesc":"Chitinase","alisqname"' \
                 ':"B6TT00","alisqfrom":26,"aliN":34},{"aliniseqs":1,"alisqacc":"F1","aliIdCount":10,"alirfline":"","' \
                 'is_included":0,"alihmmname":">Seq","bitscore":0.288663238286972,"ievalue":"13000","alisqto":169,"al' \
                 'iSim":0.571428571428571,"jali":169,"bias":"0.35","ienv":99,"cevalue":"1.3","alimline":" nyc    ++w ' \
                 ' c  g gy g g","alihmmfrom":143,"aliL":271,"is_reported":0,"alintseq":"","alisindex":"237844714","je' \
                 'nv":192,"alimmline":"","alihmmacc":"","oasc":"0.81","aliaseq":"KNYCDrnnTQW-PCQAGKGYYGRG","alihmmto"' \
                 ':163,"aliId":0.476190476190476,"alippline":"45554211344.577777777766","alimodel":"nnycc...skwgscgig' \
                 'pgycgag","aliM":171,"iali":147,"alicsline":"","aliSimCount":12,"alihmmdesc":"","alisqdesc":"Chitina' \
                 'se","alisqname":"B6TT00","alisqfrom":147,"aliN":24}],"evalue":"1.5e-07","nreported":1,"archindex":"' \
                 '113074525985478 4577"}],"algo":"phmmer","stats":{"page":1,"nhits":"36","elapsed":"0.10","Z":362042,' \
                 '"Z_setby":0,"n_past_msv":20785,"unpacked":"36","user":0,"domZ_setby":0,"nseqs":362042,"n_past_bias"' \
                 ':5893,"sys":0,"n_past_fwd":38,"total":4,"nmodels":1,"nincluded":"36","n_past_vit":541,"nreported":"' \
                 '36","domZ":36},"uuid":"B401DCCA-8999-11EC-BB2E-100EE976C163","_internal":{"lowevalue":"1.5e-07","hi' \
                 'ghevalue":"5.4e-74"}}}'
