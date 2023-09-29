# CysQuant

"CysQuant: Simultaneous quantification of cysteine oxidation and protein abundance using data dependent or independent acquisition mass spectrometry" by Jingjing Huang, An Staes, Francis Impens, Vadim Demichev, Frank Van Breusegem, Kris Gevaert and Patrick Willems. 2023 _Redox Biology_ DOI: https://doi.org/10.1016/j.redox.2023.102908

This repository provides a Python script to convert DIA-NN/plexDIA reports to a matrix of recorded cysteine oxidation degrees (%) in each sample (columns) per cysteine peptide precursor (rows).
This script is only sensible when using isotopologous IAM labeling followed by plexDIA analysis with correct parameters (see paper for detailed info).
See below for more information for each input argument and info to perform a test run.


## Input arguments
**Required arguments**

 ```-in INPUT_FILE, --input_file INPUT_FILE```
 
 DIA-NN/plexDIA main output report (default: report.tsv)
 
 ```-fasta FASTA_FILE, --fasta_file FASTA_FILE```
 UniProtKB FASTA file used in the search. Should also work with non-UniProtKB databases but not recommended.
  
 -ox {L,H}, --label_ox {L,H}
 Specify the label (either 'L' or 'H' for light/heavy) that was used to label oxidized cysteines (i.e. after DTT/TCEP reduction).
  
**Optional arguments**
The optional arguments below are there to set certain quality thresholds on the quantifications to be included in the final report. These are the provided plexDIA Q-values, as well as evidence and correlation scores. In addition, light/heavy only quantified precursors can be included. These are required to have an intensity above the 25% quantile of all quantified L/H Cys precursors in the run (besides other specified quality threshold as for L/H pairs). 

  -chq CHANNEL_QVAL, --channel_qval CHANNEL_QVAL
  Threshold of Channel Q value applied. Default is 0.01.
  
  -trq TRANSLATED_QVAL, --translated_qval TRANSLATED_QVAL
  Threshold of Translated Q value applied. Default is 0.01.
  
  -nr_chq {0,1,2}, --nr_channel_qval {0,1,2}
  Number of Channel Q values required under threshold per precursor (0, 1 [L or H], or 2 [L and H]). Default is 1.
  
  -nr_trq {0,1,2}, --nr_translated_qval {0,1,2}
  Number of Translated Q values required under threshold per precursor (0, 1 [L or H], or 2 [L and H]). Default is 1.
  
  -single {0,1}, --singleton {0,1}
  Include peptide precursors only identified as either L or H, interpreted as being 0 or 100 percent oxidized dependent on the labeling. Default is not to include (0).
  
  -ev EVIDENCE, --evidence EVIDENCE
  Channel.Evidence.Ms1/2 filter, a value above this threshold is required in both L and H precursor. Default is -1 (so no filter).
  
  -corr CORRELATION, --correlation CORRELATION
  Ms1.Profile.Correlation, a value above this threshold is required in both L and H precursor. Default is -1 (so no filter).

  **Test data**
  Run as follows a test sample:

  ```python3 CysQuant.py -in ./test_data/report.tsv -fasta ./test_data/UP000006548_3702.fasta -ox H -single 1```

  This should return an output file report_CysQuant.tsv with the resulting cysteine oxidation degree matrix.
