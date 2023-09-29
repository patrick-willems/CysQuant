#!/usr/bin/env python

"""Script to calculate the cysteine oxidation degree from the plexDIA main output report.
This script is only sensible if light and heavy iodoacetamide labels were used within samples. 
Please see 'CysQuant: Simultaneous quantification of cysteine oxidation and protein abundance using data dependent or independent acquisition mass spectrometry' (DOI: 10.1016/j.redox.2023.102908) for more information."""

__author__      = "Patrick Willems"
__email__       = "willems533@gmail.com"
__copyright__   = "Copyright 2023, Ghent University"

import re, sys, os, argparse, warnings
import statistics as stat

#Argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input_file", required=True, help="Output file of labeled DIA-NN run (Cys IAM labels have to be set!), default name of report is 'report.tsv'")
parser.add_argument("-fasta", "--fasta_file", required=True, help="UniProtKB FASTA file used in the search.")
parser.add_argument("-ox", "--label_ox", required=True, choices=['L', 'H'], help="Specify the label (either 'L' or 'H' for light/heavy) that was used to label oxidized cysteines (i.e. after DTT/TCEP reduction) .")
parser.add_argument("-chq", "--channel_qval", default=0.01, type=float, help="Threshold of Channel Q value applied. Default is 0.01.")
parser.add_argument("-trq", "--translated_qval", default=0.01, type=float, help="Threshold of Translated Q value applied. Default is 0.01.")
parser.add_argument("-nr_chq", "--nr_channel_qval", default=1, choices=range(0,3), type=int, help="Number of Channel Q values required under threshold per precursor (0, 1 [L/H], or 2 [L+H]). Default is 1.")
parser.add_argument("-nr_trq", "--nr_translated_qval", default=1, choices=range(0,3), type=int, help="Number of Translated Q values required under threshold per precursor (0, 1 [L/H], or 2 [L+H]). Default is 1.")
parser.add_argument("-single", "--singleton", default=0, choices=range(0,2), type=int, help="Include peptide precursors only identified as either L or H, indicated then as being 0 or 100 percent oxidized dependent on the labeling. Default is not to include (0).")
parser.add_argument("-ev", "--evidence", default=-1, type=float, help="Channel.Evidence.Ms1/2 filter, a value above this threshold is required in both L and H precursor. Default is -1 (so no filter).")
parser.add_argument("-corr", "--correlation", default=-1, type=float, help="Ms1.Profile.Correlation, a value above this threshold is required in both L and H precursor. Default is -1 (so no filter).")
args = parser.parse_args()

#Store UniProt protein sequences
print('#Reading FASTA')
prot = {}
try: fasta = open(args.fasta_file, 'r')
except OSError: sys.exit("Could not open/read fasta file provided to script (-fasta): "+args.fasta_file)
for line in fasta:
    line = line.rstrip()
    if '>' in line:
        try: protein = re.match(r'>\w+\|(\S+)\|\S+_',line).group(1)             #UniProtKB entry (recommended!)
        except AttributeError: protein = re.match(r'>(\S+)',line).group(1)      #Non-UniProtKB
        prot.update({protein: ''})
    else: prot[protein] += line
print('  Stored '+str(len(prot))+' proteins from '+args.fasta_file+'.')

#Read the main plexDIA/DIA-NN output report (default: report.tsv)
print('#Reading main DIA-NN TSV report')
try: f = open(args.input_file, 'r')
except OSError: sys.exit("Could not open DIA-NN main report provided to script (-in): "+args.input_file)

#Get index of required columns from the header line, lengthy but robust
ind = {}       
header = next(f)
headers = header.split('\t')
for i in range(len(headers)):
    if headers[i] == 'Run': ind.update({'run': i})                                 #Sample
    if headers[i] == 'Stripped.Sequence': ind.update({'seq': i})                   #Plain peptide
    if headers[i] == 'Precursor.Id': ind.update({'prec': i})                       #Peptide precursor
    if headers[i] == 'Ms1.Area': ind.update({'ms1': i})                            #MS1 quantification
    if headers[i] == 'Precursor.Translated': ind.update({'ms2': i})                #MS2 quantification
    if headers[i] == 'Channel.Q.Value': ind.update({'ch_q': i})                    #Channel Q value
    if headers[i] == 'Translated.Q.Value': ind.update({'tr_q': i})                 #Translated Q value
    if headers[i] == 'Protein.Ids': ind.update({'prot': i})                        #Protein IDs
    if headers[i] == 'First.Protein.Description': ind.update({'descr': i})         #Description
    if headers[i] == 'Channel.Evidence.Ms1': ind.update({'evid_ms1': i})           #Ms1 channel evidence
    if headers[i] == 'Channel.Evidence.Ms2': ind.update({'evid_ms2': i})           #Ms2 channel evidence
    if headers[i] == 'Ms1.Profile.Corr': ind.update({'corr_ms1': i})               #Ms1 profile correlation
    

data = {}       #Here all actual cysteine data will be stored
samples = {}  #Oxidation statistics per sample
line_nr = 0        
for line in f:
    line_nr += 1
    if line_nr % 100000 == 0: print('  Read '+str(line_nr)+ ' lines..')
    
    #All required data is stored to variables based on the header line parsed above
    line = line.rstrip()
    col = line.split('\t')  
    run, seq, prec, ms1, ch_q, tr_q, protein_id, descr, evidence_ms1, evidence_ms2, corr_ms1 = col[ind['run']],col[ind['seq']],col[ind['prec']],float(col[ind['ms1']]),float(col[ind['ch_q']]),float(col[ind['tr_q']]),col[ind['prot']],col[ind['descr']],float(col[ind['evid_ms1']]),float(col[ind['evid_ms2']]),float(col[ind['corr_ms1']])
    if run not in samples: samples.update({run: {'ox': [], 'sn_ox': 0,'sn_red': 0, 'int': []}})    #Store all samples

    #Only store single cysteine peptides with quantified MS1 area
    if seq.count('C') == 1 and ms1 > 0:
        if 'IAM-C-H' in prec: label = 'H'
        else: label = 'L'
        precursor = re.sub("IAM-C-\S","IAM",prec)
        if precursor not in data:
            data.update({precursor:{'prot': protein_id, 'descr': descr, 'seq': seq, 'samples': {}}})
            for sample in samples: data[precursor]['samples'].update({sample: {'tr_q': 0, 'ch_q': 0, 'evid_ms1': 0, 'evid_ms2': 0, 'corr_ms1': 0}}) 

        #Store Q value, evidence and correlation values - these can be filtered upon for the final report (see arguments)
        if tr_q < args.channel_qval: data[precursor]['samples'][run]['tr_q'] += 1  
        if ch_q < args.translated_qval: data[precursor]['samples'][run]['ch_q'] += 1
        if evidence_ms1 > args.evidence: data[precursor]['samples'][run]['evid_ms1'] += 1
        if evidence_ms2 > args.evidence: data[precursor]['samples'][run]['evid_ms2'] += 1
        if corr_ms1 > args.correlation: data[precursor]['samples'][run]['corr_ms1'] += 1

        #Store the MS1 intensity for the L/H precursor
        data[precursor]['samples'][run].update({label: ms1})

        #If L/H MS1 areas are stored and all the specified score thresholds are satisfied: calculate the oxidation degree
        if all(k in data[precursor]['samples'][run] for k in ('L','H')) and data[precursor]['samples'][run]['tr_q'] >= args.nr_translated_qval and data[precursor]['samples'][run]['ch_q'] >= args.nr_channel_qval and data[precursor]['samples'][run]['corr_ms1'] > 1: 
            if args.label_ox == 'L':    #Light IAM = oxidized Cys
                ratio = data[precursor]['samples'][run]['L'] / (data[precursor]['samples'][run]['H'] + data[precursor]['samples'][run]['L'])
            else:                       #Heavy IAM = oxidized Cys
                ratio = ratio = data[precursor]['samples'][run]['H'] / (data[precursor]['samples'][run]['H'] + data[precursor]['samples'][run]['L'])
            data[precursor]['samples'][run].update({'OX': ratio})
            samples[run]['ox'].append(ratio)
            samples[run]['int'].append(data[precursor]['samples'][run]['H'] + data[precursor]['samples'][run]['L'])

#Print out recorded ratios per sample
print('  Recorded cysteine oxidation degree per sample:')
for run in samples:
    average = 100*round((sum(samples[run]['ox'])/len(samples[run]['ox'])),5)
    print('  -> '+run+': '+str(len(samples[run]['ox']))+' Cys quantified, '+str(average)[0:6]+' average % Cys OX')

#Extracting L/H only cys precursors if matching tresholds and above 25% quantile of summed intensities of quantified precursors recorded above
if args.singleton == 1:
    print('#Extracting Cys only identified in light / heavy channels')
    for precursor in data:
        for run in data[precursor]['samples']:
            #Need to make sure the precursor is not present in both L/H channels and that it matches the specified score thresholds
            if all(k in data[precursor]['samples'][run] for k in ('L','H')) or data[precursor]['samples'][run]['tr_q'] < 1 or data[precursor]['samples'][run]['ch_q'] < 1 or data[precursor]['samples'][run]['corr_ms1'] == 0 or data[precursor]['samples'][run]['evid_ms1'] == 0 or data[precursor]['samples'][run]['evid_ms2'] == 0: continue
            elif 'L' in data[precursor]['samples'][run] and data[precursor]['samples'][run]['L'] > stat.quantiles(samples[run]['int'],n=4)[0]:
                if args.label_ox == 'L':
                    data[precursor]['samples'][run].update({'OX': 1})
                    samples[run]['sn_ox'] += 1
                else:
                    data[precursor]['samples'][run].update({'OX': 0})
                    samples[run]['sn_red'] += 1
            elif 'H' in data[precursor]['samples'][run] and data[precursor]['samples'][run]['H'] > stat.quantiles(samples[run]['int'],n=4)[0]:
                if args.label_ox == 'H':
                    data[precursor]['samples'][run].update({'OX': 1})
                    samples[run]['sn_ox'] += 1
                else:
                    data[precursor]['samples'][run].update({'OX': 0})
                    samples[run]['sn_red'] += 1

    #Print out samples
    print('  Recorded L/H only cysteine peptide precursors per sample:')
    for run in samples: print('  -> '+str(samples[run]['sn_ox'])+' Cys 100% OX, '+str(samples[run]['sn_red'])+' Cys 0% OX')

#Print out CysQuant report
outfile = args.input_file.replace('.tsv','_CysQuant.tsv')
print('#Writing out MS1 cysteine oxidation report to '+outfile)
out = open(outfile,'w')
out.write('Peptide\tPrecursor\tProteins\tDescription_first\tCysteine_pos\t' + '\t'.join(samples) + '\n')

#Print quantifications
for precursor in data:
    
    #First check if oxidation degrees are recorded
    valid_values = 0
    for sample in samples:
        if 'OX' in data[precursor]['samples'][sample]: valid_values += 1
    if valid_values == 0: continue  #None recorded, skip this peptide precursor

    #Plain peptide sequence, protein and cysteine position
    seq = data[precursor]['seq']
    proteins = data[precursor]['prot'].split(';')
    position = ''
    for protein in proteins:
        if protein not in prot:
            warnings.warn("The protein "+protein+" was identified in the search results but not found in the provided FASTA.")
            continue
        posCys = str(seq.find('C') + prot[protein].find(seq) + 1)
        position += posCys + ';'    #Trailing semicolon will be removed later

    #Write out all fields per peptide precursor
    out.write(seq+'\t'+precursor+'\t'+data[precursor]['prot']+'\t'+data[precursor]['descr']+'\t'+position[:-1])
    for sample in samples: 
        if 'OX' in data[precursor]['samples'][sample]: out.write('\t' + str(data[precursor]['samples'][sample]['OX']))
        else: out.write('\tNA')
    out.write('\n')
