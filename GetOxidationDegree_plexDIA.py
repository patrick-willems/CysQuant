#!/usr/bin/env python

"""Script to calculate the cysteine oxidation degree from the plexDIA main output report.
This script is only sensible if light and heavy iodoacetamide labels were used within samples. 
Please see xxx for more information."""

__author__      = "Patrick Willems"
__email__       = "willems533@gmail.com"
__copyright__   = "Copyright 2023, Ghent University"

import re, sys, os, argparse
import statistics as stat
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
        protein = re.match(r'>\w+\|(\S+)\|\S+_',line).group(1)
        prot.update({protein: ''})
    else: prot[protein] += line
print('  Stored '+str(len(prot))+' proteins from '+args.fasta_file+'.')

#Store disulfides
print('#Reading disulfides')
ss = {}
txt = open('../Disulfides.txt','r')
for line in txt:
    line = line.rstrip()
    protein,site,res = line.split('\t')
    id = protein+'.'+site
    ss.update({id: res})
print('  Stored '+str(len(ss))+' AlphaFold2 predictions')

#Store SUBA
print('#Reading SUBA')
suba = {}
txt = open('../SUBA.txt','r')
for line in txt:
    line = line.rstrip()
    protein,loc = line.split('\t')
    suba.update({protein: loc})
print('  Stored '+str(len(suba))+' SUBA5 consensus locations')

print('#Reading main DIA-NN TSV report')
try: f = open(args.input_file, 'r')
except OSError: sys.exit("Could not open DIA-NN main report provided to script (-in): "+args.input_file)
#Get index of required columns from the header line - quite long but should be more robust for long term
ind = {}       
header = next(f)
headers = header.split('\t')
for i in range(len(headers)):
    if headers[i] == 'Run': ind.update({'run': i})                                 #Sample
    if headers[i] == 'Stripped.Sequence': ind.update({'seq': i})                   #Plain peptide
    if headers[i] == 'Precursor.Id': ind.update({'prec': i})                       #Peptide precursor
    if headers[i] == 'Ms1.Area': ind.update({'ms1': i})                #MS1 quantification Ms1.Area Precursor.Translated Precursor.Normalised
    if headers[i] == 'Precursor.Translated': ind.update({'ms2': i})                #MS2 quantification
    if headers[i] == 'Channel.Q.Value': ind.update({'ch_q': i})                    #Channel Q value
    if headers[i] == 'Translated.Q.Value': ind.update({'tr_q': i})                 #Translated Q value
    if headers[i] == 'Protein.Ids': ind.update({'prot': i})                        #Protein IDs
    if headers[i] == 'First.Protein.Description': ind.update({'descr': i})         #Description
    if headers[i] == 'Channel.Evidence.Ms1': ind.update({'evid_ms1': i})           #Ms1 channel evidence
    if headers[i] == 'Channel.Evidence.Ms2': ind.update({'evid_ms2': i})           #Ms2 channel evidence
    if headers[i] == 'Ms1.Profile.Corr': ind.update({'corr_ms1': i})               #Ms1 profile correlation
    

samples = {}    #Keep track of all samples in the report
data = {}       #Here all actual cysteine data will be stored
summaries = {}  #Oxidation statistics per sample
groupings = ['ctrl','10min','60min','180min']
        
for line in f:
    line = line.rstrip()
    
    #All required data is stored to variables based on the header indices parsed above
    col = line.split('\t')  
    run, seq, prec, ms1, ch_q, tr_q, protein_id, descr, evidence_ms1, evidence_ms2, corr_ms1 = col[ind['run']],col[ind['seq']],col[ind['prec']],float(col[ind['ms1']]),float(col[ind['ch_q']]),float(col[ind['tr_q']]),col[ind['prot']],col[ind['descr']],float(col[ind['evid_ms1']]),float(col[ind['evid_ms2']]),float(col[ind['corr_ms1']])
    samples.update({run: True}) #Store all samples
    sampleNr = int(re.search(r'_(\d+)_DIA',run).group(1))
    if sampleNr < 5: sample_group = "ctrl"
    elif sampleNr < 9: sample_group = "10min"
    elif sampleNr < 13: sample_group = "60min"
    elif sampleNr < 17: sample_group = "180min"
    if run not in summaries: summaries.update({run: {'nr': 0, 'ox': 0, 'sn_ox': 0,'sn_red': 0, 'int': []}})

    #Only consider peptides with a single Cysteine with MS1 level quantification
    if seq.count('C') == 1 and ms1 > 0: 
        if 'IAM-C-H' in prec: label = 'H'
        else: label = 'L'
        precursor = re.sub("IAM-C-\S","IAM",prec)
        if precursor not in data:
            data.update({precursor:{'prot': protein_id, 'descr': descr, 'seq': seq, 'samples': {}, 'groups': {}}})
            for id in samples: data[precursor]['samples'].update({id: {'tr_q': 0, 'ch_q': 0, 'evid_ms1': 0, 'evid_ms2': 0, 'corr_ms1': 0}})
            for id in groupings: data[precursor]['groups'].update({id: []})
        if tr_q < args.channel_qval: data[precursor]['samples'][run]['tr_q'] += 1   #
        if ch_q < args.translated_qval: data[precursor]['samples'][run]['ch_q'] += 1
        if evidence_ms1 > args.evidence: data[precursor]['samples'][run]['evid_ms1'] += 1
        if evidence_ms2 > args.evidence: data[precursor]['samples'][run]['evid_ms2'] += 1
        if corr_ms1 > args.correlation: data[precursor]['samples'][run]['corr_ms1'] += 1
        data[precursor]['samples'][run].update({label: ms1})   #Store the MS1 intensity for the correct L/H label per precursor
        #If meets the Channel Q and Translated Q thresholds and both L/H areas are stored - calculate the oxidation ratio
        if all(k in data[precursor]['samples'][run] for k in ('L','H')) and data[precursor]['samples'][run]['tr_q'] >= args.nr_translated_qval and data[precursor]['samples'][run]['ch_q'] >= args.nr_channel_qval and data[precursor]['samples'][run]['corr_ms1'] > 1: #data[precursor]['samples'][run]['corr_ms1'] > 1 (data[precursor]['samples'][run]['evid_ms1'] > 1 or data[precursor]['samples'][run]['evid_ms2'] > 1):
            if args.label_ox == 'L': ratio = data[precursor]['samples'][run]['L'] / (data[precursor]['samples'][run]['H'] + data[precursor]['samples'][run]['L'])
            else: ratio = ratio = data[precursor]['samples'][run]['H'] / (data[precursor]['samples'][run]['H'] + data[precursor]['samples'][run]['L'])
            data[precursor]['groups'][sample_group].append(ratio)
            data[precursor]['samples'][run].update({'OX': ratio})
            summaries[run]['ox'] += ratio
            summaries[run]['nr'] += 1
            sumInt = (data[precursor]['samples'][run]['H']+ data[precursor]['samples'][run]['L'])
            summaries[run]['int'].append(sumInt)
    else: continue

if args.singleton == 1:
    for precursor in data:
        for run in data[precursor]['samples']:
            sampleNr = int(re.search(r'_(\d+)_DIA',run).group(1))
            if sampleNr < 5: sample_group = "ctrl"
            elif sampleNr < 9: sample_group = "10min"
            elif sampleNr < 13: sample_group = "60min"
            elif sampleNr < 17: sample_group = "180min"
            if all(k in data[precursor]['samples'][run] for k in ('L','H')) or data[precursor]['samples'][run]['tr_q'] < 1 or data[precursor]['samples'][run]['ch_q'] < 1 or data[precursor]['samples'][run]['corr_ms1'] == 0 or data[precursor]['samples'][run]['evid_ms1'] == 0 or data[precursor]['samples'][run]['evid_ms2'] == 0: continue
            elif 'L' in data[precursor]['samples'][run] and data[precursor]['samples'][run]['L'] > stat.quantiles(summaries[run]['int'],n=4)[0]: #(stat.mean(summaries[run]['int'])*0.15):
                if args.label_ox == 'L':
                    data[precursor]['samples'][run].update({'OX': 1})
                    data[precursor]['groups'][sample_group].append(1)
                    summaries[run]['sn_ox'] += 1
                else:
                    data[precursor]['samples'][run].update({'OX': 0})
                    data[precursor]['groups'][sample_group].append(0)
                    summaries[run]['sn_red'] += 1
            elif 'H' in data[precursor]['samples'][run] and data[precursor]['samples'][run]['H'] > stat.quantiles(summaries[run]['int'],n=4)[0]: #(stat.mean(summaries[run]['int'])*0.15):
                if args.label_ox == 'H':
                    data[precursor]['samples'][run].update({'OX': 1})
                    data[precursor]['groups'][sample_group].append(1)
                    summaries[run]['sn_ox'] += 1
                else:
                    data[precursor]['samples'][run].update({'OX': 0})
                    data[precursor]['groups'][sample_group].append(0)
                    summaries[run]['sn_red'] += 1

#Print out summaries
print('#Writing out cysteine oxidation quantification summaries per sample:')
ave_ratios = 0
ave_ox = 0
ave_sn_red = 0
ave_sn_ox = 0
for run in summaries:
    average = 100*round((summaries[run]['ox']/summaries[run]['nr']),5)
    ave_ox += average
    ave_ratios += summaries[run]['nr']
    print('  '+run+': '+str(summaries[run]['nr'])+' ratios, '+str(average)[0:6]+'% Cys OX', end='')
    if args.singleton == 1:
        print(', '+str(summaries[run]['sn_ox'])+' singletons 100% OX, '+str(summaries[run]['sn_red'])+' singleton 0% OX', end='')
        ave_sn_red += summaries[run]['sn_red']
        ave_sn_ox += summaries[run]['sn_ox']
    print()
print('  --> All averaged: '+str((ave_ratios/len(summaries)))+' ratios, '+str((ave_ox/len(summaries)))[0:6]+'% Cys OX')
if args.singleton == 1: print('  --> All singletons: '+str((ave_sn_ox/len(summaries)))+' singletons 100% OX, '+str((ave_sn_red/len(summaries)))+' singleton 0% OX')

#Print out reports
outfile = 'ms1quant_' + args.input_file
print('#Writing out MS1 cysteine oxidation report to '+outfile)
out = open(outfile,'w')
#samples_s = {k: samples[k] for k in sorted(samples)}    #Sort the sample from 1 to 16
samples_s = ['B2432_Bp_IonOpt_CMB-1251_KRGEV_Willems_1_DIA.raw','B2450_Bp_IonOpt_CMB-1251_KRGEV_Willems_2_DIA.raw','B2466_Bp_IonOpt_CMB-1251_KRGEV_Willems_3_DIA.raw','B2482_Bp_IonOpt_CMB-1251_KRGEV_Willems_4_DIA.raw','B2436_Bp_IonOpt_CMB-1251_KRGEV_Willems_5_DIA.raw','B2454_Bp_IonOpt_CMB-1251_KRGEV_Willems_6_DIA.raw','B2470_Bp_IonOpt_CMB-1251_KRGEV_Willems_7_DIA.raw','B2486_Bp_IonOpt_CMB-1251_KRGEV_Willems_8_DIA.raw','B2440_Bp_IonOpt_CMB-1251_KRGEV_Willems_9_DIA.raw','B2458_Bp_IonOpt_CMB-1251_KRGEV_Willems_10_DIA.raw','B2474_Bp_IonOpt_CMB-1251_KRGEV_Willems_11_DIA.raw','B2490_Bp_IonOpt_CMB-1251_KRGEV_Willems_12_DIA.raw','B2444_Bp_IonOpt_CMB-1251_KRGEV_Willems_13_DIA.raw','B2462_Bp_IonOpt_CMB-1251_KRGEV_Willems_14_DIA.raw','B2478_Bp_IonOpt_CMB-1251_KRGEV_Willems_15_DIA.raw','B2494_Bp_IonOpt_CMB-1251_KRGEV_Willems_16_DIA.raw']
group_header = ''
for group in groupings: group_header += group + "_validRatio\t" + group + "_mean\t" + group + "_CV\t"
out.write("Peptide\tPrecursor\tProteins\tDescription_first\tProteins_position\tAlphaFold2\tSUBA5\t" + group_header + '\t'.join(samples_s) + '\n')

#Print quantifications
for precursor in data:
    seq = data[precursor]['seq']
    proteins = data[precursor]['prot'].split(';')

    #Determine cysteine position(s) in protein
    position,classes,loc = '','',''
    for protein in proteins:
        if protein not in prot: continue
        if protein in suba: loc += suba[protein] + ','
        posCys = str(seq.find('C') + prot[protein].find(seq) + 1)
        position += posCys + ';'
        id = protein + '.' + posCys
        if id in ss: classes += ss[id]
    if 'disulfide' in classes: classes = 'disulfide'
    elif 'metal' in classes: classes = 'metal'
    else: classes = 'none'
    locs = loc[:-1].split(',')
    loc = str(set(locs))

    #First general peptide precursor fields
    to_write = seq+'\t'+precursor+'\t'+data[precursor]['prot']+'\t'+data[precursor]['descr']+'\t'+position[:-1]+'\t'+classes+'\t'+loc
    
    #Per group MS1 quantifiation
    for group in groupings: 
        if len(data[precursor]['groups'][group]) > 0:
            ratios = data[precursor]['groups'][group]
            valid = len(ratios)
            mean = stat.mean(ratios)
            if valid > 1:
                stdev = stat.stdev(ratios)
                if mean == 0: cv = 0
                else: cv = format((stdev/mean)*100,'.3f')
            else: (stdev,cv) = ('NA','NA')
            to_write += '\t' + str(valid) + '\t' + str(mean) + '\t' + str(cv)
        else: to_write += '\t0\tNA\tNA'

    #Per sample MS1 quantification
    valid_values = 0
    for sample in samples_s: 
        if 'OX' in data[precursor]['samples'][sample]:
            to_write += '\t' + str(data[precursor]['samples'][sample]['OX'])
            valid_values += 1
        else: to_write += '\tNA'
    
    if valid_values > 0: out.write(to_write+'\n')